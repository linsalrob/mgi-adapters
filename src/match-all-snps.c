
/*
 * Two more additions to our search strategy. Now we:
 * 1. Take in all primers, and match maximially upto 31 bp. So we match on different length primers
 * 2. Include all possible primers, so we read those from a fasta file
 * 3. Include all SNPs at all positions.
 * 4. Match the R1 and R2 files
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>
#include "colours.h"
#include "compare-seqs.h"
#include "create-snps.h"
#include "hash.h"
#include "kseq.h"
#include "match-all-snps.h"
#include "read_primers.h"
#include "rob_dna.h"
#include "seqs_to_ints.h"

KSEQ_INIT(gzFile, gzread);


void search_all_pairwise_snps(struct options *opt) {
	/*
	 * opt contains our variables for this search:
	 * 	opt->primers = name of a file of primers
	 * 	opt->reverse = include the reverse complement of the primers
	 * 	opt->verbose = sometimes more output!
	 * 	opt->tablesize = the size of the table to store R1 reads. Should be resonably large to avoid O(n) behaviour
	 */	
	
	typedef struct COUNTS {
		int R1_seqs;
		int R2_seqs;
		int R1_found;
		int R2_found;
		int R1_adjusted;
		int R2_adjusted;
		int R1_trimmed;
		int R2_trimmed;
		int same;
	} COUNTS;

	COUNTS counts = {};

	// create an array of kmer_bsts and malloc them
	kmer_bst_t *all_primers[MAXKMER+1];

	for (int i = 0; i<=MAXKMER; i++) {
		all_primers[i] = malloc(sizeof(kmer_bst_t));

		if (all_primers[i] == NULL) {
			fprintf(stderr, "Can't malloc memory for all_primers %d\n", i);
			exit(1);
		}
		all_primers[i]->bigger = NULL;
		all_primers[i]->smaller = NULL;
		all_primers[i]->value = 0;
		all_primers[i]->id = "";
	}

	// read the primer file
	read_primers_create_snps(opt->primers, all_primers, opt->reverse, opt->verbose);
	
	if (opt->debug) {	
		fprintf(stderr, "%sWe have read the primers%s\n", GREEN, ENDC);
		for (int i = 0; i<=MAXKMER; i++)
			print_all_primers(all_primers[i], i);
	}

	// A place to store our encdoded kmers
	// 
	// First, we make an array with just the lengths of the kmers we need to encode
	// and then we encode those kmers. This reduces our search from ~30 to ~3
	// (depending on how many primer lengths we have).
	// We also do this longer primers to shorter, so that we initially trim off the
	// longest possible primers

	int kmer_lengths[MAXKMER];
	int unique_kmer_count = 0;
	for (int i=MAXKMER; i>=0; i--) 
		if (all_primers[i]->value > 0) 
			kmer_lengths[unique_kmer_count++] = i; 	// we need to remember this kmer length

	// now we just add those kmers to this array and then we test them each time
	uint64_t encoded_kmers[unique_kmer_count];

	struct R1_read **reads;
	reads = malloc(sizeof(*reads) * opt->tablesize);

	if (reads == NULL) {
		fprintf(stderr, "%sERROR: We can not allocate memory for a table size of %d. Please try a smaller value for -t%s\n", RED, opt->tablesize, ENDC);
		exit(2);
	}

	if( access( opt->R1_file, R_OK ) == -1 ) {
		// file doesn't exist
		fprintf(stderr, "%sERROR: The file %s can not be found. Please check the file path%s\n", RED, opt->R1_file, ENDC);
		return;
	}
	if( access( opt->R2_file, R_OK ) == -1 ) {
		// file doesn't exist
		fprintf(stderr, "%sERROR: The file %s can not be found. Please check the file path%s\n", RED, opt->R2_file, ENDC);
		return;
	}

	// Step 1. Read the R1 file and find the matches to any primer
	
	gzFile fp1 = gzopen(opt->R1_file, "r");
	if (fp1 == NULL) {
		fprintf(stderr, "%sERROR: Can not open %s%s\n", RED, opt->R1_file, ENDC);
		exit(3);
	}
	kseq_t *seq = kseq_init(fp1);
	int l;

	FILE *match_out;
	if (opt->R1_matches)
		match_out = fopen(opt->R1_matches, "w");
	
	bool warning_printed = false;

	while ((l = kseq_read(seq)) >= 0) {
		counts.R1_seqs++;
		
		// housekeeping warnings and definitions
		if (!warning_printed && has_n(seq->seq.s)) {
			fprintf(stderr, "%sWARNING: sequences have an N but we don't deal with them. They are encoded as A%s\n", BLUE, ENDC);
			warning_printed = true;
		}
		struct R1_read *R1read;
		R1read = (struct R1_read *) malloc(sizeof(*R1read));
		if (R1read == NULL) {
			fprintf(stderr, "Can't allocate memory for new ID pointer\n");
			return 0;
		}
		R1read->trim = -1;
		R1read->id = strdup(seq->name.s);
		R1read->next = NULL;

		bool read_matched = false;
		// end housekeeping warnings and definitions

		//  encode the first kmers in the sequence
		for (int i = 0; i<unique_kmer_count; i++)
			encoded_kmers[i] = kmer_encoding(seq->seq.s, 0, kmer_lengths[i]);

		// test for the first kmers in our data structure. We need to iterate the 
		// array of encoded primers and see if we find an encoding in all_primers[k]
		// for k being the length of the encoding
		for (int i=0; i<unique_kmer_count; i++) {
			uint64_t enc  = encoded_kmers[i];
			kmer_bst_t *ks = find_primer(enc, all_primers[kmer_lengths[i]]);
			if (ks) {
				if (opt->R1_matches)
					fprintf(match_out, "R1\t%s\t%s\t0\t-%ld\n", ks->id, seq->name.s, strlen(seq->seq.s));
				counts.R1_found++;
				R1read->trim = 0;
				unsigned hashval = hash(R1read->id) % opt->tablesize;
				R1read->next = reads[hashval];
				reads[hashval] = R1read;
				read_matched = true;
			} 
		}

		if (read_matched)
			continue; // no point continuing if there is an adapter match at position 0!

		for (int posn=1; posn<seq->seq.l - MAXKMER + 1; posn++) {
			for (int i=0; i<unique_kmer_count; i++) {
				// calculate the next encoding for this kmer length
				uint64_t enc  = next_kmer_encoding(seq->seq.s, posn, kmer_lengths[i], encoded_kmers[i]);
				encoded_kmers[i] = enc; // remember it for next time!
				kmer_bst_t *ks = find_primer(enc, all_primers[kmer_lengths[i]]);
				if (ks) {
					if (opt->R1_matches)
						fprintf(match_out, "R1\t%s\t%s\t%d\t-%ld\n", ks->id, seq->name.s, posn, strlen(seq->seq.s)-posn);
					counts.R1_found++;
					R1read->trim = posn;
					read_matched = true;
				}
				if (read_matched)
					break;
			}
			if (read_matched)
				break;
		}

		unsigned hashval = hash(R1read->id) % opt->tablesize;
		R1read->next = reads[hashval];
		reads[hashval] = R1read;
	}

	// I am going to reset kseq so we have to initiate it again later
	kseq_destroy(seq);
	gzclose(fp1);

	if (opt->R1_matches)
		fclose(match_out);


	// Step 2. Read the R2 file and find the locations of any of the primers.

	// Open R2 for reading
	gzFile fp2 = gzopen(opt->R2_file, "r");
	if (fp2 == NULL) {
		fprintf(stderr, "%sERROR: Can not open %s%s\n", RED, opt->R2_file, ENDC);
		exit(3);
	}
	seq = kseq_init(fp2);

	// if we want to write the files, we open a pipe
	// otherwise it is null. so we just need to check before writing
	// Open R2 for writing
	FILE *pipe = NULL;
	if (opt->R2_output) {
		char* pipe_file = malloc(sizeof(char) * (strlen(opt->R2_output) + 10));
		strcpy(pipe_file, "gzip - > ");
		strcat(pipe_file, opt->R2_output);
		pipe = popen(pipe_file, "w");
		free(pipe_file);
	}

	// open our log files
	FILE *adjust;
	if (opt->R2_matches)
		match_out = fopen(opt->R2_matches, "w");

	if (opt->adjustments) {
		adjust = fopen(opt->adjustments, "w");
		fprintf(adjust, "R1/R2\tSeq ID\tFrom\tTo\n");
	}

	while ((l = kseq_read(seq)) >= 0) {
		counts.R2_seqs++;
		bool read_matched = false;
		//  encode the first kmers in the sequence
		for (int i = 0; i<unique_kmer_count; i++)
			encoded_kmers[i] = kmer_encoding(seq->seq.s, 0, kmer_lengths[i]);
		
		int trim = -1;
		
		for (int i=0; i<unique_kmer_count; i++) {
			uint64_t enc  = encoded_kmers[i];
			kmer_bst_t *ks = find_primer(enc, all_primers[kmer_lengths[i]]);

			if (ks) {
				if (opt->R2_matches)
					fprintf(match_out, "R2\t%s\t%s\t0\t-%ld\n", ks->id, seq->name.s, strlen(seq->seq.s));
				counts.R2_found++;
				trim = 0;
				read_matched = true;
			} 
		}

		if (!read_matched) {
			for (int posn=1; posn<seq->seq.l - MAXKMER + 1; posn++) {
				for (int i=0; i<unique_kmer_count; i++) {
					// calculate the next encoding for this kmer length
					uint64_t enc  = next_kmer_encoding(seq->seq.s, posn, kmer_lengths[i], encoded_kmers[i]);
					encoded_kmers[i] = enc; // remember it for next time!
					kmer_bst_t *ks = find_primer(enc, all_primers[kmer_lengths[i]]);
					if (ks) {
						if (opt->R2_matches)
							fprintf(match_out, "R2\t%s\t%s\t%d\t-%ld\n", ks->id, seq->name.s, posn, strlen(seq->seq.s)-posn);
						counts.R2_found++;
						trim = posn;
						read_matched = true;
					}
					if (read_matched)
						break;
				}
				if (read_matched)
					break;
			}
		}
		// we either have a value or -1 for trim.
		// Now find the matching R1
		unsigned hashval = hash(seq->name.s) % opt->tablesize;
		struct R1_read *R1 = reads[hashval];
		bool matched = false;
		while (R1 != NULL) {
			if (strcmp(R1->id, seq->name.s) == 0) {
				// we found a match
				matched = true;
				if (trim == R1->trim) {
					// nothing to do, we can process both reads
					counts.same++;
					break;
				}
				if (trim == -1 && R1->trim > -1) {
					if (opt->adjustments)
						fprintf(adjust, "R2\t%s\t%d\t%d\n", seq->name.s, trim, R1->trim);
					trim = R1->trim;
					counts.R2_adjusted++;
					break;
				}
				if (R1->trim == -1 && trim > -1) {
					if (opt->adjustments)
						fprintf(adjust, "R1\t%s\t%d\t%d\n", R1->id, R1->trim, trim);
					R1->trim = trim;
					counts.R1_adjusted++;
					break;
				}

				if (opt->verbose)
					fprintf(stderr, "%sWe want to trim starting at %d from R1 and %d from R2 in %s. We went with the shorter%s\n", BLUE, R1->trim, trim, seq->name.s, ENDC);
				if (trim < R1->trim) {
					if (opt->adjustments)
						fprintf(adjust, "R1\t%s\t%d\t%d\n", R1->id, R1->trim, trim);
					R1->trim = trim;
					counts.R1_adjusted++;
					break;
				}
				if (trim > R1->trim) {
					if (opt->adjustments)
						fprintf(adjust, "R2\t%s\t%d\t%d\n", seq->name.s, trim, R1->trim);
					trim = R1->trim;
					counts.R2_adjusted++;
					break;
				}
			}
			R1 = R1->next;
		}
		if (!matched) {
			fprintf(stderr, "%s We did not find an R1 that matches %s%s\n", PINK, seq->name.s, ENDC);
		}
		if (trim > -1) {
			seq->seq.s[trim] = '\0';
			seq->qual.s[trim] = '\0';
			counts.R2_trimmed++;
		}
		if (pipe)
			fprintf(pipe, "@%s %s\n%s\n+\n%s\n", seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);


	}
	if (pipe)
		pclose(pipe);

	
	// do we need to write to R1
	if (opt->R1_output) {
		// Step 3. Reread R1 and write the left reads, trimming at (strcmp(id, seq->name.s) == 0) -> trim
		// We only need to do this if we are going to write to the file.

		fp1 = gzopen(opt->R1_file, "r");
		if (fp1 == NULL) {
			fprintf(stderr, "%sERROR: Can not open %s%s\n", RED, opt->R1_file, ENDC);
			exit(3);
		}
		seq = kseq_init(fp1);

		char* pipe_file = malloc(sizeof(char) * (strlen(opt->R1_output) + 10));
		strcpy(pipe_file, "gzip - > ");
		strcat(pipe_file, opt->R1_output);

		pipe = popen(pipe_file, "w");

		while ((l = kseq_read(seq)) >= 0) {
			unsigned hashval = hash(seq->name.s) % opt->tablesize;
			struct R1_read *R1 = reads[hashval];
			while (R1 != NULL) {
				if (strcmp(R1->id, seq->name.s) == 0) {
					if (R1->trim > -1) {
						seq->seq.s[R1->trim] = '\0';
						seq->qual.s[R1->trim] = '\0';
						counts.R1_trimmed++;
					}
				}
				R1 = R1->next;
			}
			fprintf(pipe, "@%s %s\n%s\n+\n%s\n", seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);
		}

		pclose(pipe);

	}


	if (opt->adjustments)
		fclose(adjust);
	if (opt->R2_matches)
		fclose(match_out);
	kseq_destroy(seq);
	gzclose(fp2);

	for (int i=0; i<=MAXKMER; i++)
		free(all_primers[i]);

	free(reads);

	fprintf(stderr, "Total sequences: R1 %d R2 %d\n", counts.R1_seqs, counts.R2_seqs);
	fprintf(stderr, "Primer found: R1 %d R2 %d\n", counts.R1_found, counts.R2_found);
	fprintf(stderr, "Same Offset: %d (includes no adapter)\n", counts.same);
	fprintf(stderr, "Adjusted offset: R1 %d R2 %d\n", counts.R1_adjusted, counts.R2_adjusted);
	fprintf(stderr, "Sequences trimmed: R1 %d R2 %d\n", counts.R1_trimmed, counts.R2_trimmed);

}

