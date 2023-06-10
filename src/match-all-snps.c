
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
	
	int debugger = 0;
	fprintf(stderr, "Debug: %i\n", debugger++);

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

	fprintf(stderr, "Debug: %i\n", debugger++);
	for (int i = 0; i<=MAXKMER; i++) {
		all_primers[i] = malloc(sizeof(kmer_bst_t));

		if (all_primers[i] == NULL) {
			fprintf(stderr, "Can't create entry %d\n", i);
			exit(1);
		}
		all_primers[i]->bigger = NULL;
		all_primers[i]->smaller = NULL;
		all_primers[i]->value = 0;
		all_primers[i]->id = "";
	}

	printf("All primers 0 is %ld\n", all_primers[0]->value);
	fprintf(stderr, "Debug: %i\n", debugger++);
	// read the primer file
	read_primers_create_snps(opt->primers, all_primers, opt->reverse, opt->verbose);
	fprintf(stderr, "Debug: %i\n", debugger++);
	
	fprintf(stderr, "%sWe have read the primers%s\n", GREEN, ENDC);
	for (int i = 0; i<=MAXKMER; i++)
		print_all_primers(all_primers[i], i);
	

	struct R1_read **reads;
	reads = malloc(sizeof(*reads) * opt->tablesize);

	if (reads == NULL) {
		fprintf(stderr, "%sERROR: We can not allocate memory for a table size of %d. Please try a smaller value for -t%s\n", RED, opt->tablesize, ENDC);
		exit(2);
	}

	// Step 1. Read the R1 file and find the matches to any primer
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
	else
		match_out = stdout;
	
	bool warning_printed = false;

	struct encoded_kmers {
		uint64_t enc;
		*encoded_kmers next;
	}

	// A place to store our encdoded kmers
	// 
	// This should be an array of uint64_t with only
	// those set that are needed.
	encoded_kmers *enc_kmers[MAXKMER+1];
	for (int i=0; i<=MAXKMER; i++)
		enc_kmers[i] = malloc(sizeof(char) * i);
	



	while ((l = kseq_read(seq)) >= 0) {
		//  encode the first kmers in the sequence
		for (int i = 0; i<=MAXKMER; i++) {
			if (all_primers[i]->value > 0) {
				fprintf(stderr, "Encoded a kmer for %d\n", i);
				enc_kmers[
		uint64_t enc = kmer_encoding(seq->seq.s, 0, r1kmer);
		counts.R1_seqs++;
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


		kmer_bst_t *ks = find_primer(enc, i7l_primers);
                if (ks) {
			fprintf(match_out, "R1\t%s\t%s\t0\t-%ld\n", ks->id, seq->name.s, strlen(seq->seq.s));
			counts.R1_found++;
			R1read->trim = 0;
			unsigned hashval = hash(R1read->id) % opt->tablesize;
			R1read->next = reads[hashval];
			reads[hashval] = R1read;
			continue;
		} 

		for (int i=1; i<seq->seq.l - r1kmer + 1; i++) {
			enc = next_kmer_encoding(seq->seq.s, i, r1kmer, enc);
			kmer_bst_t *ks = find_primer(enc, i7l_primers);
			if (ks) {
				fprintf(match_out, "R1\t%s\t%s\t%d\t-%ld\n", ks->id, seq->name.s, i, strlen(seq->seq.s)-i);
				counts.R1_found++;
				R1read->trim = i;
				break;
			}
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

	// we don't need the I7left primers any more!
	free(i7l_primers);

	// Step 2. Read the R2 file and find the locations of I5right
	// In this case, we need to rc the string so we can make all possible combinations of it
	int r2kmer = strlen(opt->I5right);
	char* i5right_rc = malloc(sizeof(char) * strlen(opt->I5right)+1);
	rc(i5right_rc, opt->I5right);
	if (r2kmer >= MAXKMER) {
		fprintf(stderr, "%sKmer is %d. It needs to be less than %d, so we trimmed it to fit!%s\n", YELLOW, r2kmer, MAXKMER, ENDC);
		r2kmer = MAXKMER;
		i5right_rc[r2kmer]='\0';
	}

	//uint64_t i5r = kmer_encoding(opt->I5right, 0, r2kmer);
	// i5r = reverse_complement(i5r, r2kmer);
	
	// malloc for possibilities for the bst
	kmer_bst_t *i5r_primers;
	i5r_primers = malloc(sizeof(*i5r_primers) * ((3*r2kmer)+1));
	i5r_primers->bigger = NULL;
	i5r_primers->smaller = NULL;
	i5r_primers->value = 0;
	i5r_primers->id = "";

	// i5r_primers = (kmer_bst_t *) malloc(sizeof(*i5r_primers) * ((4*r2kmer)+1));
	if (i5r_primers == NULL) {
		fprintf(stderr, "%sCANNOT malloc for primer bst%s\n", RED, ENDC);
		exit(3);
	}


	// Open R2 for reading
	gzFile fp2 = gzopen(opt->R2_file, "r");
	if (fp2 == NULL) {
		fprintf(stderr, "%sERROR: Can not open %s%s\n", RED, opt->R2_file, ENDC);
		exit(3);
	}
	seq = kseq_init(fp2);
	
	// open our log files
	FILE *adjust;
	if (opt->R2_matches)
		match_out = fopen(opt->R2_matches, "w");
	else
		match_out = stdout;
	
	if (opt->adjustments)
		adjust = fopen(opt->adjustments, "w");
	else
		adjust = stdout;

	fprintf(adjust, "R1/R2\tSeq ID\tFrom\tTo\n");

	while ((l = kseq_read(seq)) >= 0) {
		// find the location of i5r if it is in this sequence
		counts.R2_seqs++;
		uint64_t enc = kmer_encoding(seq->seq.s, 0, r2kmer);
		int trim = -1;
		kmer_bst_t *ks = find_primer(enc, i5r_primers);
		if (ks) {
			fprintf(match_out, "R2\t%s\t%s\t0\t-%ld\n", ks->id, seq->name.s, strlen(seq->seq.s));
			counts.R2_found++;
			trim = 0;
		} else {
			for (int i=1; i<seq->seq.l - r2kmer + 1; i++) {
				enc = next_kmer_encoding(seq->seq.s, i, r2kmer, enc);
				kmer_bst_t *ks = find_primer(enc, i5r_primers);
				if (ks) {
					fprintf(match_out, "R2\t%s\t%s\t%d\t-%ld\n", ks->id, seq->name.s, i, strlen(seq->seq.s)-i);
					counts.R2_found++;
					trim = i;
				}
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
					fprintf(adjust, "R2\t%s\t%d\t%d\n", seq->name.s, trim, R1->trim);
					trim = R1->trim;
					counts.R2_adjusted++;
					break;
				}
				if (R1->trim == -1 && trim > -1) {
					fprintf(adjust, "R1\t%s\t%d\t%d\n", R1->id, R1->trim, trim);
					R1->trim = trim;
					counts.R1_adjusted++;
					break;
				}

				fprintf(stderr, "%sWe want to trim starting at %d from R1 and %d from R2 in %s. We went with the shorter%s\n", BLUE, R1->trim, trim, seq->name.s, ENDC);
				if (trim < R1->trim) {
					fprintf(adjust, "R1\t%s\t%d\t%d\n", R1->id, R1->trim, trim);
					R1->trim = trim;
					counts.R1_adjusted++;
					break;
				}
				if (trim > R1->trim) {
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


	}

	if (opt->adjustments)
		fclose(adjust);
	if (opt->R2_matches)
		fclose(match_out);
	kseq_destroy(seq);
	gzclose(fp2);

	// we don't need the I5right primers any more!
	free(i5r_primers);


	fprintf(stderr, "Total sequences: R1 %d R2 %d\n", counts.R1_seqs, counts.R2_seqs);
	fprintf(stderr, "Primer found: R1 %d R2 %d\n", counts.R1_found, counts.R2_found);
	fprintf(stderr, "Same Offset: %d (includes no adapter)\n", counts.same);
	fprintf(stderr, "Adjusted offset: R1 %d R2 %d\n", counts.R1_adjusted, counts.R2_adjusted);
	fprintf(stderr, "Sequences trimmed: R1 %d R2 %d\n", counts.R1_trimmed, counts.R2_trimmed);

*/
}

