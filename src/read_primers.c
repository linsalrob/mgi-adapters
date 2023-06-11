/*
 * Read the primer file
 *
 * TODO: How do we want to store multiple primers? We should probably use an array of pointers to kmer_bst_t's then
 * we can just store at each position that we need to search for i=0; i<100. We need to malloc this for 101*sizeof(kmer_bst_t)
 */



#include <float.h>
#include <getopt.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "colours.h"
#include "compare-seqs.h"
#include "create-snps.h"
#include "kseq.h"
#include "print-sequences.h"
#include "read_primers.h"
#include "rob_dna.h"
#include "seqs_to_ints.h"
#include "version.h"

KSEQ_INIT(gzFile, gzread);

void read_primers(char* primerfile, kmer_bst_t** all_primers, bool reverse, int verbose) {
	/*
	 * encode the primers in primerfile 
	 *
	 * We expect all_primers to be an array of 0-MAXKMER kmer_bst_t's
	 *
	 * Note: we need to take the correct substring of the sequence to reverse complement! We need the rightmose k-bases 
	 * otherwise we have an offset of length(string) - kmer to where the match should be.
	 *
	 */
	if( access( primerfile, R_OK ) == -1 ) {
		// file doesn't exist
		fprintf(stderr, "%sERROR: The file %s can not be found. Please check the file path%s\n", RED, primerfile, ENDC);
		return;
	}

	gzFile fp;
	kseq_t *seq;

	fp = gzopen(primerfile, "r");
	seq = kseq_init(fp);
	int l;
	while ((l = kseq_read(seq)) >= 0) {
		int kmer = strlen(seq->seq.s);
		if (seq->seq.l > MAXKMER) {
			fprintf(stderr, "%sWARNING: Length of %s is longer than our maximum (%ld bp), so we had to truncate it%s\n", RED, seq->name.s, seq->seq.l, ENDC);
			kmer = MAXKMER;
		}
		uint64_t enc = kmer_encoding(seq->seq.s, 0, kmer);
		
		add_primer(enc, seq->name.s, all_primers[kmer]);
		if (reverse) {
			uint64_t rcenc = kmer_encoding(seq->seq.s, strlen(seq->seq.s) - kmer, kmer);
			enc = reverse_complement(rcenc, kmer);
			char revname[strlen(seq->name.s)+3];
			strcpy(revname, seq->name.s);
			strcat(revname, " rc");
			if (verbose)
				fprintf(stderr, "%sAdded a rc primer: %s %s\n", GREEN, revname, ENDC);
			add_primer(enc, revname, all_primers[kmer]);
		}
			
		if (verbose)
			fprintf(stderr, "%sEncoding %s with length %ld using k-mer %d%s\n", GREEN, seq->seq.s, seq->seq.l, kmer, ENDC);
	}
	kseq_destroy(seq);
	gzclose(fp);
}

void read_primers_create_snps(char* primerfile, kmer_bst_t** all_primers, bool reverse, int verbose) {
	/*
	 * encode the primers in primerfile and create all snps for all primers.
	 *
	 * We expect all_primers to be an array of 0-MAXKMER kmer_bst_t's
	 *
	 * Note: we need to take the correct substring of the sequence to reverse complement! We need the rightmose k-bases 
	 * otherwise we have an offset of length(string) - kmer to where the match should be.
	 *
	 */
	if( access( primerfile, R_OK ) == -1 ) {
		// file doesn't exist
		fprintf(stderr, "%sERROR: The file %s can not be found. Please check the file path%s\n", RED, primerfile, ENDC);
		return;
	}

	gzFile fp;
	kseq_t *seq;

	if (verbose)
		fprintf(stderr, "%sREADING %s%s\n", PINK, primerfile, ENDC);

	fp = gzopen(primerfile, "r");
	seq = kseq_init(fp);
	int l;
	while ((l = kseq_read(seq)) >= 0) {
		int kmer = strlen(seq->seq.s);
		if (seq->seq.l > MAXKMER) {
			if (verbose)
				fprintf(stderr, "%sWARNING: Length of %s is longer than our maximum (%ld bp), so we had to truncate it%s\n", RED, seq->name.s, seq->seq.l, ENDC);
			kmer = MAXKMER;
		}

		create_all_snps(seq->seq.s, kmer, seq->name.s, all_primers[kmer]);
		
		if (reverse) {
			char revname[strlen(seq->name.s)+3];
			strcpy(revname, seq->name.s);
			strcat(revname, " rc");
			if (verbose)
				fprintf(stderr, "%sAdded a rc primer: %s %s\n", GREEN, revname, ENDC);
			
			char* rcseq = malloc(strlen(seq->seq.s) + 1);
			rc(rcseq, seq->seq.s);
			create_all_snps(rcseq, kmer, revname, all_primers[kmer]);

		}
			
		if (verbose)
			fprintf(stderr, "%sEncoding %s with length %ld using k-mer %d%s\n", GREEN, seq->seq.s, seq->seq.l, kmer, ENDC);
	}
	kseq_destroy(seq);
	gzclose(fp);
}
