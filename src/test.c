#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "colours.h"
#include "compare-seqs.h"
#include "print-sequences.h"
#include "rob_dna.h"
#include "seqs_to_ints.h"

void test_primers() {
	// populate the bst
	
	kmer_bst_t *kmers;
	kmers = (kmer_bst_t *) malloc(sizeof(*kmers));
	uint64_t is[] = {10, 12, 11, 9, 5, 8, 14, 13};
	// uint64_t is[] = {5, 8, 9, 10, 11, 12, 13, 14};
	for (int i=0; i<= 7; i++) {
		char s[100];
		sprintf(s, "NAME: %ld", is[i]);
		add_primer(is[i], s, kmers);
		printf("Added: %ld\n", is[i]);
	}

	uint64_t is2[] = {10, 12, 11, 9, 5, 8, 14, 13, 17, 2, 99};
	for (int i=0; i<= 10; i++) {
		kmer_bst_t* ks = find_primer(is2[i], kmers);
		if ( ks == NULL )
			printf("NOT found: %ld\n", is2[i]);
		else
			printf("Found: %ld:id '%s'\n", is2[i], ks->id);
	}
}

void test_comparisons() {
	
	char * seq = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
	// char * seq = "CTCTCTCTCTCTCT";
	for (int i =1; i<32; i++) {
		uint64_t enc = kmer_encoding(seq, 0, i);
		char* bin = int_to_binary(enc);
		printf("%d gives: %li : %s\n", i, kmer_encoding(seq, 0, i), bin);
	}

	printf("\n\n\nNew hashes\n\n");
	seq = "CCCGCCCCTCCactgCCCCAAAAATTTT";
	uint64_t enc = kmer_encoding(seq, 0, 3);
	printf("i: 0 str: %s enc: %ld\t%s\n", substr(seq, 0, 3), enc, int_to_binary(enc));
	for (int i = 1; i<=20; i++) {
		enc = next_kmer_encoding(seq, i, 3, enc);
		char* ss = substr(seq, i, 3);
		printf("i: %d str: %s enc: %ld\t%s\n", i, ss, enc, int_to_binary(enc));
	}

}



void test_reverse_complement() {
	char* seq1 = "ATGC";
	char* seq2 = "GCAT";
	uint64_t enc1 =  kmer_encoding(seq1, 0, 4);
	uint64_t enc2 =  kmer_encoding(seq2, 0, 4);
	if (reverse_complement(enc1, 4) == enc2) 
		printf("Success! We did a rc\n");
	else
		printf("Can't reverse complement %s\n", seq1);
	seq1 = "ATGCATCAGCTAGCATACGTACGTA";
	seq2 = "TACGTACGTATGCTAGCTGATGCAT";
	enc1 =  kmer_encoding(seq1, 0, 25);
	enc2 =  kmer_encoding(seq2, 0, 4);
	printf("Testing encoding rc:\n");
	if (reverse_complement(enc1, 4) == enc2) 
		printf("\tSuccess! We did a rc of the encoding\n");
	else
		printf("\tCan't reverse complement the encoding of %s\n", seq1);
	printf("Testing char* rc:\n");

	char* seq1rc = malloc(sizeof(char) * strlen(seq1)+1);
	rc(seq1rc, seq1);
	if (strcmp(seq1rc, seq2) == 0)
		printf("\tSuccess! We did a rc of the char*\n");
	else
		printf("\t%sCan't reverse complement the char* \n%s\n%s%s\n", RED, seq1, seq1rc, ENDC);

}

void test_has_n() {
	char* non = "ATGCATGCTACGATCGACT";
	char* withn = "AGCTAGCATNAGCTAGCAT";
	if (has_n(non)) 
		printf("WRONG: %s doesn't have an N but we think it does\n", non);
	else
		printf("Correct: %s does not have an N\n", non);
	if (has_n(withn)) 
		printf("Correct: %s DOES have an N\n", withn);
	else
		printf("WRONG: %s DOES have an N but we think it does not\n", withn);

}

void test_kmer_enc_dec() {
	char* seq = "TGAAATGCTACGATCGACT";
	uint64_t enc = kmer_encoding(seq, 0, 19);
	printf("for %s we got %s\n", seq, int_to_binary(enc));
	char* dseq = kmer_decoding(enc, 18);
	if (strcmp(seq, dseq) == 0)
		printf("%s%s and %s are the same%s\n", GREEN, seq, dseq, ENDC);
	else 
		printf("%sStarted with %s ended with %s and they are not the same!%s\n", RED, seq, dseq, ENDC);
}

void test_enc() {
	char *test = "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
	uint64_t enc = kmer_encoding(test, 0, strlen(test)-1);
	fprintf(stderr, "enc of %s is %ld\n", test, enc);
}


int main(int argc, char *argv[]) {
	// test_primers();
	// test_comparisons();
	// test_reverse_complement();
	// test_has_n();
	test_kmer_enc_dec();
	test_enc();
}

