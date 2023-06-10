/*
 * abstracting out the method to create snps because we use it in a couple of places
 */

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "create-snps.h"
#include "seqs_to_ints.h"



void create_all_snps(char *adapter, int kmer, char* seqid, kmer_bst_t *primers) {
	/*
	 * calculate a SNP at every position. This is O(3*k) where k is length of kmer
	 */

	uint64_t enc = kmer_encoding(adapter, 0, kmer);
	add_primer(enc, seqid, primers);

	char *base = "ACGT";
	for (int i = 0; i < kmer; i++) {
		char* snp = malloc(sizeof(char) * strlen(adapter)+1);
		strcpy(snp, adapter);
		for (int j = 0; j<4; j++) {
			if (adapter[i] == base[j])
				continue;
			snp[i] = base[j];
			char* name = malloc(sizeof(char) * (strlen(seqid) + 10));
			sprintf(name, "%s %d %c->%c", seqid, i, adapter[i], snp[i]);
			uint64_t encs = kmer_encoding(snp, 0, kmer);
			add_primer(encs, name, primers);
		}
	}
}




