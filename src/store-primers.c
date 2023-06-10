#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include "compare-seqs.h"
#include "print-sequences.h"
#include "seqs_to_ints.h"


void add_primer(uint64_t encoding, char* primerid, kmer_bst_t* ks) {
	/*
	 * Add a new integer to the BST. Either at the current location, and
	 * then set bigger/smaller to be new empty structs
	 * else at bigger or smaller.
	 *
	 * We use bigger && smaller == NULL to indicate we are at a new
	 * kmer_bst_t since the value is also not set
	 */
	
	if (ks->bigger == NULL && ks->smaller == NULL) {
		ks->value = encoding;
		ks->id = malloc(strlen(primerid) + 1);
		strcpy(ks->id, primerid);
		ks->bigger = (kmer_bst_t *) malloc(((2 * sizeof(*ks->bigger)) + sizeof(uint64_t) + sizeof(char *)));
		ks->smaller = (kmer_bst_t *)  malloc(((2 * sizeof(*ks->bigger)) + sizeof(uint64_t) + sizeof(char *)));
		ks->bigger->bigger = NULL;
		ks->bigger->smaller = NULL;
		ks->smaller->bigger = NULL;
		ks->smaller->smaller= NULL;
		ks->bigger->value = 0;
		ks->smaller->value = 0;
		ks->bigger->id = "";
		ks->smaller->id = "";
		return;
	}
	if (encoding > ks->value)
		return add_primer(encoding, primerid, ks->bigger);
	
	return add_primer(encoding, primerid, ks->smaller);
}

kmer_bst_t* find_primer(uint64_t encoding, kmer_bst_t* ks) {
	/*
	 * Find the int in the bst
	 * Note, as above we use bigger && smaller == NULL to define
	 * not set because value = 0 could be AAAAAAA
	 */

	if (ks->value == encoding)
		return ks;
	if (ks->bigger == NULL && ks->smaller == NULL) 
		return NULL;
	if (encoding > ks->value)
		return find_primer(encoding, ks->bigger);
	return find_primer(encoding, ks->smaller);
}


void print_all_primers(kmer_bst_t* ks, int kmer) {
	/*
	 * print all the primers in this data structure
	 * it does this recursively
	 *
	 * ks is our data structure, kmer is the length of the kmers
	 */


	if (ks->bigger)
		print_all_primers(ks->bigger, kmer);

	if (ks->value > 0)
		printf("%s: %s (%ld)\n", ks->id, kmer_decoding(ks->value, kmer), ks->value);
	
	if (ks->smaller)
		print_all_primers(ks->smaller, kmer);
}



