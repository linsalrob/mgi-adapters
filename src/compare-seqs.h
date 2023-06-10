
#ifndef COMPARE_SEQS_H
#define COMPARE_SEQS_H

#include <stdbool.h>

/*
 * This is an unbalanced binary search tree, and so could devolve into O(n)
 * performance, however with random ints it should be ~O(log n)
 */
typedef struct kmer_bst {
    uint64_t value;
    char* id;
    struct kmer_bst *bigger;
    struct kmer_bst *smaller;
} kmer_bst_t;

/* 
 * Add a sequence encoding (a long long int) to a binary
 * search tree.
 */
void add_primer(uint64_t, char*, kmer_bst_t*);

/*
 * Find an encoding in a bst
 */
kmer_bst_t* find_primer(uint64_t, kmer_bst_t*);

/*
 * recursively print all primers
 */

void print_all_primers(kmer_bst_t*, int);

#endif
