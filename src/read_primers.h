
#ifndef READ_PRIMERS_H
#define READ_PRIMERS_H

// maximum kmer length
#define MAXKMER 31

// read the priemrs and populate the kmer_bst_t
void read_primers(char*, kmer_bst_t**, bool, int);

// read the priemrs and populate the kmer_bst_t with a snp in every position
void read_primers_create_snps(char*, kmer_bst_t**, bool, int);

#endif

