#ifndef MATCH_PAIRED_SNPS_H  
#define MATCH_PAIRED_SNPS_H  

#include <stdbool.h>

/*
 * The options that we need
 */

struct options {
	char* R1_file;
	char* R2_file;
	char* R1_output;
	char* R2_output;
	char* R1_matches;
	char* R2_matches;
	char* adjustments;
	char* I7left;
	char* I7right;
	int tablesize;
	bool debug;
};

/*
 * R1_read is a struct with the R1 name, the id string for the sequence, and whether there was a match to I7left
 * next is a pointer to the next R1_read element in the hash.
 */
struct R1_read {
	int trim;
	char *id;
	struct R1_read *next;
};


// how long should our lines be. This is a 64k buffer
#define MAXLINELEN 65536


/* this version tries to get all single bp snps in the adapters */

void trim_pairwise_snps(struct options *opt);

/* search but don't trim */

void search_all_pairwise_snps(struct options *opt);
#endif
