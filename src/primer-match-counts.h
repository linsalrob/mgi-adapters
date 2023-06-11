
#ifndef PRIMER_MATCH_COUNTS_H
#define PRIMER_MATCH_COUNTS_H

/*
 * Some counts and information about primer matches. Currently this is done
 * as an O(n) search each time we find a match because (a) we don't have a 
 * lot of primers (by name), and (b) we don't have a lot of matches.
 *
 * id: primer name
 * count: occurrence
 * before: array of counts of A, C, G, T as the preceeding base
 * after:  array of counts of A, C, G, T as the following base
 * next_primer: the next one in the list
 */

typedef struct primer_counts {
	char* id;
	int count;
	int before[5];
	int after[5];
	struct primer_counts *next_primer;
} primer_counts_t;

/*
 * Add a primer match
 * id: primer id
 * before: the base preceeding (or NULL)
 * after: the base following (or NULL)
 */

void count_primer_occurrence(primer_counts_t *pc, char * id, char before, char after);

void print_primers(primer_counts_t *pc, int);
#endif
