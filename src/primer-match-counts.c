#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "primer-match-counts.h"


void count_primer_occurrence(primer_counts_t *pc, char * id, char before, char after) {
	// do we already have id
	bool found = false;

	primer_counts_t *n = pc;

	while (n->id != NULL && !found) {
		if (strcmp(n->id, id) == 0) {
			n->count++;
			found = true;
			switch(before) {
				case 'A':
					n->before[0]++;
					break;
				case 'C':
					n->before[1]++;
					break;
				case 'G':
					n->before[2]++;
					break;
				case 'T':
					n->before[3]++;
					break;
				default:
					n->before[4]++;
			}

			switch(after) {
				case 'A':
					n->after[0]++;
					break;
				case 'C':
					n->after[1]++;
					break;
				case 'G':
					n->after[2]++;
					break;
				case 'T':
					n->after[3]++;
					break;
				default:
					n->after[4]++;
			}
		}
		n = n->next_primer;
	}

	if (found)
		return;
	n->id = id;
	n->count++;
	switch(before) {
		case 'A':
			n->before[0]++;
			break;
		case 'C':
			n->before[1]++;
			break;
		case 'G':
			n->before[2]++;
			break;
		case 'T':
			n->before[3]++;
			break;
		default:
			n->before[4]++;
	}

	switch(after) {
		case 'A':
			n->after[0]++;
			break;
		case 'C':
			n->after[1]++;
			break;
		case 'G':
			n->after[2]++;
			break;
		case 'T':
			n->after[3]++;
			break;
		default:
			n->after[4]++;
	}


	primer_counts_t *new;
	new = malloc(sizeof(primer_counts_t *) + sizeof(char*) + sizeof(int) + (2 * sizeof(int[5])));
	new->id = NULL;
	new->count = 0;
	for (int i=0; i<5; i++) {
		new->before[i] = 0;
		new->after[i] = 0;
	}
	n->next_primer = new;
	return;
}


float fraction(int n, int tot) {
	return ((float) n / tot) * 100;
}

void print_primers(primer_counts_t *pc, int min_occurrences) {
	primer_counts_t *n  = pc;
	int skippedn = 0;
	int skippedc = 0;
	while (n->id != NULL) {
		if (n->count > min_occurrences) {
			// count the befores/afters
			int tb = 0;
			int ta = 0;
			for (int i = 0; i<5; i++) {
				tb += n->before[i];
				ta += n->after[i];
			}

			printf("Adapter %s. Occurrence: %d\n", n->id, n->count);
			printf("Bases before: A: %d (%.1f%%) C: %d (%.1f%%) G: %d (%.1f%%) T: %d (%.1f%%) Other: (%d) (%.1f%%)\n",
					n->before[0], fraction(n->before[0], tb), n->before[1], fraction(n->before[1], tb),
					n->before[2], fraction(n->before[2], tb), n->before[3], fraction(n->before[3], tb),
					n->before[4], fraction(n->before[4], tb));

			printf("Bases after: A: %d (%.1f%%) C: %d (%.1f%%) G: %d (%.1f%%) T: %d (%.1f%%) Other: (%d) (%.1f%%)\n",
					n->after[0], fraction(n->after[0], ta), n->after[1], fraction(n->after[1], ta),
					n->after[2], fraction(n->after[2], ta), n->after[3], fraction(n->after[3], ta),
					n->after[4], fraction(n->after[4], ta));
		} else {
			skippedc += n->count;
			skippedn++;
		}
		n = n->next_primer;
	}
	printf("There were %d adapters that matched a total of %d times, but we didn't report because they were found fewer than %d times each.\n", skippedn, skippedc, min_occurrences);
}


