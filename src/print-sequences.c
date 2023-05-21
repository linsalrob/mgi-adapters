#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

char* int_to_binary(uint64_t enc) {
	/*
	 * print a binary representation of the ulli
	 */

	char* bt = (char*) malloc(sizeof(char) * 64);
	int k = 0;
	uint64_t i = (1ULL << 62);
	for (uint64_t i = (1ULL << 62); i > 0; i = i / 2) {
		bt[k++] = (enc & i) ? '1' : '0';
	}
	bt[k] = '\0';
	return bt;
}

char* substr(char* seq, int start, int k) {
	if ((start + k - 1) > strlen(seq)) {
		fprintf(stderr, "Can't get a substring longer than the string\n");
		k = strlen(seq);
	}
	char* ss = malloc(sizeof(char) * strlen(seq));
	int j = 0;
	for (int i = start; i <= (start + k - 1); i++)
		ss[j++] = seq[i];
	ss[j] = '\0';
	return ss;
}
