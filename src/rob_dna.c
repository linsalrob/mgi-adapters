#include "rob_dna.h"
#include "colours.h"
#include <stdint.h>
#include <string.h>
#include <stdio.h>

uint64_t reverse_complement(uint64_t enc, int len) {
	/*
	 * Given a primer encoding and the length of the primer, 
	 * calculate its reverse complement. Note that we need
	 * the length otherwise all the leading 00's become 11's
	 */

        uint64_t r = 0;
        for (int i=0; i<len; i++) {
                //left shift, first time, its all 0's
		r<<=2;
                // get the value of the two lsbs
                int l = (enc & 1);
                enc >>= 1;
                int h = (enc & 1);
                enc >>= 1;
                // are low (l) and high (h) set, 
		// if l and h: (11) T -> need A (0)
                // if l not h: (01) C -> need G (2)
                // if h not l: (10) G -> need C (1)
                // otherwise need (00) T (3)
                int t = l & h ? 0 : l ? 2 : h ? 1 : 3;
                r += t;
        }
        return r;
}

char* rc(char* to, char* from) {
	/*
	 * Reverse complement the sequence in from and put it in to
	 */

	//            A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
	int rc[26] = {3, 3, 2, 3, 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 3};
	char*bases = "ACGT";

	int i;
	for (i = 0; i<strlen(from); i++)
		to[strlen(from)-i-1] = bases[rc[(int) from[i] - (int)'A']];
	to[i]='\0';
	return to;
}



int has_n(char * seq) {
	/*
	 * Does a sequence contain an N?
	 */

	char * i;

	for (i=seq; *i; i++) {
		// i points successively to a[0], a[1], ... until *i is '\0' 
		if (*i == 78 || *i == 110)
			return 1;
	}
	return 0;
}
