#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "print-sequences.h"
#include "colours.h"

// ALternate DNA encoding lookup table
//                         A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z     
static int dnaEncodeTable [26] = {0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0};
static int dnaDecodeTable [26] = {65, 67, 71, 84};



int encode_base(int base) {
	/*
	 * Convert a base (A, G, C, T) to a number.
	 * Note, you can use seq[i]
	 *
	 * Encoding:
	 * 	A : 0 : 00
	 * 	C : 1 : 01
	 * 	G : 2 : 10 
	 * 	T : 3 : 11
	 *
	 */
	if ((base >= (int)'a') && (base <= (int)'z')) {
		return dnaEncodeTable[base - (int)'a'];
	}
	if ((base >= (int)'A') && (base <= (int)'Z')) {
		return dnaEncodeTable[base - (int)'A'];
	}
	fprintf(stderr, "We can't encode a base that is not [a-z][A-Z]. We have |%c|\n", (char) base);
	return 0;
}

int decode_base(int val) {

	if (val < 0 || val > 3) {
		fprintf(stderr, "We can't decode a value of %d\n", val);
		return 78;
	}
	return dnaDecodeTable[val];
}



uint64_t kmer_encoding(char * seq, int start_position, int k) {
	/*
	 * Given a sequence, seq, start at start_position (0 indexed), and read k characters. 
	 * Convert that to a 64-bit int using 2-bit encoding
	 *
	 * The maximum k-mer length (k-start_position is 32)
	 */

	if ((k - start_position) > 32) {
		fprintf(stderr, "%s We can only encode k<=32 strings in 64 bits. Please reduce k %s\n", PINK, ENDC);
		exit(-1);
	}

	uint64_t enc = 0;


	for (int i=start_position; i < start_position+k; i++)
		enc = (enc << 2) + encode_base(seq[i]);

	return enc;
}

char* kmer_decoding(uint64_t enc, int k) {
	/* 
	 * convert an encoded string back to a base
	 */

	char* seq = malloc(sizeof(char *) * k);
	int posn = k;
	while (posn >= 0) {
		int l = (enc & 1);
		enc >>= 1;
		int h = (enc & 1);
		enc >>= 1;
		h = (h << 1) + l;
		seq[posn--] = (char) decode_base(h);
	}
	return seq;
}



uint64_t next_kmer_encoding(char* seq, int start_position, int k, uint64_t enc) {
	/*
	 * Given a sequence, a start position, k-mer size and a previous encoding
	 * calculate the encoding that starts at start_position but remove the base
	 * at seq[start_position-1]
	 */

	if ((start_position + k - 1) > strlen(seq)) {
		fprintf(stderr, "%s Can't calculate a k-mer beyond the end of the sequence. Start: %d, k-mer %d, sequence length %ld %s\n", RED, start_position, k, strlen(seq), ENDC);
		exit(2);
	}

	enc = enc - ((uint64_t) encode_base(seq[start_position - 1]) << (2 * (k-1)));
	enc = (enc << 2) + encode_base(seq[start_position + k - 1]);
	return enc;
}





