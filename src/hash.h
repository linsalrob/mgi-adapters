#ifndef HASH_H
#define HASH_H


/*
 * calculate the hash for a fastq sequence
 *
 * This is a simple hash but widely used!
 *
 * we use an unsigned here so that the answer is > 0
 *
 * You still need to mod this on the table size
 */

unsigned hash (char *s);

#endif
