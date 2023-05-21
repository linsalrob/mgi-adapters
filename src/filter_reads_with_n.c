/* 
 * Filter fastq sequences that don't have an N
 *
 * I am also trying to implementa faster gzip writer
 */


#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <zlib.h>
#include <string.h>
#include "colours.h"
#include "rob_dna.h"
#include "kseq.h"
#include "version.h"

KSEQ_INIT(gzFile, gzread);

void print_usage() {
	printf("Usage: filter_reads_with_n -f sequences.fastq.gz -o output.fastq.gz\n");
}


int main(int argc, char *argv[]) {
	char* fqfile = NULL;
	char* outputfile = NULL;
	static struct option long_options[] = {};
	int opt = 0;
        int option_index = 0;
	bool errors = false;
        while ((opt = getopt_long(argc, argv, "f:o:ev", long_options, &option_index)) != -1) {
                switch (opt) {
                        case 'f' :
				fqfile = strdup(optarg);
                                break;
                        case 'o' :
				outputfile = strdup(optarg);
                                break;
			case 'e' :
				errors = true;
				break;
                        case 'v':
                                printf("Version: %f\n", __version__);
                                return 0;
                        default:
                                print_usage();
                                exit(EXIT_FAILURE);
                }
        }

	if (outputfile == NULL && fqfile == NULL) {
		print_usage();
		exit(EXIT_FAILURE);
	}


	if( access( fqfile, R_OK ) == -1 ) {
		// file doesn't exist
		fprintf(stderr, "%sERROR: The file %s can not be found. Please check the file path%s\n", RED, fqfile, ENDC);
		return;
	}

	gzFile fp;
	kseq_t *seq;

	fp = gzopen(fqfile, "r");
	seq = kseq_init(fp);
	int l;

	// We are using a pipe to gzip rather than libz. The advantage is that we gain 
	// parallel processing because the gzip happens on a separate thread.
	// create a file name that includes the pipe
	char* pipe_file = malloc(sizeof(char) * (strlen(outputfile) + 10));
	strcpy(pipe_file, "gzip - > ");
	strcat(pipe_file, outputfile);

	FILE *pipe = popen(pipe_file, "w");
	while ((l = kseq_read(seq)) >= 0) {
		if (has_n(seq->seq.s)) {
			if (errors)
				printf("%s has an N\n", seq->name.s);
			continue;
		}
		fprintf(pipe, "@%s %s\n%s\n+\n%s\n", seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);
	}
	pclose(pipe);


	kseq_destroy(seq);
	gzclose(fp);
}

