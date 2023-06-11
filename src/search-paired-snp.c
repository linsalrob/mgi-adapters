/*
 * Search for adapters using paired end files
 *
 * (c) Rob. 2023
 */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <getopt.h>
#include <string.h>
#include "match-paired-files.h"
#include "colours.h"
#include "version.h"

void help() {
	printf("USAGE: search-paired-snp -l -r -s -f -p -q\n");
	printf("\nThis version attempts to look for SNPs at each position of your adapters\n\n");
	printf("-1 --R1 R%s1%s file\n", GREEN, ENDC);
	printf("-2 --R2 R%s2%s file\n", GREEN, ENDC);
	printf("-s --I7left I7left primer sequence (%ss%seven). Default 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'\n", GREEN, ENDC);
	printf("-f --I5right I5right primer sequence (%sf%sive). Default 'ACACTCTTTCCCTACACGACGCTCTTCCGATC'\n", GREEN, ENDC);
	printf("-j --matchesR1 Write the R1 matches to this file. Default: stdout\n");
	printf("-k --matchesR2 Write the R2 matches to this file. Default: stdout\n");
	printf("-a --adjustments Write the trimming adjustments here\n");
	printf("--verbose more output\n");
	printf("--debug lots of output\n");
}


int main(int argc, char* argv[]) {
	if (argc == 2 && (strcmp(argv[1], "-V") == 0)) {
                printf("Version: %f\n", __version__);
		exit(0);
	}

	if (argc < 3) {
		help(argv[0]);
		exit(0);
	}


	struct options *opt;
	opt = malloc(sizeof(struct options));

	opt->tablesize = 100003;
	opt->R1_file = NULL;
	opt->R2_file = NULL;
	opt->R1_output = NULL;
	opt->R2_output = NULL;
	opt->R1_matches = NULL;
	opt->R2_matches = NULL;
	opt->I7left = NULL;
	opt->I5right = NULL;
	opt->debug = false;
	opt->verbose = false;
	opt->adjustments = NULL;

    int gopt = 0;
    static struct option long_options[] = {
            {"R1",  required_argument, 0, 'l'},
            {"R2",  required_argument, 0, 'r'},
            {"outputR1",  required_argument, 0, 'p'},
            {"outputR2",  required_argument, 0, 'q'},
            {"matchesR1",  required_argument, 0, 'j'},
            {"matchesR2",  required_argument, 0, 'k'},
            {"I7left",  required_argument, 0, 's'},
            {"I5right", required_argument, 0, 'f'},
	    {"adjustments", required_argument, 0, 'a'},
            {"debug", no_argument, 0, 'd'},
	    {"verbose", no_argument, 0, 'b'},
            {"version", no_argument, 0, 'v'},
            {0, 0, 0, 0}
    };
    int option_index = 0;
    while ((gopt = getopt_long(argc, argv, "1:2:p:q:s:f:j:k:a:bdv", long_options, &option_index )) != -1) {
        switch (gopt) {
            case '1' :
                opt->R1_file = strdup(optarg);
                break;
            case '2':
                opt->R2_file = strdup(optarg);
                break;
            case 'p' :
                opt->R1_output = strdup(optarg);
                break;
            case 'q':
                opt->R2_output = strdup(optarg);
                break;
            case 'j' :
                opt->R1_matches = strdup(optarg);
                break;
            case 'k':
                opt->R2_matches = strdup(optarg);
                break;
            case 's' :
                opt->I7left = strdup(optarg);
                break;
            case 'f':
                opt->I5right = strdup(optarg);
                break;
            case 'a':
                opt->adjustments = strdup(optarg);
                break;
            case 'b': opt->verbose = true;
                break;
            case 'd': opt->debug = true;
                break;
            case 'v':
                printf("Version: %f\n", __version__);
                return 0;
	    default: help();
	        exit(EXIT_FAILURE);
	}
    }

    if (opt->I7left == NULL)
	    opt->I7left = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
    if (opt->I5right == NULL)
	    opt->I5right = "ACACTCTTTCCCTACACGACGCTCTTCCGATC";
    if (opt->R1_output == NULL)
	    opt->R1_output = "R1.trimmed.fastq.gz";
    if (opt->R2_output == NULL)
	    opt->R2_output = "R2.trimmed.fastq.gz";

    search_pairwise_snps(opt);
}

