/*
 * Search for all adapters using a fasta file of adapters and paired end files
 *
 * (c) Rob. 2023
 */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <getopt.h>
#include <string.h>
#include "match-all-snps.h"
#include "colours.h"
#include "version.h"

void help() {
	printf("USAGE: search-paired-snp -l -r -s -f -p -q\n");
	printf("\nSearch for primers listed in %s--primers%s, allowing for 1-bp mismatches, against all the reads in %s--R1%s and %s--R2%s\n", 
			GREEN, ENDC, GREEN, ENDC, GREEN, ENDC);
	printf("-1 --R1 R%s1%s file (%srequired%s)\n", GREEN, ENDC, RED, ENDC);
	printf("-2 --R2 R%s2%s file (%srequired%s)\n", GREEN, ENDC, RED, ENDC);
	printf("-f --primers fasta file of primers (%srequired%s)\n", RED, ENDC);
	printf("-p --outputR1 R1 output fastq file (will be gzip compressed)\n");
	printf("-q --outputR2 R2 output fastq file (will be gzip compressed)\n");
	printf("-j --matchesR1 Write the R1 matches to this file. Default: stdout\n");
	printf("-k --matchesR2 Write the R2 matches to this file. Default: stdout\n");
	printf("-n --noreverse Do not reverse the sequences\n");
	printf("-a --adjustments Write the trimming adjustments here\n");
	printf("--verbose more output (but less than --debug)\n");
	printf("--debug more more output\n");
	printf("-v --version print the version and exit\n");
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
	opt->reverse = true;
	opt->primers = NULL;
	opt->debug = false;
	opt->verbose = false;
	opt->adjustments = NULL;

    int gopt = 0;
    static struct option long_options[] = {
            {"R1",  required_argument, 0, '1'},
            {"R2",  required_argument, 0, '2'},
            {"primers",  required_argument, 0, 'f'},
            {"outputR1",  required_argument, 0, 'p'},
            {"outputR2",  required_argument, 0, 'q'},
            {"matchesR1",  required_argument, 0, 'j'},
            {"matchesR2",  required_argument, 0, 'k'},
	    {"noreverse", required_argument, 0, 'n'},
	    {"adjustments", required_argument, 0, 'a'},
            {"debug", no_argument, 0, 'd'},
            {"version", no_argument, 0, 'v'},
	    {"verbose", no_argument, 0, 'b'},
            {0, 0, 0, 0}
    };
    int option_index = 0;
    while ((gopt = getopt_long(argc, argv, "1:2:p:q:f:j:k:a:bndv", long_options, &option_index )) != -1) {
        switch (gopt) {
            case '1' :
                opt->R1_file = strdup(optarg);
                break;
            case '2':
                opt->R2_file = strdup(optarg);
                break;
            case 'p':
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
            case 'n' :
                opt->reverse = false;
                break;
            case 'f':
                opt->primers = strdup(optarg);
                break;
            case 'a':
                opt->adjustments = strdup(optarg);
                break;
            case 'd': 
		opt->debug = true;
                break;
	    case 'b':
		opt->verbose = true;
		break;
            case 'v':
                printf("Version: %f\n", __version__);
                return 0;
	    default: help();
	        exit(EXIT_FAILURE);
	}
    }

    if (opt->R1_file == NULL) {
	    fprintf(stderr, "Please provide an R1 read file\n");
	    help();
	    exit(EXIT_FAILURE);
    }
    if (opt->R2_file == NULL) {
	    fprintf(stderr, "Please provide an R22 read file\n");
	    help();
	    exit(EXIT_FAILURE);
    }
    if (opt->primers == NULL) {
	    fprintf(stderr, "Please provide a primer file\n");
	    help();
	    exit(EXIT_FAILURE);
    }


    if (opt->R1_output == NULL)
	    opt->R1_output = "R1.trimmed.fastq.gz";
    if (opt->R2_output == NULL)
	    opt->R2_output = "R2.trimmed.fastq.gz";

    search_all_pairwise_snps(opt);
    free(opt);
}

