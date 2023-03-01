//---------------------------------------------------------
// Copyright 2019 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// Modified by kriemo 2023 to extract by tag rather than by readname
//
// bti - simple utility to provide random access to
//       bam records by tag value
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include "bti_index.h"

//
// Getopt
//
enum {
    OPT_HELP = 1,
};

static const char* shortopts = ""; // placeholder
static const struct option longopts[] = {
    { "help",                      no_argument,       NULL, OPT_HELP },
    { NULL, 0, NULL, 0 }
};

//
void print_usage_show()
{
    fprintf(stderr, "usage: bti show <index_filename.bti>\n");
}

int bam_read_idx_show_main(int argc, char** argv)
{
    int die = 0;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        switch (c) {
            case OPT_HELP:
                print_usage_show();
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 1) {
        fprintf(stderr, "bti show: not enough arguments\n");
        die = 1;
    }

    if(die) {
        print_usage_show();
        exit(EXIT_FAILURE);
    }

    char* input_bti = argv[optind++];
    bam_read_idx* bti = bam_read_idx_load(NULL, input_bti);

    for(size_t i = 0; i < bti->record_count; ++i) {
        printf("%s %zu\n", bti->records[i].read_name.ptr, bti->records[i].n_aln);
    }

    return 0;
}
