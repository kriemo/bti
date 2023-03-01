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
#include "bti_index.h"
#include "bti_get.h"
#include "bti_show.h"
#include "bti_test.h"

#define bti_VERSION "0.3"

void print_version()
{
    printf("%s\n", bti_VERSION);
}

int main(int argc, char** argv)
{
    if(argc < 2) {
        fprintf(stderr, "[bti] usage bti <subprogram> [...]\n");
        exit(EXIT_FAILURE);
    }

    if(strcmp(argv[1], "index") == 0) {
        bam_read_idx_index_main(argc - 1, argv + 1);
    } else if(strcmp(argv[1], "get") == 0) {
       bam_read_idx_get_main(argc - 1, argv + 1);
    } else if(strcmp(argv[1], "show") == 0) {
       bam_read_idx_show_main(argc - 1, argv + 1);
    } else if(strcmp(argv[1], "test") == 0) {
        bam_read_idx_test_main(argc - 1, argv + 1);
    } else if(strcmp(argv[1], "version") == 0) {
        print_version();
    } else {
        fprintf(stderr, "[bti] unrecognized subprogram: %s\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    return 0;
}

