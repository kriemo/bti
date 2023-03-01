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
#include "bti_get.h"

enum {
    OPT_HELP = 1,
};

static const char* shortopts = ":i:"; // placeholder
static const struct option longopts[] = {
    { "help",                      no_argument,       NULL, OPT_HELP },
    { NULL, 0, NULL, 0 }
};

void print_usage_test()
{
    fprintf(stderr, "usage: bti test [-i <index_filename.bti>] <tag> <input.bam>\n");
}

//
int bam_read_idx_test_main(int argc, char** argv)
{
    char* input_bti = NULL;

    int die = 0;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        switch (c) {
            case OPT_HELP:
                print_usage_test();
                exit(EXIT_SUCCESS);
            case 'i':
                input_bti = optarg;
        }
    }
    
    if (argc - optind < 2) {
        fprintf(stderr, "bti test: not enough arguments\n");
        die = 1;
    }

    if(die) {
        print_usage_test();
        exit(EXIT_FAILURE);
    }
    char* tag = argv[optind++];
    char* input_bam = argv[optind++];

    // open files
    bam_read_idx* bti = bam_read_idx_load(input_bam, input_bti);
    htsFile* bam_fp = hts_open(input_bam, "r");
    bam_hdr_t* h = sam_hdr_read(bam_fp);
    bam1_t* b = bam_init1();

    // iterate over each record and run get on each readname
    const char* prev_readname = NULL;
    for(size_t ri = 0; ri < bti->record_count; ++ri) {
        const char* readname = bti->records[ri].read_name.ptr;

        // skip if same as previous readname
        if(readname == prev_readname) {
            continue;
        }

        bam_read_idx_record* start;
        bam_read_idx_record* end;
        bam_read_idx_get_range(bti, readname, &start, &end);
        while(start != end) {
        
            bam_read_idx_get_by_record(bam_fp, h, b, start);
            char* cur_tag ;
            uint8_t* aux_info = bam_aux_get(b, tag) ;
            if (aux_info) {
              cur_tag = bam_aux2Z(aux_info) ;
            } else {
              continue;
            }

            fprintf(stderr, "[bti-test] %s %s\n", readname, cur_tag);
            assert(strcmp(readname, cur_tag) == 0);

            // mark this record as used so we can make sure every record is present in the bam
            // this destroys the index
            start->file_offset = 0;
            start++;
        }

        prev_readname = readname;
    }

    // check that all records were accessed
    for(size_t ri = 0; ri < bti->record_count; ++ri) {
        assert(bti->records[ri].file_offset == 0);
    }
    
    bam_destroy1(b);
    bam_hdr_destroy(h);
    hts_close(bam_fp);
    bam_read_idx_destroy(bti);
    bti = NULL;

    return 0;
}
