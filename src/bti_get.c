//---------------------------------------------------------
// Copyright 2019 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
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

static const char* shortopts = ":i::"; // placeholder
static const struct option longopts[] = {
    { "help",                      no_argument,       NULL, OPT_HELP },
    { "index",               required_argument,       NULL,      'i' },
    { NULL, 0, NULL, 0 }
};

void print_usage_get()
{
    fprintf(stderr, "usage: bti get [-i <index_filename.bti>] <input.bam> <tag values>  \n");
}

// comparator used by bsearch, direct strcmp through the name pointer
int compare_records_by_tagvalue_ptr(const void* r1, const void* r2)
{
    const char* n1 = ((bam_read_idx_record*)r1)->read_name.ptr;
    const char* n2 = ((bam_read_idx_record*)r2)->read_name.ptr;
    return strcmp(n1, n2);
}

//
void bam_read_idx_get_range(const bam_read_idx* bti, const char* tagvalue, bam_read_idx_record** start, bam_read_idx_record** end)
{
    // construct a query record to pass to bsearch
    bam_read_idx_record query;
    query.read_name.ptr = tagvalue;
    query.file_offset = 0;

    // if rec is NULL then tagvalue does not appear in index
    bam_read_idx_record* rec = 
        bsearch(&query, bti->records, bti->record_count, sizeof(bam_read_idx_record), compare_records_by_tagvalue_ptr);
    if(rec == NULL) {
        *start = NULL;
        *end = NULL;
        return;
    }

    // rec points to a valid record, but it can be an arbitrary record in the range
    // move start to the first record in the range, and end to be one past the end
    size_t sri, eri;
    sri = eri = rec - bti->records;
    assert(bti->records[sri].file_offset == rec->file_offset);

    while(sri > 0 && bti->records[sri].read_name.ptr == bti->records[sri - 1].read_name.ptr) {
        sri -= 1;
    }
    assert(strcmp(bti->records[sri].read_name.ptr, tagvalue) == 0);

    do {
        eri += 1;
    } while(eri < bti->record_count && bti->records[eri].read_name.ptr == bti->records[sri].read_name.ptr);
    assert(eri == bti->record_count || strcmp(bti->records[eri].read_name.ptr, tagvalue) != 0);
    //fprintf(stderr, "r: %zu sri: %zu eri: %zu\n", rec - bti->records, sri, eri);
    *start = &bti->records[sri];
    *end = &bti->records[eri];
}

//
void bam_read_idx_get_by_record(htsFile* fp, bam_hdr_t* hdr, bam1_t* b, bam_read_idx_record* bti_record)
{
    int ret = bgzf_seek(fp->fp.bgzf, bti_record->file_offset, SEEK_SET);
    if(ret != 0) {
        fprintf(stderr, "[bti] bgzf_seek failed\n");
        exit(EXIT_FAILURE);
    }

    ret = sam_read1(fp, hdr, b);
    if(ret < 0) {
        fprintf(stderr, "[bti] sam_read1 failed\n");
        exit(EXIT_FAILURE);
    }
}

//
int bam_read_idx_get_main(int argc, char** argv)
{
    char* input_bti = NULL;

    int die = 0;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        switch (c) {
            case OPT_HELP:
                print_usage_get();
                exit(EXIT_SUCCESS);
            case 'i':
                input_bti = optarg;
        }
    }
    
    if (argc - optind < 2) {
        fprintf(stderr, "bti get: not enough arguments\n");
        die = 1;
    }

    if(die) {
        print_usage_get();
        exit(EXIT_FAILURE);
    }

    char* input_bam = argv[optind++];
    
    bam_read_idx* bti = bam_read_idx_load(input_bam, input_bti);
    
    htsFile *bam_fp = sam_open(input_bam, "r");
    bam_hdr_t *h = sam_hdr_read(bam_fp);
    htsFile *out_fp ;
    out_fp = hts_open("-", "wb");
    
    sam_hdr_write(out_fp, h);
    
    bam_read_idx_record* start;
    bam_read_idx_record* end;

    for(int i = optind; i < argc; i++) {
        char* tagvalue = argv[i];
        bam_read_idx_get_range(bti, tagvalue, &start, &end);
        
        bam1_t *b = bam_init1();
        int n_rec = 0;
        while(start != end) {
            int ret = bgzf_seek(bam_fp->fp.bgzf , start->file_offset, SEEK_SET);
            if(ret != 0) {
                fprintf(stderr, "[bti] bgzf_seek failed\n");
                exit(EXIT_FAILURE);
            }
            
            while(n_rec < start->n_aln){
                ret = sam_read1(bam_fp, h, b);
                if(ret < 0) {
                    fprintf(stderr, "[bti] sam_read1 failed\n");
                    exit(EXIT_FAILURE);
                }
                int ret = sam_write1(out_fp, h, b);
                if(ret < 0) {
                    fprintf(stderr, "[bti] sam_write1 failed\n");
                    exit(EXIT_FAILURE);
                }
                n_rec += 1;

            }
            start++;
        }
        bam_destroy1(b);
    }

    hts_close(out_fp);
    bam_hdr_destroy(h);
    hts_close(bam_fp);
    bam_read_idx_destroy(bti);
    bti = NULL;

    return 0;
}
