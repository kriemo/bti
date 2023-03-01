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
#ifndef BAM_READ_IDX_INDEX
#define BAM_READ_IDX_INDEX

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/bgzf.h>

// An entry record in the index, storing
// either an offset or pointer to the tagvalue
// and a position in the bgzf-compressed bam file.
// Also stores the # of alignments with a given tag/tagvalue
typedef struct bam_read_idx_record
{
    // When building the index and storing it on disk
    // the read name for this record is stored as an
    // offset into tagvalues. When it is loaded from
    // disk the size of the index is fixed and we
    // convert the offset into a direct pointer.
    union read_name {
        size_t offset;
        const char* ptr;
    } read_name;

    size_t file_offset;
    size_t n_aln;
} bam_read_idx_record;

//
// The index itself consists of two parts,
//  1) a memory block containing the names of every unique indexed read/tag
//  2) records (see above) describing the position in the file for each alignment
//
typedef struct bam_read_idx
{
    // read names are stored contiguously in one large block
    // of memory as null terminated strings
    size_t name_capacity_bytes;
    size_t name_count_bytes;
    char* tagvalues;

    // records giving the offset into the bam file
    size_t record_capacity;
    size_t record_count;
    bam_read_idx_record* records;
} bam_read_idx;

// load the index for input_bam file
// returns a pointer to the index, which must be deallocated by
// the caller using bam_read_idx_destroy.
bam_read_idx* bam_read_idx_load(const char* input_bam, const char* input_bti);

// construct the index for input_bam and save it to disk
// to use the created index bam_read_idx_load should be called
void bam_read_idx_build(const char* input_bam, const char* output_bti, const char* tag);

// cleanup the index by deallocating everything
void bam_read_idx_destroy(bam_read_idx* bti);

// main of subprogram - build an index for a file
// given in the command line arguments
int bam_read_idx_index_main(int argc, char** argv);

#endif
