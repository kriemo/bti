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
#ifndef BAM_READ_IDX_GET
#define BAM_READ_IDX_GET

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/bgzf.h>

// retrieve pointers to the range of records for tagvalue
// start and end will be NULL if tagvalue is not in the index
// otherwise start will point at the first record with tagvalue 
// and end will point to one-past the last record with tagvalue
void bam_read_idx_get_range(const bam_read_idx* bti, 
                            const char* tagvalue, 
                            bam_read_idx_record** start, 
                            bam_read_idx_record** end);

// fill in the bam record (b) by seeking to the right offset in fp using the information stored in bti_record
void bam_read_idx_get_by_record(htsFile* fp, bam_hdr_t* hdr, bam1_t* b, bam_read_idx_record* bti_record);

// main of the "get" subprogram
int bam_read_idx_get_main(int argc, char** argv);

#endif
