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
#ifndef BAM_READ_IDX_TEST
#define BAM_READ_IDX_TEST

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/bgzf.h>

// main of the "test" subprogram
int bam_read_idx_test_main(int argc, char** argv);

#endif
