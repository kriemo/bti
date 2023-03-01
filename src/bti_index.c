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

// avoid warnings in qsort_r
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include "bti_index.h"
#include "sort_r.h"

//#define bti_INDEX_DEBUG 1
char verbose = 0;

// make the index filename based on the name of input_bam
// caller must free the returned pointer
char* generate_index_filename(const char* input_bam, const char* input_bti) 
{
    char* out_fn;

    if(input_bti != NULL) {
        out_fn = malloc(strlen(input_bti));
        if(out_fn == NULL) {
            exit(EXIT_FAILURE);
        }
        strcpy(out_fn, input_bti);
    } else {
        out_fn = malloc(strlen(input_bam) + 5);
        if(out_fn == NULL) {
            exit(EXIT_FAILURE);
        }
        strcpy(out_fn, input_bam);
        strcat(out_fn, ".bti");
    }

    return out_fn;
}

//
bam_read_idx* bam_read_idx_init()
{
    bam_read_idx* bti = (bam_read_idx*)malloc(sizeof(bam_read_idx));

    bti->name_capacity_bytes = 0;
    bti->name_count_bytes = 0;
    bti->readnames = NULL;

    bti->record_capacity = 0;
    bti->record_count = 0;
    bti->records = NULL;

    return bti;
}

//
void bam_read_idx_destroy(bam_read_idx* bti)
{
#ifdef bti_INDEX_DEBUG
    fprintf(stderr, "[bti-destroy] %zu name bytes %zu records\n", bti->name_count_bytes, bti->record_count);
#endif
    free(bti->readnames);
    bti->readnames = NULL;

    free(bti->records);
    bti->records = NULL;

    free(bti);
}

// comparison function used by quicksort, names pointers to a block of C strings
// and allows us to indirectly sorted records by their offset.
int compare_records_by_readname_offset(const void* r1, const void* r2, void* names)
{
    const char* cnames = (const char*)names;
    const char* n1 = cnames + ((bam_read_idx_record*)r1)->read_name.offset;
    const char* n2 = cnames + ((bam_read_idx_record*)r2)->read_name.offset;
    return strcmp(n1, n2);
}

//
void bam_read_idx_save(bam_read_idx* bti, const char* filename)
{
    FILE* fp = fopen(filename, "wb");

    // Sort records by readname
    sort_r(bti->records, bti->record_count, sizeof(bam_read_idx_record), compare_records_by_readname_offset, bti->readnames);
    
    // write header, containing file version, the size (in bytes) of the read names
    // and the number of records. The readnames size is a placeholder and will be
    // corrected later.
    
    // version
    size_t FILE_VERSION = 1;
    fwrite(&FILE_VERSION, sizeof(FILE_VERSION), 1, fp);
    
    // readname length
    size_t readname_bytes = 0;
    fwrite(&readname_bytes, sizeof(readname_bytes), 1, fp);

    // num records
    fwrite(&bti->record_count, sizeof(bti->record_count), 1, fp);

    // Pass 1: count up the number of non-redundant read names, write them to disk
    // Also store the position in the file where the read name for each record was written
    size_t* disk_offsets_by_record = malloc(bti->record_count * sizeof(size_t));
    const char* rn = bti->readnames; // for convenience 

    for(size_t i = 0; i < bti->record_count; ++i) {
        
        int redundant = i > 0 && 
            strcmp(rn + bti->records[i].read_name.offset, rn + bti->records[i - 1].read_name.offset) == 0;
        
        if(!redundant) {
            disk_offsets_by_record[i] = readname_bytes; // current position in file
            size_t len = strlen(rn + bti->records[i].read_name.offset) + 1;
            fwrite(rn + bti->records[i].read_name.offset, len, 1, fp);
            readname_bytes += len;
        } else {
            disk_offsets_by_record[i] = disk_offsets_by_record[i - 1];
        }

#ifdef bti_INDEX_DEBUG
        fprintf(stderr, "record %zu name: %s redundant: %d do: %zu offset: %zu\n", 
            i, bti->readnames + bti->records[i].read_name.offset, redundant, disk_offsets_by_record[i], bti->records[i].file_offset);
#endif
    }

    // Pass 2: write the records, getting the read name offset from the disk offset (rather than
    // the memory offset stored)
    for(size_t i = 0; i < bti->record_count; ++i) {
        bam_read_idx_record btir = bti->records[i];
        btir.read_name.offset = disk_offsets_by_record[i];
        fwrite(&btir, sizeof(btir), 1, fp);
#ifdef bti_INDEX_DEBUG
        fprintf(stderr, "[bti-save] record %zu %s name offset: %zu file offset: %zu\n", 
            i, bti->readnames + bti->records[i].read_name.offset, disk_offsets_by_record[i], bti->records[i].file_offset);
#endif
    }
    
    // finish by writing the actual size of the read name segment
    fseek(fp, sizeof(FILE_VERSION), SEEK_SET);
    fwrite(&readname_bytes, sizeof(readname_bytes), 1, fp);

    free(disk_offsets_by_record);
    fclose(fp);
}

// add a record to the index, growing the dynamic arrays as necessary
void bam_read_idx_add(bam_read_idx* bti, const char* readname, size_t offset)
{
    // 
    // add readname to collection if not seen, otherwise keep original offset
    //

    if(bti->record_count > 0){
        const char* orn = bti->readnames + bti->records[bti->record_count - 1].read_name.offset ;
        if(strcmp(orn, readname) == 0) {
            bti->records[bti->record_count - 1].n_aln += 1; 
          //  bti->records[bti->record_count - 1].file_offset = offset;
           // printf("%s %s\n", bti->readnames + bti->records[bti->record_count - 1].read_name.offset, readname);
            return;
        }
    }
    
    size_t len = strlen(readname) + 1;
    if(bti->name_capacity_bytes <= bti->name_count_bytes + len) {

        // if already allocated double size, if initialization start with 1Mb
        bti->name_capacity_bytes = bti->name_capacity_bytes > 0 ? 2 * bti->name_capacity_bytes : 1024*1024;
#ifdef bti_INDEX_DEBUG
        fprintf(stderr, "[bti] allocating %zu bytes for names\n", bti->name_capacity_bytes);
#endif
    
        bti->readnames = realloc(bti->readnames, bti->name_capacity_bytes);
        if(bti->readnames == NULL) {
            fprintf(stderr, "[bti] malloc failed\n");
            exit(EXIT_FAILURE);
        }
    }
    
    // in principle this incoming name can be so larger than doubling the size
    // doesn't allow it to fit. This really shouldn't happen so we'll just exit here
    if(bti->name_capacity_bytes <= bti->name_count_bytes + len) {
        fprintf(stderr, "[bti] incoming name with length %zu is too large (%zu %zu)\n", len, bti->name_count_bytes + len, bti->name_capacity_bytes);
        exit(EXIT_FAILURE);
    }
    
    // copy name
    size_t name_offset = bti->name_count_bytes;
    strncpy(bti->readnames + bti->name_count_bytes, readname, len);
    bti->name_count_bytes += len;
    assert(bti->readnames[bti->name_count_bytes - 1] == '\0');

    //
    // add record
    //
    if(bti->record_count == bti->record_capacity) {
        bti->record_capacity = bti->record_capacity > 0 ? 2 * bti->record_capacity : 1024;
        bti->records = realloc(bti->records, sizeof(bam_read_idx_record) * bti->record_capacity);
        if(bti->records == NULL) {
            fprintf(stderr, "[bti] malloc failed\n");
            exit(EXIT_FAILURE);
        }
    }
    bti->records[bti->record_count].n_aln = 1;
    bti->records[bti->record_count].read_name.offset = name_offset;
    bti->records[bti->record_count].file_offset = offset;
    bti->record_count += 1;
    
}

//
void bam_read_idx_build(const char* filename, const char* output_bti, const char* tag)
{
    htsFile *fp = hts_open(filename, "r");
    if(fp == NULL) {
        fprintf(stderr, "[bti] could not open %s\n", filename);
        exit(EXIT_FAILURE);
    }

    bam_read_idx* bti = bam_read_idx_init();

    bam1_t* b = bam_init1();
    bam_hdr_t *h = sam_hdr_read(fp);
    int ret = 0;
    int niter = 0;
    size_t file_offset = bgzf_tell(fp->fp.bgzf);
    while ((ret = sam_read1(fp, h, b)) >= 0) {
        char* readname ;
        uint8_t* aux_info = bam_aux_get(b, tag) ;
        if (aux_info) {
            readname = bam_aux2Z(aux_info) ;
        } else {
            file_offset = bgzf_tell(fp->fp.bgzf);
            continue;
        }
        
        bam_read_idx_add(bti, readname, file_offset);

        bam_read_idx_record btir = bti->records[bti->record_count - 1];
        if(verbose && (niter == 1 || niter % 100000 == 0)) {
            fprintf(stderr, "[bti-build] record %zu [%zu %zu] %s\n",
                bti->record_count,
                btir.read_name.offset,
                btir.file_offset,
                bti->readnames + btir.read_name.offset
            );
        }
        niter += 1;
        // update offset for next record
        file_offset = bgzf_tell(fp->fp.bgzf);
    }

    bam_hdr_destroy(h);
    bam_destroy1(b);
    hts_close(fp);

    // save to disk and cleanup
    if(verbose) {
        fprintf(stderr, "[bti-build] writing to disk...\n");
    }

    char* out_fn = generate_index_filename(filename, output_bti);
    bam_read_idx_save(bti, out_fn);

    if(verbose) {
        fprintf(stderr, "[bti-build] wrote index for %zu records.\n", bti->record_count);
    }

    free(out_fn);
    bam_read_idx_destroy(bti);
}

void print_error_and_exit(const char* msg)
{
    fprintf(stderr, "[bti] %s\n", msg);
    exit(EXIT_FAILURE);
}

//
bam_read_idx* bam_read_idx_load(const char* input_bam, const char* input_bti)
{
    char* index_fn = generate_index_filename(input_bam, input_bti);
    FILE* fp = fopen(index_fn, "rb");
    if(fp == NULL) {
        fprintf(stderr, "[bti] index file not found for %s\n", input_bam);
        exit(EXIT_FAILURE);
    }

    bam_read_idx* bti = bam_read_idx_init();
    size_t file_version;
    // currently ignored
    int bytes_read = fread(&file_version, sizeof(file_version), 1, fp);
    if(bytes_read <= 0) {
        print_error_and_exit("read error");
    }

    // read size of readames segment
    bytes_read = fread(&bti->name_count_bytes, sizeof(bti->name_count_bytes), 1, fp);
    if(bytes_read <= 0) {
        print_error_and_exit("read error");
    }
    bti->name_capacity_bytes = bti->name_count_bytes;

    // read number of records on disk
    bytes_read = fread(&bti->record_count, sizeof(bti->record_count), 1, fp);
    if(bytes_read <= 0) {
        print_error_and_exit("read error");
    }

    bti->record_capacity = bti->record_count;

    // allocate filenames
    bti->readnames = malloc(bti->name_capacity_bytes);
    if(bti->readnames == NULL) {
        fprintf(stderr, "[bti] failed to allocate %zu bytes for read names\n", bti->name_capacity_bytes);
        exit(EXIT_FAILURE);
    }

    // allocate records
    bti->records = malloc(bti->record_capacity * sizeof(bam_read_idx_record));

    // read the names
    bytes_read = fread(bti->readnames, bti->name_count_bytes, 1, fp);
    if(bytes_read <= 0) {
        print_error_and_exit("read error");
    }

    // read the records
    bytes_read = fread(bti->records, bti->record_count, sizeof(bam_read_idx_record), fp);
    if(bytes_read <= 0) {
        print_error_and_exit("read error");
    }

    // convert read name offsets to direct pointers
    for(size_t i = 0; i < bti->record_count; ++i) {
        bti->records[i].read_name.ptr = bti->readnames + bti->records[i].read_name.offset;
#ifdef bti_INDEX_DEBUG
        fprintf(stderr, "[bti-load] record %zu %s %zu\n", i, bti->records[i].read_name.ptr, bti->records[i].file_offset);
#endif
    }

    fclose(fp);
    free(index_fn);
    return bti;
}

//
// Getopt
//
enum {
    OPT_HELP = 1,
};

static const char* shortopts = ":i:v"; // placeholder
static const struct option longopts[] = {
    { "help",                      no_argument,       NULL, OPT_HELP },
    { "index",               required_argument,       NULL,      'i' },
    { "verbose",                   no_argument,       NULL,      'v' },
    { NULL, 0, NULL, 0 }
};

//
void print_usage_index()
{
    fprintf(stderr, "usage: bti index [-v] [-i <index_filename.bti>] <tag> <input.bam>\n");
}

//
int bam_read_idx_index_main(int argc, char** argv)
{
    char* output_bti = NULL;
    
    int die = 0;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        switch (c) {
            case OPT_HELP:
                print_usage_index();
                exit(EXIT_SUCCESS);
            case 'i':
                output_bti = optarg;
                break;
            case 'v':
                verbose = 1;
                break;
        }
    }
    
    if (argc - optind < 2) {
        fprintf(stderr, "bti index: not enough arguments\n");
        die = 1;
    }

    if(die) {
        print_usage_index();
        exit(EXIT_FAILURE);
    }
    
    char* tag = argv[optind++]; 
    char* input_bam = argv[optind++];
    bam_read_idx_build(input_bam, output_bti, tag);

    return 0;
}
