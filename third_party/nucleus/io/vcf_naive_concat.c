/*  vcfconcat.c -- Concatenate or combine VCF/BCF files.

    Copyright (C) 2013-2023 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


    Copyright 2024 Google LLC.

Note(lucasbrambrink): This is a trimmed down and edited version of
bcftools/vcfconcat.c to only implement --naive concatenation of VCF files.
*/

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <sys/time.h>
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/tbx.h" // for hts_get_bgzfp()
#include "htslib/hts_endian.h"
#include "htslib/hts_log.h"

#define MAGIC_LEN 9  // 5 + 4 = "Magic" string + header length


int print_vcf_gz_header(BGZF *fp, BGZF *bgzf_out, int print_header, kstring_t *tmp)
{
    char *buffer = (char*) fp->uncompressed_block;

    // Read the header and find the position of the data block
    if ( buffer[0]!='#' ) hts_log_error("Could not parse the header, expected '#', found '%c'\n", buffer[0]);

    int nskip = 1;     // end of the header in the current uncompressed block
    while (1)
    {
        if ( buffer[nskip]=='\n' )
        {
            nskip++;
            if ( nskip>=fp->block_length )
            {
                kputsn(buffer,nskip,tmp);
                if ( bgzf_read_block(fp) != 0 ) return -1;
                if ( !fp->block_length ) break;
                nskip = 0;
            }
            // The header has finished
            if ( buffer[nskip]!='#' )
            {
                kputsn(buffer,nskip,tmp);
                break;
            }
        }
        nskip++;
        if ( nskip>=fp->block_length )
        {
            kputsn(buffer,fp->block_length,tmp);
            if ( bgzf_read_block(fp) != 0 ) return -1;
            if ( !fp->block_length ) break;
            nskip = 0;
        }
    }
    if ( print_header )
    {
        if ( bgzf_write(bgzf_out,tmp->s,tmp->l) != tmp->l ) hts_log_error("Failed to write %"PRIu64" bytes\n", (uint64_t)tmp->l);
        tmp->l = 0;
    }
    return nskip;
}

inline int unpackInt16(const uint8_t *buffer)
{
    return buffer[0] | buffer[1] << 8;
}
int check_header(const uint8_t *header)
{
    if ( header[0] != 31 || header[1] != 139 || header[2] != 8 ) return -2;
    return ((header[3] & 4) != 0
            && unpackInt16((uint8_t*)&header[10]) == 6
            && header[12] == 'B' && header[13] == 'C'
            && unpackInt16((uint8_t*)&header[14]) == 2) ? 0 : -1;
}
void _check_hrecs(const bcf_hdr_t *hdr0, const bcf_hdr_t *hdr, const char *fname0, const char *fname)
{
    int j;
    for (j=0; j<hdr0->nhrec; j++)
    {
        bcf_hrec_t *hrec0 = hdr0->hrec[j];
        if ( hrec0->type!=BCF_HL_FLT && hrec0->type!=BCF_HL_INFO && hrec0->type!=BCF_HL_FMT && hrec0->type!=BCF_HL_CTG ) continue;    // skip fields w/o IDX
        int itag = bcf_hrec_find_key(hrec0, "ID");
        bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, hrec0->type, "ID", hrec0->vals[itag], NULL);

        char *type = NULL;
        if ( hrec0->type==BCF_HL_FLT ) type = "FILTER";
        if ( hrec0->type==BCF_HL_INFO ) type = "INFO";
        if ( hrec0->type==BCF_HL_FMT ) type = "FORMAT";
        if ( hrec0->type==BCF_HL_CTG ) type = "contig";

        if ( !hrec )
            hts_log_error("Cannot use --naive, incompatible headers, the tag %s/%s not present in %s\n",type,hrec0->vals[itag],fname);

        int idx0 = bcf_hrec_find_key(hrec0, "IDX");
        int idx  = bcf_hrec_find_key(hrec,  "IDX");
        if ( idx0<0 || idx<0 )
            hts_log_error("Unexpected IDX<0 for %s/%s in %s or %s\n",type,hrec0->vals[itag],fname0,fname);
        if ( strcmp(hrec0->vals[idx0],hrec->vals[idx]) )
            hts_log_error("Cannot use --naive, use --naive-force instead: different order the tag %s/%s in %s vs %s\n",type,hrec0->vals[itag],fname0,fname);
    }
}
void naive_concat_check_headers(int nfnames, const char **fnames)
{
    fprintf(stderr,"Checking the headers of %d files.\n", nfnames);
    bcf_hdr_t *hdr0 = NULL;
    int i,j;
    for (i=0; i< nfnames; i++)
    {
        htsFile *fp = hts_open(fnames[i], "r"); if ( !fp ) hts_log_error("Failed to open: %s\n", fnames[i]);
        bcf_hdr_t *hdr = bcf_hdr_read(fp); if ( !hdr ) hts_log_error("Failed to parse header: %s\n", fnames[i]);
        htsFormat type = *hts_get_format(fp);
        hts_close(fp);

        if ( i==0 )
        {
            hdr0 = hdr;
            continue;
        }

        // check the samples
        if ( bcf_hdr_nsamples(hdr0)!=bcf_hdr_nsamples(hdr) )
            hts_log_error("Cannot concatenate, different number of samples: %d vs %d in %s vs %s\n",bcf_hdr_nsamples(hdr0),bcf_hdr_nsamples(hdr),fnames[0],fnames[i]);
        for (j=0; j<bcf_hdr_nsamples(hdr0); j++)
            if ( strcmp(hdr0->samples[j],hdr->samples[j]) )
                hts_log_error("Cannot concatenate, different samples in %s vs %s\n",fnames[0],fnames[i]);

        // if BCF, check if tag IDs are consistent in the dictionary of strings
        if ( type.compression!=bgzf )
            hts_log_error("The --naive option works only for compressed BCFs or VCFs, sorry :-/\n");

        _check_hrecs(hdr0,hdr,fnames[0],fnames[i]);
        _check_hrecs(hdr,hdr0,fnames[i],fnames[0]);

        bcf_hdr_destroy(hdr);
    }
    if ( hdr0 ) bcf_hdr_destroy(hdr0);
    fprintf(stderr,"Done, the headers are compatible.\n");
}
void naive_concat(const char *output_fname, const int nfnames, const char **fnames)
{
    naive_concat_check_headers(nfnames, fnames);

    // only compressed BCF atm
    BGZF *bgzf_out = bgzf_open(output_fname,"w");;

    htsFormat output_type;
    output_type.format = bcf;
    output_type.compression = bgzf;

    const size_t page_size = BGZF_MAX_BLOCK_SIZE;
    uint8_t *buf = (uint8_t*) malloc(page_size);
    kstring_t tmp = {0,0,0};
    int i, file_types = 0;
    for (i=0; i<nfnames; i++)
    {
        htsFile *hts_fp = hts_open(fnames[i],"r");
        if ( !hts_fp ) hts_log_error("\nFailed to open: %s\n", fnames[i]);
        htsFormat type = *hts_get_format(hts_fp);

        file_types |= type.format==vcf ? 1 : 2;

        BGZF *fp = hts_get_bgzfp(hts_fp);
        if ( !fp || bgzf_read_block(fp) != 0 || !fp->block_length )
            hts_log_error("\nFailed to read %s: %s\n", fnames[i], strerror(errno));

        int nskip;
        if ( type.format==bcf )
        {
            const size_t magic_len = 5 + 4; // "Magic" string + header length
            uint8_t magic[MAGIC_LEN];
            if ( bgzf_read(fp, magic, magic_len) != magic_len ) hts_log_error("\nFailed to read the BCF header in %s\n", fnames[i]);
            // First five bytes are the "Magic" string
            if (strncmp((char*)magic, "BCF\2\2", 5) != 0) hts_log_error("\nInvalid BCF magic string in %s\n", fnames[i]);
            // Next four are the header length (little-endian)
            tmp.l = le_to_u32(magic + 5);
            hts_expand(char,tmp.l,tmp.m,tmp.s);
            if ( bgzf_read(fp, tmp.s, tmp.l) != tmp.l ) hts_log_error("\nFailed to read the BCF header in %s\n", fnames[i]);

            // write only the first header
            if ( i==0 )
            {
                if ( bgzf_write(bgzf_out, magic, magic_len) != magic_len ) hts_log_error("\nFailed to write %zu bytes to %s\n", magic_len, output_fname);
                if ( bgzf_write(bgzf_out, tmp.s, tmp.l) != tmp.l) hts_log_error("\nFailed to write %"PRId64" bytes to %s\n", (uint64_t)tmp.l, output_fname);
            }
            nskip = fp->block_offset;
        }
        else
        {
            nskip = print_vcf_gz_header(fp, bgzf_out, i==0?1:0, &tmp);
            if ( nskip==-1 ) hts_log_error("\nhts_log_error reading %s\n", fnames[i]);
        }

        // Output all non-header data that were read together with the header block
        if ( fp->block_length - nskip > 0 )
        {
            if ( bgzf_write(bgzf_out, (char *)fp->uncompressed_block+nskip, fp->block_length-nskip)<0 ) hts_log_error("\nhts_log_error: %d\n",fp->errcode);
        }
        if ( bgzf_flush(bgzf_out)<0 ) hts_log_error("\nhts_log_error: %d\n",bgzf_out->errcode);


        // Stream the rest of the file as it is, without recompressing, but remove BGZF EOF blocks
        // The final bgzf eof block will be added by bgzf_close.
        ssize_t nread, nblock, nwr;
        const int nheader = 18, neof = 28;
        const uint8_t *eof = (uint8_t*) "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0";
        while (1)
        {
            nread = bgzf_raw_read(fp, buf, nheader);
            if ( !nread ) break;
            if ( nread != nheader || check_header(buf)!=0 ) hts_log_error("\nCould not parse the header of a bgzf block: %s\n", fnames[i]);
            nblock = unpackInt16(buf+16) + 1;
            assert( nblock <= page_size && nblock >= nheader );
            nread += bgzf_raw_read(fp, buf+nheader, nblock - nheader);
            if ( nread!=nblock ) hts_log_error("\nCould not read %"PRId64" bytes: %s\n",(uint64_t)nblock, fnames[i]);
            if ( nread==neof && !memcmp(buf,eof,neof) ) continue;
            nwr = bgzf_raw_write(bgzf_out, buf, nread);
            if ( nwr != nread ) hts_log_error("\nWrite failed, wrote %"PRId64" instead of %d bytes.\n", (uint64_t)nwr,(int)nread);
        }
        if (hts_close(hts_fp)) hts_log_error("\nClose failed: %s\n", fnames[i]);
    }
    free(buf);
    free(tmp.s);
    if (bgzf_close(bgzf_out) < 0) hts_log_error("hts_log_error: %d\n",bgzf_out->errcode);
}
