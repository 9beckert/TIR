/*
  Copyright (c) 2012 Dorle Osterode <9osterod@informatik.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "core/arraydef.h"
#include "core/encseq.h"
#include "core/str_api.h"
#include "extended/genome_node.h"
#include "extended/node_stream_api.h"
#include "match/esa-maxpairs.h"
#include "match/esa-seqread.h"


typedef struct
{
  unsigned long pos1;         /* first position of maximal repeat (seed) */
  unsigned long offset;       /* second position = pos1 + offset */
  unsigned long len;          /* length of maximal repeat  */
  unsigned long contignumber; /* number of contig for this repeat */
} Seed;

GT_DECLAREARRAYSTRUCT(Seed);

typedef enum {
  GT_TIR_STREAM_STATE_START,
  GT_TIR_STREAM_STATE_REGIONS,
  GT_TIR_STREAM_STATE_COMMENTS,
  GT_TIR_STREAM_STATE_FEATURES
} GtTIRStreamState;

struct GtTIRStream
{
  const GtNodeStream parent_instance;
  GtStr *str_indexname;
  unsigned long minseedlength;
  const GtEncseq *encseq;
  Sequentialsuffixarrayreader *ssar;
  GtArraySeed seedarray;
  GtTIRStreamState state;
};

/*
 * processmaxpairs call-back-function, which stores the seeds in an array
 */
static int gt_store_seeds(void *info,
                          const GtEncseq *encseq,
                          unsigned long len,
                          unsigned long pos1,
                          unsigned long pos2,
                          GT_UNUSED GtError *err)
{
  unsigned long offset;
  unsigned long length_encseq;
  unsigned long contignumber = 0,
                seqnum1,
                seqnum2;
  bool samecontig = false;
  bool diffstrands = false;
  GtEncseq *encseq = encseq;
  GtArraySeed *seeds = (GTArraySeed *) info;

  gt_error_check(err);
  
  /* different strands */
  length_encseq = gt_encseq_seqlength(encseq, contignumber);
  if(pos1 <= length_encseq/2 && pos2 > length_encseq/2)
  {
    diffstrands = true;
  }
  
  /* same contig */
  seqnum1 = gt_encseq_seqnum(encseq, pos1);
  seqnum2 = gt_encseq_seqnum(encseq, pos2);
  
  if(seqnum1 == seqnum2)
  {
    samecontig = true;
    contignumber = seqnum1;
  }
  
  /* calculate offset */
  offset = length_encseq - pos2 - pos1 + 1;
  
  /* like this, cause we want to have other constraints later */
  if (samecontig && diffstrands)
  {
    Seed *nextfreeseedpointer;
    
    //TODO was bedeutet die 10??
    GT_GETNEXTFREEINARRAY(nextfreeseedpointer, &seeds,
                       Seed, 10);
    
    nextfreeseedpointer->pos1 = pos1;
    nextfreeseedpointer->offset = offset;
    nextfreeseedpointer->len = len;
    nextfreeseedpointer->contignumber = contignumber;
  }
  return 0;
}


/*
 * saves the next node of the annotation graph in gn
 */
static int gt_tir_stream_next(GT_UNUSED GtNodeStream *ns,
                              GtGenomeNode **gn,
                              GtError *err)
{
  // TODO magic has to be understood
  GtTIRStream *tir_stream;
  int had_err = 0;
  int seedcount;
  gt_error_check(err);
  tir_stream = gt_node_stream_cast(gt_tir_stream_class(), ns);
  
  /* find seeds */
  if(tir_stream->state == GT_TIR_STREAM_STATE_START)
  {
    GT_INITARRAY(&tir_stream->seedarray, Seed);
    
    if (!had_err && gt_enumeratemaxpairs(tir_stream->ssar,
                      tir_stream->encseq,
                      gt_readmodeSequentialsuffixarrayreader(tir_stream->ssar),
                      (unsigned int) tir_stream->minseedlength,
                      gt_store_seeds,
                      &tir_stream->seedarray,
                      err) != 0)
    {
      had_err = -1;
    }
    //tir_stream->state = GT_TIR_STREAM_STATE_REGIONS;
  }
  
  /* output */
  
  printf("Seeds found:\n");
  for(seedcount = 0; seedcount < tir_stream->seedarray.nextfreeSeed; seedcount++)
  {
    printf("contig %lu, position %lu, length %lu, offset %lu\n", 
      tir_stream->seedarray.spaceSeed[seedcount].contignumber,
      tir_stream->seedarray.spaceSeed[seedcount].pos1,
      tir_stream->seedarray.spaceSeed[seedcount].len,
      tir_stream->seedarray.spaceSeed[seedcount].offset);
  }
  
  /* magic has to happen here */
  
}

/*
 * frees stuff
 */
static void gt_tir_stream_free(GtNodeStream *ns)
{
  GtTIRStream *tir_stream = gt_node_stream_cast(gt_tir_stream_class() , ns);
  
  if (tir_stream->ssar != NULL)
    gt_freeSequentialsuffixarrayreader(&tir_stream->ssar);
}

/*
 * this function is needed to create a GtNodeStream with the function gt_node_stream_create(class, bool)
 */
const GtNodeStreamClass* gt_tir_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtTIRStream),
                                   gt_tir_stream_free,
                                   gt_tir_stream_next);
  }
  return nsc;
}

/*
 * creates a new tir-stream
 */
GtNodeStream* gt_tir_stream_new(GtStr *str_indexname,
                                unsigned long minseedlength,
                                GtError *err)
{
  int had_err = 0;
  GtNodeStream *ns = gt_node_stream_create(gt_tir_stream_class(), false);
  GtTIRStream *tir_stream = gt_node_stream_cast(gt_tir_stream_class(), ns);
  tir_stream->str_indexname = str_indexname;
  tir_stream->minseedlength = minseedlength;
  tir_stream->state = GT_TIR_STREAM_STATE_START;
  tir_stream->ssar = 
      gt_newSequentialsuffixarrayreaderfromfile(gt_str_get(str_indexname),
                                                SARR_LCPTAB | SARR_SUFTAB |
                                                SARR_ESQTAB | SARR_DESTAB |
                                                SARR_SSPTAB | SARR_SDSTAB,
                                                SEQ_mappedboth,
                                                NULL,
                                                err);
  if (tir_stream->ssar == NULL)
  {
    gt_node_stream_delete(ns);
    return NULL;
  }
  
  tir_stream->encseq = gt_encseqSequentialsuffixarrayreader(tir_stream->ssar);
  
  if(!had_err)
  {
    return ns;
  }
}
