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
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/str_api.h"
#include "extended/genome_node.h"
#include "extended/node_stream_api.h"
#include "ltr/ltr_xdrop.h"
#include "match/esa-maxpairs.h"
#include "match/esa-seqread.h"
#include "match/greedyedist.h"
#include "match/spacedef.h"
#include "tools/gt_tir_stream.h"


typedef struct
{
  unsigned long pos1;         /* position of first seed */
  unsigned long pos2;         /* position of seconf seed (other contig) */
  unsigned long offset;       /* distance between them IRL */
  unsigned long len;          /* length of maximal seed  */
  unsigned long contignumber; /* number of contig for this seed */
} Seed;

GT_DECLAREARRAYSTRUCT(Seed);

typedef struct
{
  unsigned long contignumber,
                left_tir_start,  /* first position of TIR on forward strand */ 
                left_tir_end,    /* last position of TIR on forward strand */   // TODO sollten sich auf den nicht gespiegelten teil  beziehen!!!
                right_tir_start,  /* first position of TIR on reverse strand */  // dafür gibts ein makro: GT_REVERSEPOS(gesamtlänge,position)
                right_tir_end;    /* last position of TIR on reverse strand */
  double        similarity;         /* similarity of the two TIRs */
} TIRPair;

GT_DECLAREARRAYSTRUCT(TIRPair);


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
  GtArrayTIRPair tir_pairs;
  Arbitraryscores arbitscores;
  int xdrop_belowscore;
  double similarity_threshold;
};

/*
 * processmaxpairs call-back-function, which stores the seeds in an array
 */
static int gt_store_seeds(void *info,
                          const GtEncseq *encseq,
                          unsigned long len,
                          unsigned long pos1,
                          unsigned long pos2,
                          GtError *err)
{
  unsigned long offset = 0;
  unsigned long contignumber = 0,
                seqnum1,
                seqnum2;
  unsigned long num_of_contigs;
  unsigned long contig_length;
  unsigned long offset_pos1;
  unsigned long offset_pos2;
  unsigned long startpos_contig1;
  unsigned long startpos_contig2;
  unsigned long tmp = 0;
  bool samecontig_and_diffstrands = false;
  bool direction = false;
  GtArraySeed *seeds = (GtArraySeed *) info;

  gt_error_check(err);
  
  /* check if pos1 is lower pos2 */
  if(pos1 > pos2)
  {
    tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
  
  /* same contig and different strands */
  num_of_contigs = gt_encseq_num_of_sequences(encseq);
  
  seqnum1 = gt_encseq_seqnum(encseq, pos1);
  seqnum2 = gt_encseq_seqnum(encseq, pos2);
  
  if(seqnum2 == num_of_contigs - seqnum1 -1)
  {
    contignumber = seqnum1;
    samecontig_and_diffstrands = true;
  }
  
  /* offset is the distance between the first position of the seed on the 
  5'3'-strand and the first position of the seed on the 3'5'-strand */
  
  contig_length = gt_encseq_seqlength(encseq, seqnum1);
  startpos_contig1 = gt_encseq_seqstartpos(encseq, seqnum1);
  startpos_contig2 = gt_encseq_seqstartpos(encseq, seqnum2);
  
  offset_pos1 = pos1 - startpos_contig1;
  offset_pos2 = pos2 - startpos_contig2;
  offset_pos2 = offset_pos2 + len -1;
  offset_pos2 = contig_length - offset_pos2 -1;
  
  if(offset_pos1 < offset_pos2)
  {
    offset = offset_pos2 - offset_pos1;
    direction = true;
  }
  

  
  /* like this, cause we want to have other constraints later */
  if (samecontig_and_diffstrands && direction)
  {
    Seed *nextfreeseedpointer;
    
    GT_GETNEXTFREEINARRAY(nextfreeseedpointer, seeds,
                       Seed, 10);
    
    nextfreeseedpointer->pos1 = pos1;
    nextfreeseedpointer->pos2 = pos2;
    nextfreeseedpointer->offset = offset;
    nextfreeseedpointer->len = len;
    nextfreeseedpointer->contignumber = contignumber;
    
    printf("Seeds found:\n");
    printf("pos1: %lu\n", pos1);
    printf("pos2: %lu\n", pos2);
    printf("len: %lu\n", len);
    printf("offset: %lu\n\n", offset);
  }
  return 0;
}

/*
 * Extends the Seeds and searches for TIRs.
 */
static int gt_searchforTIRs(GtTIRStream *tir_stream,
                             GtArrayTIRPair *tir_pairs,
                             const GtEncseq *encseq,
                             GtError *err)
{
  unsigned long seedcounter = 0;
  GtArrayMyfrontvalue fronts; /* needed to use xdrop */
  Myxdropbest xdropbest_left; /* return parameters for xdrop */
  Myxdropbest xdropbest_right;
  //unsigned long alignment_length = 0,
  unsigned long total_length = 0,
                left_tir_length = 0, 
                right_tir_length = 0,
                max_left_length = 0,    /* maximal length of left TIR*/ 
                max_right_length = 0;
  GtUchar *left_tir_char = NULL,        /* next character to align */
          *right_tir_char = NULL;
  unsigned long edist = 0;              /* edit distance */
  Seed *seedptr;
  TIRPair *pair;
  bool had_err = false;
  
  /* Did we have errors so far? */
  gt_error_check(err);
  
  /* Iterating over seeds */
  for(seedcounter = 0; seedcounter < tir_stream->seedarray.nextfreeSeed;
      seedcounter++)
  {
    /* Getting the seed entry at position seedcounter */
    seedptr = &(tir_stream->seedarray.spaceSeed[seedcounter]);
    
    /* Initialize array for extending seeds */
    GT_INITARRAY (&fronts, Myfrontvalue);
    
    /* Left alignment */
    gt_evalxdroparbitscoresleft(&tir_stream->arbitscores,
                               &xdropbest_left,
                               &fronts,
                               encseq,
                               encseq,
                               seedptr->pos1,
                               seedptr->pos2,
                               (int) seedptr->pos1,
                               (int) seedptr->pos2,
                               (Xdropscore)tir_stream->xdrop_belowscore);
                               
    /* Re-Initialize fronts */
    GT_FREEARRAY (&fronts, Myfrontvalue);
    GT_INITARRAY (&fronts, Myfrontvalue);

    
    
    /* init total_length */
    total_length = gt_encseq_total_length(encseq);
    
    /* Right alignment */
    gt_evalxdroparbitscoresright(&tir_stream->arbitscores,
                                &xdropbest_right,
                                &fronts,
                                encseq,
                                encseq,
                                seedptr->pos1 + seedptr->len,
                                seedptr->pos2 + seedptr->len,
                                (int) (total_length - (seedptr->pos1 + seedptr->len)),
                                (int) (total_length - (seedptr->pos2 + seedptr->len)),
                                (Xdropscore)tir_stream->xdrop_belowscore);
    GT_FREEARRAY (&fronts, Myfrontvalue);
    
    GT_GETNEXTFREEINARRAY(pair, tir_pairs, TIRPair, 5);
    
    /* Store positions for the found TIR */
    pair->contignumber = seedptr->contignumber;
    pair->left_tir_start = seedptr->pos1 - xdropbest_left.ivalue;
    pair->left_tir_end = seedptr->pos1 + seedptr->len - 1 + xdropbest_right.ivalue;
    
    // Wir wollen die tatsächlichen Positionen (dh nicht auf mirrored bezogen)
    unsigned long right_tir_start = seedptr->pos2 + seedptr->len - 1 + xdropbest_right.jvalue; // end bezogen auf mirrored
    unsigned long right_tir_end   = seedptr->pos2 - xdropbest_left.jvalue;  // start bezogen auf mirrored
    
    // Dafür benutzen wir GT_REVERSEPOS(gesamtlänge,position)
    pair->right_tir_start = GT_REVERSEPOS(total_length,right_tir_start);
    pair->right_tir_end = GT_REVERSEPOS(total_length,right_tir_end); 
    pair->similarity = 0.0;
    
    /* Check similarity */
    left_tir_length = pair->left_tir_end - pair->left_tir_start + 1;
    right_tir_length = pair->right_tir_end - pair->right_tir_start + 1;
    
    /* Realloc tir_char if length too high */
    if(left_tir_length > max_left_length)
    {
      max_left_length = left_tir_length;
      ALLOCASSIGNSPACE(left_tir_char, left_tir_char, GtUchar, max_left_length);
    }
    if(right_tir_length > max_right_length)
    {
      max_right_length = right_tir_length;
      ALLOCASSIGNSPACE(right_tir_char, right_tir_char, GtUchar, max_right_length);
    }
    
    /* Store encoded substring */
    gt_encseq_extract_encoded(encseq, left_tir_char,
                                         pair->left_tir_start,
                                         pair->left_tir_end);
    gt_encseq_extract_encoded(encseq, right_tir_char,
                                         right_tir_end,
                                         right_tir_start);
    
    /* Get edit distance */
    edist = greedyunitedist(left_tir_char,(unsigned long) left_tir_length, 
                            right_tir_char,(unsigned long) right_tir_length);
                            
    /* Determine similarity */
    pair->similarity = 100.0 * (1 - (((double) edist)/
                      (MAX(left_tir_length, right_tir_length))));
                            
    /* Discard this candidate if similarity too small */
    if(gt_double_smaller_double(pair->similarity, 
        tir_stream->similarity_threshold))                                              
    {
      tir_pairs->nextfreeTIRPair--;
    }      

  }
  
  FREESPACE(left_tir_char);
  FREESPACE(right_tir_char);
  
  // TODO schaun ob array zu free'n ist
  
  return had_err;
  
}


/*
 * saves the next node of the annotation graph in gn
 */
static int gt_tir_stream_next(GtNodeStream *ns,
                              GT_UNUSED GtGenomeNode **gn,
                              GtError *err)
{
  GtTIRStream *tir_stream;
  int had_err = 0;
  int count;
  gt_error_check(err);
  tir_stream = gt_node_stream_cast(gt_tir_stream_class(), ns);
  
  /* find seeds */
  if(tir_stream->state == GT_TIR_STREAM_STATE_START)
  {
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
    
    if(!had_err && gt_searchforTIRs(tir_stream, &tir_stream->tir_pairs, tir_stream->encseq, err) != 0)
    {
      had_err = -1;
    }
    
    //tir_stream->state = GT_TIR_STREAM_STATE_REGIONS;
  }
  
  /* output */
  
  printf("TIRs found:\n");
  for(count = 0; count < tir_stream->tir_pairs.nextfreeTIRPair; count++)
  {
    printf("contig %lu\n left: start %lu, end %lu\n right: start %lu, end %lu\n similarity: %f\n\n", 
      tir_stream->tir_pairs.spaceTIRPair[count].contignumber,
      tir_stream->tir_pairs.spaceTIRPair[count].left_tir_start,
      tir_stream->tir_pairs.spaceTIRPair[count].left_tir_end,
      tir_stream->tir_pairs.spaceTIRPair[count].right_tir_start,
      tir_stream->tir_pairs.spaceTIRPair[count].right_tir_end,
      tir_stream->tir_pairs.spaceTIRPair[count].similarity);
  }
  
  /* magic has to happen here */
  
  return 0;
}

/*
 * frees stuff
 */
static void gt_tir_stream_free(GtNodeStream *ns)
{
  GtTIRStream *tir_stream = gt_node_stream_cast(gt_tir_stream_class(), ns);
  
  // gt_str_delete damit Referenzcounter erniedrigt wird (str_indexname)
  
  gt_str_delete(tir_stream->str_indexname);
  
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
                                Arbitraryscores arbitscores,
                                int xdrop_belowscore,
                                double similarity_threshold,
                                GtError *err)
{
  int had_err = 0;
  GtNodeStream *ns = gt_node_stream_create(gt_tir_stream_class(), false);
  GtTIRStream *tir_stream = gt_node_stream_cast(gt_tir_stream_class(), ns);
  tir_stream->str_indexname = gt_str_ref(str_indexname);  // gt_str_ref damit auch der Referenzcounter erhoeht wird
  tir_stream->minseedlength = minseedlength;
  tir_stream->arbitscores = arbitscores;
  tir_stream->xdrop_belowscore = xdrop_belowscore;
  tir_stream->similarity_threshold = similarity_threshold;
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
  GT_INITARRAY(&tir_stream->seedarray, Seed);
  GT_INITARRAY(&tir_stream->tir_pairs, TIRPair);
  
  
  if(!had_err)
  {
    return ns;
  }
  
  // TODO ueberlegen, ob had_err benutzt werden müsste
  return NULL;
}
