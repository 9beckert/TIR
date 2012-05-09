/*
  Copyright (c) 2012 Manuela Beckert <9beckert@informatik.uni-hamburg.de>
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

/* A stream, which finds Terminal Inverted Repeats in the encoded sequence
   and returns a pair of them with every next-call. */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "core/arraydef.h"
#include "core/encseq.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/str_api.h"
#include "core/undef_api.h"
#include "extended/genome_node.h"
#include "extended/node_stream_api.h"
#include "extended/region_node_api.h"
#include "ltr/ltr_xdrop.h"
#include "match/esa-maxpairs.h"
#include "match/esa-seqread.h"
#include "match/greedyedist.h"
#include "match/spacedef.h"
#include "tools/gt_tir_stream.h"

/* A pair of Seeds (identical inverted parts in the sequence) */
typedef struct
{
  unsigned long pos1;         /* position of first seed */
  unsigned long pos2;         /* position of seconf seed (other contig) */
  unsigned long offset;       /* distance between them related to the actual sequence 
				 (not mirrored) */
  unsigned long len;          /* length of the seed  */
  unsigned long contignumber; /* number of contig for this seed */
} Seed;

GT_DECLAREARRAYSTRUCT(Seed);

/* A pair of TIRs (extended Seeds) */
typedef struct
{
  unsigned long contignumber,
                left_tir_start,  /* first position of TIR on forward strand */ 
                left_tir_end,    /* last position of TIR on forward strand */    
                right_tir_start, /* first position of TIR on reverse strand */  
                right_tir_end;   /* last position of TIR on reverse strand */   
  double        similarity;      /* similarity of the two TIRs */
} TIRPair;

GT_DECLAREARRAYSTRUCT(TIRPair);

/* currently not really in use */
typedef enum {
  GT_TIR_STREAM_STATE_START,
  GT_TIR_STREAM_STATE_REGIONS,
  GT_TIR_STREAM_STATE_COMMENTS,
  GT_TIR_STREAM_STATE_FEATURES
} GtTIRStreamState;

/* The arguments */
struct GtTIRStream
{
  const GtNodeStream parent_instance; /* node which could be before us to call next on */
  GtStr *str_indexname;               /* index name of suffix array */
  unsigned long minseedlength;        /* minimal seed length */
  const GtEncseq *encseq;             /* encoded sequence */
  Sequentialsuffixarrayreader *ssar;  /* suffix array reader */
  GtArraySeed seedarray;              /* array containing our seeds */
  GtTIRStreamState state;             /* current state of output */
  GtArrayTIRPair tir_pairs;           /* array contianing our TIRs */
  Arbitraryscores arbitscores;        /* cost of match, mismatch, insertion, deletion */
  int xdrop_belowscore;               /* threshold for xdop algorithm when to discard extensions */
  double similarity_threshold;        /* decides whether a TIR is accepted or not */
  unsigned long num_of_tirs;          /* number of TIRs found */
  unsigned long cur_elem_index;       /* index of the TIR to be put out next */
  unsigned long prev_seqnum;          /* number of previous contig */
};

/*
 * Stores the seeds in an array.
 * (processmaxpairs call-back-function)
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
    
  }
  return 0;
}

static int gt_compare_TIRs(TIRPair *pair1, TIRPair *pair2)
{
  if(pair1->contignumber < pair2->contignumber)
  {
    return -1;
  }
  else if(pair1->contignumber == pair2->contignumber)
  {
    if(pair1->left_tir_start < pair2->left_tir_start)
    {
      return -1;
    }
    else if(pair1->left_tir_start == pair2-> left_tir_start)
    {
      if(pair1->right_tir_start < pair2->right_tir_start)
      {
        return -1;
      }
      else if(pair1->right_tir_start == pair2->right_tir_start)
      {
        return 0;
      }
    }
  }
  return 1;
}

/* swaps the elements at pos1 and pos2 in tir_pairs */
static void gt_swap_TIRs(GtArrayTIRPair *tir_pairs, unsigned long pos1, unsigned long pos2)
{
  TIRPair tmp;
  tmp = tir_pairs->spaceTIRPair[pos1];
  tir_pairs->spaceTIRPair[pos1] = tir_pairs->spaceTIRPair[pos2];
  tir_pairs->spaceTIRPair[pos2] = tmp;
}

/* sorts an array of tirs with inplace quicksort*/
static void gt_sort_TIRs(GtArrayTIRPair *tir_pairs,unsigned long start, unsigned long end)
{
  TIRPair pivot;
  unsigned long l;
  unsigned long r;
  
  if(end > start)
  {
    pivot = tir_pairs->spaceTIRPair[start];
    l = start + 1;
    r = end;
    while(l < r)
    {
      if(gt_compare_TIRs(&tir_pairs->spaceTIRPair[l], &pivot) <= 0)
      {
        l++;
      }
      else
      {
        r--;
        gt_swap_TIRs(tir_pairs, l, r);
      }
    }
    l--;
    gt_swap_TIRs(tir_pairs, start, l);
    gt_sort_TIRs(tir_pairs, start, l);
    gt_sort_TIRs(tir_pairs, r, end);
    
  }
  
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
  int had_err = 0;
  
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
    
    /* We want the actual positions (not mirrored) */
    unsigned long right_tir_start = seedptr->pos2 + seedptr->len - 1 + xdropbest_right.jvalue; // end of mirrored is start of actual TIR
    unsigned long right_tir_end   = seedptr->pos2 - xdropbest_left.jvalue;  // start of mirrored is end of actual TIR
    
    /* We can get the corresponding position of a mirrored one with GT_REVERSEPOS(total length,position) */
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
                            
    /* Discard this candidate if similarity too small or increase number of TIRs*/
    if(gt_double_smaller_double(pair->similarity, 
        tir_stream->similarity_threshold))                                              
    {
      tir_pairs->nextfreeTIRPair--;

    }
    else
    {
      tir_stream->num_of_tirs++;
    }      

  }
  
  FREESPACE(left_tir_char);
  FREESPACE(right_tir_char);
  
  /* sort the tir_pairs */
  gt_sort_TIRs(tir_pairs, 0, tir_stream->num_of_tirs);
  
  return had_err;
  
}


/*
 * Saves the next node of the annotation graph in gn.
 */
static int gt_tir_stream_next(GtNodeStream *ns,
                              GT_UNUSED GtGenomeNode **gn,
                              GtError *err)
{
  GtTIRStream *tir_stream;
  int had_err = 0;
  //int count;
  gt_error_check(err);
  tir_stream = gt_node_stream_cast(gt_tir_stream_class(), ns);
  
  /* generate and check seeds */
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
    
    /* extend seeds to TIRs and check TIRs */
    if(!had_err && gt_searchforTIRs(tir_stream, &tir_stream->tir_pairs, tir_stream->encseq, err) != 0)
    {
      had_err = -1;
    }

    /* free the seed array since we don't need it any longer */
    GT_FREEARRAY(&tir_stream->seedarray, Seed);
    
    // TODO apply further filters like removing duplicates, sort elements!!!!!!

    tir_stream->state = GT_TIR_STREAM_STATE_REGIONS;
  }
  
  /* output on console */
  /*printf("TIRs found:\n");
  for(count = 0; count < tir_stream->tir_pairs.nextfreeTIRPair; count++)
  {
    printf("contig %lu\n left: start %lu, end %lu\n right: start %lu, end %lu\n similarity: %f\n\n", 
      tir_stream->tir_pairs.spaceTIRPair[count].contignumber,
      tir_stream->tir_pairs.spaceTIRPair[count].left_tir_start,
      tir_stream->tir_pairs.spaceTIRPair[count].left_tir_end,
      tir_stream->tir_pairs.spaceTIRPair[count].right_tir_start,
      tir_stream->tir_pairs.spaceTIRPair[count].right_tir_end,
      tir_stream->tir_pairs.spaceTIRPair[count].similarity);
  }*/

  /* stream out the region nodes */
  if (!had_err && tir_stream->state == GT_TIR_STREAM_STATE_REGIONS) 
  {

    bool skip = false;

    /* check whether index is valid */
    if(tir_stream->cur_elem_index < tir_stream->num_of_tirs)
    {
      unsigned long seqnum, seqlength;
      GtGenomeNode *rn;
      GtStr *seqid;
      seqnum = tir_stream->tir_pairs.spaceTIRPair[tir_stream->cur_elem_index].contignumber;

      /* if first time we do this */
      if(tir_stream->prev_seqnum == GT_UNDEF_ULONG)
      {
        /* use current seqnum */
        tir_stream->prev_seqnum = seqnum;
      }
      else
      {
        /* else get seqnum of next contig */
        while(tir_stream->prev_seqnum == seqnum)
        {
          tir_stream->cur_elem_index++;
          
          /* don't go on if index is out of bounds */
          if(tir_stream->cur_elem_index >= tir_stream->num_of_tirs)
          {
            skip = true;
            break;
          }
          
          seqnum = tir_stream->tir_pairs.spaceTIRPair[tir_stream->cur_elem_index].contignumber;
        }
      }
      
      /* create node */
      if(!skip)
      {
        tir_stream->prev_seqnum = seqnum; /* for next call of REGIONS */
        seqlength = gt_encseq_seqlength(tir_stream->encseq, seqnum);
        
        /* add description */
        seqid = gt_str_new_cstr("seq");
        gt_str_append_ulong(seqid, seqnum);
        rn = gt_region_node_new(seqid, 1, seqlength); // seqid, start, end
        
        gt_str_delete(seqid);
        *gn = rn; // set return variable
        tir_stream->cur_elem_index++;
      }
      else
      {
        /* we skip */
        tir_stream->cur_elem_index = 0;
        tir_stream->state = GT_TIR_STREAM_STATE_COMMENTS;
        *gn = NULL;
      }
    }
    else
    {
      /* no valid index */
      tir_stream->cur_elem_index = 0;
      tir_stream->state = GT_TIR_STREAM_STATE_COMMENTS;
      *gn = NULL;
    }
  }

  /* then stream out the comment nodes */
  if (!had_err && tir_stream->state == GT_TIR_STREAM_STATE_COMMENTS) 
  {
    bool skip = false;
    if (tir_stream->cur_elem_index < tir_stream->num_of_tirs)
    {
      const char *description;
      char description_string[BUFSIZ];
      unsigned long description_len, seqnum;
      GtGenomeNode *cn;
      
      seqnum = tir_stream->tir_pairs.spaceTIRPair[tir_stream->cur_elem_index].contignumber;
      
      /* for the first time */
      if(tir_stream->prev_seqnum == GT_UNDEF_LONG)
      {
        /* use current seqnum */
        tir_stream->prev_seqnum = seqnum;
      }
      else
      {
        /* else get seqnum of next contig */
        while(tir_stream->prev_seqnum == seqnum)
        {
          tir_stream->cur_elem_index++;
          if(tir_stream->cur_elem_index >= tir_stream->num_of_tirs)
          {
            skip = true;
            break;
          }
        }
        
        seqnum = tir_stream->tir_pairs.spaceTIRPair[tir_stream->cur_elem_index].contignumber;
      }
      
      /* create a new comment node */
      if(!skip)
      {
        tir_stream->prev_seqnum = seqnum;
        
        /* get description and descriptionlength of current contig */
        description = gt_encseq_description(tir_stream->encseq, &description_len, seqnum);
       
       /* make a \0 terminated string of description */ 
      (void) strncpy(description_string, description, (size_t) (description_len * sizeof(char)));
      description_string[description_len] = '\0';
      
      /* make a new comment node */
      cn = gt_comment_node_new(description_string);
      
      *gn = cn;
      tir_stream->cur_elem_index++;
      }
      else
      {
        /* skip */
        tir_stream->cur_elem_index = 0;
        tir_stream->state = GT_TIR_STREAM_STATE_FEATURES;
        *gn = NULL;
      }
    }
    else
    {
      /* no valid index */
      tir_stream->cur_elem_index = 0;
      tir_stream->state = GT_TIR_STREAM_STATE_FEATURES;
      *gn = NULL;
    }
  }

  /* finally stream out the features */
  if (!had_err && tir_stream->state == GT_TIR_STREAM_STATE_FEATURES)
  {
    if(tir_stream->cur_elem_index < tir_stream->num_of_tirs)
    {
      GtStr *seqid, *source;
      GtGenomeNode *node, *parent;
      const TIRPair *pair = &tir_stream->tir_pairs.spaceTIRPair[tir_stream->cur_elem_index];
      unsigned long seqstartpos;
      
      seqstartpos = gt_encseq_seqstartpos(tir_stream->encseq, pair->contignumber);
      seqid = gt_str_new_cstr("seq");
      source = gt_str_new_cstr("TIR");
      
      gt_str_append_ulong(seqid, pair->contignumber);
      
      /* features */
      /* repeat region */
      /* make the parent node */
      node = gt_feature_node_new(seqid, "repeat region", pair->left_tir_start - seqstartpos + 1, pair->right_tir_start - seqstartpos + 1, GT_STRAND_UNKNOWN);
      
      /* set stuff */
      //TODO woher kommt gt_feature_nodes_set_source funktion?
      //gt_feature_nodes_set_source((GtFeatureNode*) node, source);
      *gn = node;
      parent = node;
      
      /*  */
      
      //TODO alle anderen features müssen noch implementiert 
      
      /* clean up and get next pair */
      gt_str_delete(seqid);
      gt_str_delete(source);
      tir_stream->cur_elem_index++;
    }
    else
    {
      *gn = NULL;
    }
  } 
  
  return had_err;
}

/*
 * Frees stuff.
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
 * This function is needed to create a GtNodeStream with the function gt_node_stream_create(class, bool)
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
 * Creates a new TIR-stream
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
  tir_stream->str_indexname = gt_str_ref(str_indexname);  /* use of gt_str_ref so the reference counter will be increased */
  tir_stream->minseedlength = minseedlength;
  tir_stream->arbitscores = arbitscores;
  tir_stream->xdrop_belowscore = xdrop_belowscore;
  tir_stream->similarity_threshold = similarity_threshold;
  tir_stream->num_of_tirs = 0;
  tir_stream->cur_elem_index = 0;
  tir_stream->prev_seqnum = GT_UNDEF_ULONG;
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
