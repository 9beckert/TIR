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

#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "core/arraydef.h"
#include "core/encseq.h"
#include "core/log_api.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/str_api.h"
#include "core/undef_api.h"
#include "extended/genome_node.h"
#include "extended/node_stream_api.h"
#include "extended/region_node_api.h"
#include "ltr/ltr_xdrop.h"
#include "match/esa-maxpairs.h"
#include "match/esa-mmsearch.h"
#include "match/esa-seqread.h"
#include "match/greedyedist.h"
#include "match/spacedef.h"
#include "match/querymatch.h"
#include "tools/gt_tir_stream.h"

/* A pair of Seeds (identical inverted parts in the sequence) */
typedef struct
{
  unsigned long pos1;         /* position of first seed */
  unsigned long pos2;         /* position of seconf seed (other contig) */
  unsigned long offset;       /* distance between them related to the actual sequence (not mirrored) */
  unsigned long len;          /* length of the seed  */
  unsigned long contignumber; /* number of contig for this seed */
} Seed;
GT_DECLAREARRAYSTRUCT(Seed);

typedef struct
{
	GtArraySeed seed;
	unsigned long max_tir_length;
	unsigned long min_tir_length;
	unsigned long max_tir_distance;
	unsigned long min_tir_distance;
} SeedInfo;												/* this is needed, to check distance constraints for tirs while storing possible seeds */

/* A pair of TIRs (extended Seeds) */
typedef struct
{
  unsigned long contignumber,
                left_tir_start,  /* first position of TIR on forward strand */ 
                left_tir_end,    /* last position of TIR on forward strand */    
                right_tir_start, /* first position of TIR on reverse strand */  
                right_tir_end;   /* last position of TIR on reverse strand */   
  double        similarity;      /* similarity of the two TIRs */
  bool          skip;           /* needed to remove overlaps if wanted */
  unsigned long tsd_length;     /* length of tsd at start of left tir and end of right tir */
} TIRPair;

GT_DECLAREARRAYSTRUCT(TIRPair);

/* A little struct to store all stuff needed to process TSDs */
typedef struct
{
  unsigned long left_start_pos,   /* represents the start position for the TSD search */
                right_start_pos;
  GtArraySeed TSDs;               /* array to store the TSDs */
}TSDinfo;

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
  /* very important stuff */
  const GtNodeStream          parent_instance;      /* node which could be before us to call next on */  
  const GtEncseq              *encseq;              /* encoded sequence */
  Sequentialsuffixarrayreader *ssar;                /* suffix array reader */
	SeedInfo                 		seedinfo;            	/* struct containing infos for storing seeds and a seedarray */
  GtArrayTIRPair              tir_pairs;            /* array contianing our TIRs */
  GtTIRStreamState            state;                /* current state of output */  
  
  unsigned long               num_of_tirs,          /* number of TIRs found */
                              cur_elem_index,       /* index of the TIR to be put out next */
                              prev_seqnum;          /* number of previous contig */
  
  /* options */
  GtStr                       *str_indexname;       /* index name of suffix array */  
  unsigned long               min_seed_length,      /* minimal seed length */  
                              min_TIR_length,       /* minimal length of TIR */
                              max_TIR_length,       /* maximal length of TIR */
                              min_TIR_distance,     /* minimal distance of TIRs */
                              max_TIR_distance;     /* maximal distance of TIRs */
  Arbitraryscores             arbit_scores;         /* cost of match, mismatch, insertion, deletion */
  int                         xdrop_belowscore;     /* threshold for xdop algorithm when to discard extensions */
  double                      similarity_threshold; /* decides whether a TIR is accepted or not */
  bool                        no_overlaps;          /* if true all overlaps will be deleted */
  bool                        best_overlaps;        /* if nooverlaps false and this true, of overlapping TIRs 
                                                       only the one with best similary will remain */
  unsigned long               min_TSD_length,       /* minimal length of TSD */
                              max_TSD_length,       /* maxlimal length of TSD */
                              vicinity;             /* vicinity to be searched for TSDs */
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
  bool tir_distance = true;									/* if the distance constraints for the tirs can be satisfied */
  SeedInfo *seeds = (SeedInfo *) info;

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
  
	if(direction && samecontig_and_diffstrands)
	{
		if(pos1 - startpos_contig1 > seeds->max_tir_length - len)
		{
			if(((pos1 + offset) - (pos1 - seeds->max_tir_length + len)) < seeds->min_tir_distance)
			{
				tir_distance = false;
			}
		}
		else if(pos1 - startpos_contig1 <= seeds->max_tir_length - len)
		{
			if(pos1 + offset - startpos_contig1 < seeds->min_tir_distance)
			{
				tir_distance = false;
			}
		}
		else if((pos1 + offset - seeds->max_tir_length + len) > pos1)
		{
			if((pos1 + offset - seeds->max_tir_length + len) - pos1 > seeds->max_tir_distance)
			{
				tir_distance = false;
			}
		}
		
	}
  
  /* like this, cause we want to have other constraints later */
  if (samecontig_and_diffstrands && direction && tir_distance)
  {
    Seed *nextfreeseedpointer;
    
    GT_GETNEXTFREEINARRAY(nextfreeseedpointer, &seeds->seed, Seed, 10);
    
    nextfreeseedpointer->pos1 = pos1;
    nextfreeseedpointer->pos2 = pos2;
    nextfreeseedpointer->offset = offset;
    nextfreeseedpointer->len = len;
    nextfreeseedpointer->contignumber = contignumber;
    
  }
  return 0;
}

/* this function is a call back function to store all TSDs found */

static int gt_store_TSDs(void *info, GT_UNUSED const GtEncseq *encseq, const Querymatch *querymatch, GT_UNUSED GtError *err)
{
  Seed *nextfree;
  TSDinfo *TSDs = (TSDinfo *) info;
  
  /* store the TSD at the next free index of info */
  GT_GETNEXTFREEINARRAY(nextfree, &TSDs->TSDs, Seed, 10);
  nextfree->pos1 = TSDs->left_start_pos + gt_querymatch_dbstart(querymatch);
  nextfree->offset = TSDs->right_start_pos + gt_querymatch_querystart(querymatch) - (nextfree->pos1);
  nextfree->len = gt_querymatch_len(querymatch);
  
  return 0;
}

/* compares two tirs */
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

/* removes overlaps or saves tir with best similarity and returns the lasting number of tirs in the array*/

static unsigned long gt_remove_overlaps(GtArrayTIRPair *src, GtArrayTIRPair *dest, unsigned long size_of_array, bool remove_all)
{
  TIRPair *new_pair;
  TIRPair *pair1;
  TIRPair *pair2;
  unsigned long start_pair1;
  unsigned long end_pair1;
  unsigned long start_pair2;
  unsigned long end_pair2;
  unsigned long num_of_tirs = size_of_array;
  int i;
  int j;

  
  for(i = 0; i < size_of_array; i++)
  {
    pair1 = &src->spaceTIRPair[i];
    
    /* to prevent that a skipped one is checked twice */
    if(pair1->skip)
    {
      continue;
    }
    start_pair1 = pair1->left_tir_start;
    end_pair1 = pair1->right_tir_end;
    
    for(j = i + 1; j < size_of_array; j++)
    {
      pair2 = &src->spaceTIRPair[j];
      
      /* to prevent that a skipped one is checked twice */
      if(pair2->skip)
      {
        continue;
      }
      start_pair2 = pair2->left_tir_start;
      end_pair2 = pair2->right_tir_end;
      
      /* check if there is an overlap */
      if(!((end_pair1 < start_pair2) || (end_pair2 < start_pair1)))
      {
        /* no overlaps allowed */
        if(remove_all)
        {
          //TODO verstehen, warum LTRharvest die min_startpos und die max_endpos nimmt
          /* set skip to delete the tirs later */
          if(!pair1->skip)
          {
            num_of_tirs = num_of_tirs - 2;
          }
          else
          {
            num_of_tirs = num_of_tirs - 1;
          }
          pair1->skip = true;
          pair2->skip = true;
        }
        else
        {
          /* take the tir with best similarity */
          if(pair1->similarity >= pair2->similarity)
          {
            pair2->skip = true;
            num_of_tirs = num_of_tirs - 1;
          }
          else
          {
            pair1->skip = true;
            num_of_tirs = num_of_tirs - 1;
            break;
          }
        }
      }
    }
  }
  
  /* the tirs without overlaps store in the new array */
  for(i = 0; i < size_of_array; i++)
  {
    pair1 = &src->spaceTIRPair[i];
    
    if(!pair1->skip)
    {
      GT_GETNEXTFREEINARRAY(new_pair, dest, TIRPair, 5);
      *new_pair = *pair1;
    }
  }
  return num_of_tirs;
}

/* this function finds the best TSD */
static void gt_find_best_TSD(TSDinfo *info,
                             GtTIRStream *tir_stream,
                             TIRPair *tir_pair)
{
   int i;
   int j;
   Seed *tsd;
   unsigned long tsd_length;
   unsigned long optimal_tsd_length;
   unsigned long new_left_tir_start = tir_pair->left_tir_start;
   unsigned long new_right_tir_end = tir_pair->right_tir_end;
   unsigned long new_cost_left = 0;
   unsigned long new_cost_right = 0;
   unsigned long best_cost = ULONG_MAX;
   unsigned long new_cost = 0;
   
   for(i = 0; i < info->TSDs.nextfreeSeed; i++)
   {
      
      tsd = &info->TSDs.spaceSeed[i];
      optimal_tsd_length = tsd->len;
      
	  /* needed, because of overflow */
      if(tsd->len <= tir_stream->min_TSD_length)
      {
		continue;
	  }
      for(j = 0; j < tsd->len - tir_stream->min_TSD_length + 1; j++)
      {
        tsd_length = tsd->len - j;
        /* max len constraint */
        if(tsd_length < tir_stream->max_TSD_length)
        {
          if(tsd->pos1 + tsd_length - 1 < tir_pair->left_tir_start)
          {
            new_cost_left = tir_pair->left_tir_start - (tsd->pos1 + tsd_length - 1);
          }
          else
          {
            new_cost_left = (tsd->pos1 + tsd_length - 1) - tir_pair->left_tir_start;
          }
          
          if(tir_pair->right_tir_end < tsd->pos1 + tsd->offset)
          {
            new_cost_right = (tsd->pos1 + tsd->offset) - tir_pair->right_tir_end;
          }
          else
          {
            new_cost_right = tir_pair->right_tir_end - (tsd->pos1 + tsd->offset);
          }
          
          new_cost = new_cost_left + new_cost_right;
          
          if(new_cost < best_cost)
          {
            best_cost = new_cost;
            new_left_tir_start = tsd->pos1 + tsd_length;
            new_right_tir_end = tsd->pos1 + tsd->offset - 1;
            optimal_tsd_length = tsd_length;
          }
           
        }
      }
   }
   /* save the new borders and tsd length */
   if(info->TSDs.nextfreeSeed > 0)
   {
     tir_pair->left_tir_start = new_left_tir_start;
     tir_pair->right_tir_end = new_right_tir_end;
     tir_pair->tsd_length = optimal_tsd_length;
   }
}


/* this function searches for TSDs in the range of vicinity around the TIRs */
static int gt_search_for_TSDs(GtTIRStream *tir_stream, TIRPair *tir_pair, const GtEncseq *encseq, GtError *err)
{
  unsigned long start_left_tir,
                end_left_tir,
                start_right_tir,
                end_right_tir,
                left_length,
                right_length,
                seq_end_pos,
                seq_start_pos,
                seq_length;
  unsigned long contignumber = tir_pair->contignumber;
  TSDinfo info;
  bool haserr = false;

  gt_error_check(err);

  /* check border cases */

  /* check vicinity for left tir start */
  seq_start_pos = gt_encseq_seqstartpos(encseq, contignumber);
  seq_length = gt_encseq_seqlength(encseq, contignumber);
  
  gt_assert(tir_pair->left_tir_start >= seq_start_pos);
  gt_assert(tir_pair->left_tir_start <= tir_pair->left_tir_end);
  gt_assert(tir_pair->right_tir_start <= tir_pair->right_tir_end);
  gt_assert(tir_pair->right_tir_end <= seq_start_pos + seq_length);
  
  /* check if left tir start with vicinity aligns over sequence border */
  if(tir_pair->left_tir_start - seq_start_pos < tir_stream->vicinity)
  {
    start_left_tir = seq_start_pos;
  }
  else
  {
    start_left_tir = tir_pair->left_tir_start - tir_stream->vicinity;
  }
  
  /* do not align over end of left tir */
  if(tir_pair->left_tir_start + tir_stream->vicinity > tir_pair->left_tir_end)
  {
    end_left_tir = tir_pair->left_tir_end;
  }
  else
  {
    end_left_tir = tir_pair->left_tir_start + tir_stream->vicinity;
  }
  
  left_length = end_left_tir - start_left_tir + 1;

  /* vicinity of 3'-border of right tir
     do not align over 5'border of right tir */
  if(tir_pair->right_tir_end < tir_pair->right_tir_start + tir_stream->vicinity)
  {
    start_right_tir = tir_pair->right_tir_start;
  }
  else
  {
    start_right_tir = tir_pair->right_tir_end - tir_stream->vicinity;
  }
  
  seq_end_pos = seq_start_pos + seq_length - 1;
  
  /* do not align into next sequence in case of need decrease alignment
     length */
  if (tir_pair->right_tir_end + tir_stream->vicinity > seq_end_pos)
  {
    end_right_tir = seq_end_pos;
  }
  else
  {
    end_right_tir = tir_pair->right_tir_end + tir_stream->vicinity;
  }
  
  right_length = end_right_tir - start_right_tir + 1;

  /* search for TSDs */
  if (tir_stream->min_TSD_length > 1U)
  {
    /* dbseq (left) and query (right) are the encseqs wich will be aligned */
    GtUchar *dbseq, *query;
    ALLOCASSIGNSPACE(dbseq,NULL,GtUchar,left_length);
    ALLOCASSIGNSPACE(query,NULL,GtUchar,right_length);

		//TODO vor einem aufruf hier passiert ab einer gewissen größe ein speicherzugriffsfehler -.-
    gt_encseq_extract_encoded(encseq,dbseq,start_left_tir,end_left_tir);
    gt_encseq_extract_encoded(encseq,query,start_right_tir,end_right_tir);

    GT_INITARRAY(&info.TSDs, Seed);
    gt_assert(start_left_tir < start_right_tir);
    info.left_start_pos = start_left_tir;
    info.right_start_pos = start_right_tir; 
    
    
    if (gt_sarrquerysubstringmatch(dbseq,
                                   left_length,
                                   query,
                                   right_length,
                                   tir_stream->min_TSD_length,
                                   gt_encseq_alphabet(encseq),
                                   gt_store_TSDs,
                                   &info,
                                   NULL,
                                   err) != 0)
    {
       haserr = true;
    }

    FREESPACE(dbseq);
    FREESPACE(query);

    /* find the best TSD */
    if (!haserr)
    {
      gt_find_best_TSD(&info, tir_stream, tir_pair);
    }
    GT_FREEARRAY(&info.TSDs, Seed);
  }
  
  return haserr ? -1 : 0;
}


/*
 * Extends the Seeds and searches for TIRs.
 */
static int gt_searchforTIRs(GtTIRStream *tir_stream,
                             GtArrayTIRPair *tir_pairs,
                             const GtEncseq *encseq,
                             GtError *err)
{
  int count;
  unsigned long seedcounter = 0;
  GtArrayTIRPair new;             /* need to remove overlaps */
  unsigned long right_tir_start;  /* need to calculate the reverse position */
  unsigned long right_tir_end;    /* need to calculate the reverse position */
  GtArrayMyfrontvalue fronts; /* needed to use xdrop */
  Myxdropbest xdropbest_left; /* return parameters for xdrop */
  Myxdropbest xdropbest_right;
  unsigned long total_length = 0,
                left_tir_length = 0, 
                right_tir_length = 0,
                max_left_length = 0,    /* maximal length of left TIR*/ 
                max_right_length = 0,
                distance = 0;
  GtUchar *left_tir_char = NULL,        /* next character to align */
          *right_tir_char = NULL;
  unsigned long edist = 0;              /* edit distance */
  Seed *seedptr;
  TIRPair *pair;
  int had_err = 0;
  
  /* Did we have errors so far? */
  gt_error_check(err);
  
  /* Iterating over seeds */
  for(seedcounter = 0; seedcounter < tir_stream->seedinfo.seed.nextfreeSeed;
      seedcounter++)
  {
    /* Getting the seed entry at position seedcounter */
    seedptr = &(tir_stream->seedinfo.seed.spaceSeed[seedcounter]);
    
    /* Initialize array for extending seeds */
    GT_INITARRAY (&fronts, Myfrontvalue);
    
    /* Left alignment */
    gt_evalxdroparbitscoresleft(&tir_stream->arbit_scores,
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
    gt_evalxdroparbitscoresright(&tir_stream->arbit_scores,
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
    right_tir_start = seedptr->pos2 + seedptr->len - 1 + xdropbest_right.jvalue; // end of mirrored is start of actual TIR
    right_tir_end   = seedptr->pos2 - xdropbest_left.jvalue;  // start of mirrored is end of actual TIR
    
    /* We can get the corresponding position of a mirrored one with GT_REVERSEPOS(total length,position) */
    pair->right_tir_start = GT_REVERSEPOS(total_length,right_tir_start);
    pair->right_tir_end = GT_REVERSEPOS(total_length,right_tir_end); 
    pair->similarity = 0.0;
    pair->skip = false;
    
    had_err = gt_search_for_TSDs(tir_stream, pair, encseq, err);
    
       	
   	//TODO hier möchte überlegt werden, wo dieser fehler abgefangen werden möchte
   	if(pair->left_tir_start > pair->left_tir_end || pair->right_tir_start > pair->right_tir_end)
   	{
   		tir_pairs->nextfreeTIRPair--;
   		continue;
   	}
    
    /* Calculate TIR lengths and distance*/
    left_tir_length = pair->left_tir_end - pair->left_tir_start + 1;
    right_tir_length = pair->right_tir_end - pair->right_tir_start + 1;
    
    /* needed because of border change by search_forTSDs */
    right_tir_start = GT_REVERSEPOS(total_length, pair->right_tir_end);
    right_tir_end = GT_REVERSEPOS(total_length, pair->right_tir_start);
    
    gt_assert( pair->right_tir_end - pair->right_tir_start + 1 == (right_tir_end - right_tir_start + 1));
    distance = pair->right_tir_start - pair->left_tir_start;
    
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
                                         right_tir_start,
                                         right_tir_end);
    
    /* Get edit distance */
    edist = greedyunitedist(left_tir_char,(unsigned long) left_tir_length, 
                            right_tir_char,(unsigned long) right_tir_length);
                            
    /* Determine similarity */
    pair->similarity = 100.0 * (1 - (((double) edist)/
                      (MAX(left_tir_length, right_tir_length))));
                      
    /* Discard this candidate if lengths constrains not satisfied */
    if(left_tir_length > tir_stream->max_TIR_length ||
       right_tir_length > tir_stream->max_TIR_length ||
       left_tir_length < tir_stream->min_TIR_length ||
       right_tir_length < tir_stream->min_TIR_length)
    {
      tir_pairs->nextfreeTIRPair--;
    }    
    /* or if distance constraints not satisfies */
    else if(distance > tir_stream->max_TIR_distance ||
            distance < tir_stream->min_TIR_distance)
    {
      tir_pairs->nextfreeTIRPair--;
    }                   
    /* or if similarity too small or increase number of TIRs*/
    else if(gt_double_smaller_double(pair->similarity, 
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
  
  /* initialize array for removing overlaps */
  GT_INITARRAY(&new, TIRPair);
  
  /* remove overlaps if wanted */
  if(tir_stream->best_overlaps || tir_stream->no_overlaps)
  {
    tir_stream->num_of_tirs = gt_remove_overlaps(tir_pairs, &new, tir_stream->num_of_tirs, tir_stream->no_overlaps);
      
    /* set references to the new array and free old array, that is not needed any longer */
    GT_FREEARRAY(tir_pairs, TIRPair);
    tir_stream->tir_pairs = new;
    tir_pairs = &new; 
  }
  
  /* sort the tir_pairs */
  gt_sort_TIRs(tir_pairs, 0, tir_stream->num_of_tirs);
  
    
  /* output on console */
  gt_log_log("TIRs found:\n");
  for(count = 0; count < tir_stream->tir_pairs.nextfreeTIRPair; count++)
  {
    gt_log_log("contig %lu\n left: start %lu, end %lu\n right: start %lu, end %lu\n similarity: %f\n\n", 
      tir_stream->tir_pairs.spaceTIRPair[count].contignumber,
      tir_stream->tir_pairs.spaceTIRPair[count].left_tir_start,
      tir_stream->tir_pairs.spaceTIRPair[count].left_tir_end,
      tir_stream->tir_pairs.spaceTIRPair[count].right_tir_start,
      tir_stream->tir_pairs.spaceTIRPair[count].right_tir_end,
      tir_stream->tir_pairs.spaceTIRPair[count].similarity);
  }
  
  return had_err;
  
}


/*
 * Saves the next node of the annotation graph in gn.
 * Search for Seeds here.
 */
static int gt_tir_stream_next(GtNodeStream *ns,
                              GT_UNUSED GtGenomeNode **gn,
                              GtError *err)
{
  GtTIRStream *tir_stream;
  int had_err = 0;
  gt_error_check(err);
  tir_stream = gt_node_stream_cast(gt_tir_stream_class(), ns);
  
  /* generate and check seeds */
  if(tir_stream->state == GT_TIR_STREAM_STATE_START)
  {
    if (!had_err && gt_enumeratemaxpairs(tir_stream->ssar,
                      tir_stream->encseq,
                      gt_readmodeSequentialsuffixarrayreader(tir_stream->ssar),
                      (unsigned int) tir_stream->min_seed_length,
                      gt_store_seeds,
                      &tir_stream->seedinfo,
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
    GT_FREEARRAY(&tir_stream->seedinfo.seed, Seed);
    

    tir_stream->state = GT_TIR_STREAM_STATE_REGIONS;
  }

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
      node = gt_feature_node_new(seqid, "repeat_region", pair->left_tir_start - seqstartpos + 1, pair->right_tir_end - seqstartpos + 1, GT_STRAND_UNKNOWN);
      
      /* set stuff */
      gt_feature_node_set_source((GtFeatureNode*) node, source);
      *gn = node;
      parent = node;

      /* target site duplication */
      
      // TODO child of region
      
      
      /* terminal inverted repeat element */
      
      node = gt_feature_node_new(seqid, "terminal_inverted_repeat_element", pair->left_tir_start - seqstartpos + 1, pair->right_tir_end - seqstartpos +1, GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*)node, source);
      gt_feature_node_add_child((GtFeatureNode*)parent, (GtFeatureNode*)node);
      parent = node;
      
      /* left terminal inverted repeat */

      node = gt_feature_node_new(seqid, "terminal_inverted_repeat", pair->left_tir_start - seqstartpos + 1, pair->left_tir_end - seqstartpos + 1, GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*)node, source);
      gt_feature_node_add_child((GtFeatureNode*)parent, (GtFeatureNode*)node);
      
      
      /* right terminal inverted repeat */
      
      node = gt_feature_node_new(seqid, "terminal_inverted_repeat", pair->right_tir_start - seqstartpos + 1, pair->right_tir_end - seqstartpos + 1, GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*)node, source);
      gt_feature_node_add_child((GtFeatureNode*)parent, (GtFeatureNode*)node);
          
      
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
                                unsigned long min_seed_length,
                                unsigned long min_TIR_length,
                                unsigned long max_TIR_length,
                                unsigned long min_TIR_distance,
                                unsigned long max_TIR_distance,
                                Arbitraryscores arbit_scores,
                                int xdrop_belowscore,
                                double similarity_threshold,
                                bool best_overlaps,
                                bool no_overlaps,
                                unsigned long min_TSD_length,
                                unsigned long max_TSD_length,
                                unsigned long vicinity,
                                GtError *err)
{
  int had_err = 0;
  GtNodeStream *ns = gt_node_stream_create(gt_tir_stream_class(), false);
  GtTIRStream *tir_stream = gt_node_stream_cast(gt_tir_stream_class(), ns);
  tir_stream->num_of_tirs = 0;
  tir_stream->cur_elem_index = 0;
  tir_stream->prev_seqnum = GT_UNDEF_ULONG;
  tir_stream->state = GT_TIR_STREAM_STATE_START;
  
  tir_stream->str_indexname = gt_str_ref(str_indexname);  /* use of gt_str_ref so the reference counter will be increased */
  tir_stream->min_seed_length = min_seed_length;
  tir_stream->min_TIR_length = min_TIR_length;
  tir_stream->max_TIR_length = max_TIR_length;
  tir_stream->min_TIR_distance = min_TIR_distance;
  tir_stream->max_TIR_distance = max_TIR_distance;
  tir_stream->arbit_scores = arbit_scores;
  tir_stream->xdrop_belowscore = xdrop_belowscore;
  tir_stream->similarity_threshold = similarity_threshold;
  tir_stream->best_overlaps = best_overlaps;
  tir_stream->no_overlaps = no_overlaps;
  tir_stream->min_TSD_length = min_TSD_length;
  tir_stream->max_TSD_length = max_TSD_length;
  tir_stream->vicinity = vicinity;

	tir_stream->seedinfo.max_tir_length = max_TIR_length;
	tir_stream->seedinfo.min_tir_length = min_TIR_length;
	tir_stream->seedinfo.max_tir_distance = max_TIR_distance;
	tir_stream->seedinfo.min_tir_distance = min_TIR_distance;
	
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
  GT_INITARRAY(&tir_stream->seedinfo.seed, Seed);
  GT_INITARRAY(&tir_stream->tir_pairs, TIRPair);
  
  
  if(!had_err)
  {
    return ns;
  }
  
  // TODO ueberlegen, ob had_err benutzt werden müsste
  return NULL;
}
