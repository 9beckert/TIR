/*
  Copyright (c) 2010-2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007      David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
  Copyright (c) 2007-2011 Center for Bioinformatics, University of Hamburg

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
#include "core/array_api.h"
#include "core/arraydef.h"
#include "core/assert_api.h"
#include "core/encseq_api.h"
#include "core/error_api.h"
#include "core/log.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/str_api.h"
#include "core/types_api.h"
#include "core/undef_api.h"
#include "extended/genome_node.h"
#include "extended/node_stream_api.h"
#include "match/esa-seqread.h"
#include "match/esa-maxpairs.h"
#include "match/esa-mmsearch.h"
#include "match/greedyedist.h"
#include "match/spacedef.h"
#include "ltr/ltrharvest_stream.h"
#include "ltr/ltr_xdrop.h"

#define GT_LTRHARVEST_NAME "LTRharvest"

typedef struct
{
  unsigned long pos1;         /* first position of maximal repeat (seed) */
  unsigned long offset;       /* second position = pos1 + offset */
  unsigned long len;          /* length of maximal repeat  */
  unsigned long contignumber; /* number of contig for this repeat */
} Repeat;

GT_DECLAREARRAYSTRUCT(Repeat);

/* The datatype RepeatInfo stores all maximal repeats (seeds) and */
/* information about the length and distance constraints. */
typedef struct
{
  GtArrayRepeat repeats; /* array of maximal repeats (seeds) */
  unsigned long lmin;        /* minimum allowed length of a LTR */
  unsigned long lmax;        /* maximum allowed length of a LTR */
  unsigned long dmin;        /* minimum distance between LTRs */
  unsigned long dmax;        /* maximum distance between LTRs */
  const GtEncseq *encseq;
  GtRange ltrsearchseqrange; /* if start and end are 0, then no range */
} RepeatInfo;

/* The datatype SubRepeatInfo stores information about the maximal repeats */
/* for the TSD detection. */
typedef struct
{
  GtArrayRepeat repeats; /* array of maximal repeats for TSDs */
  unsigned long lmin;   /* minimal length of TSD */
  unsigned long lmax;   /* maximal length of TSD */
  unsigned long offset1;      /* offset1 for absolute position 1 in sequence */
  unsigned long offset2;      /* offset2 for absolute position 2 in sequence */
                       /* pos1 < pos2 */
} SubRepeatInfo;

typedef struct
{
  unsigned long contignumber;   /* ordinal number of sequence in encseq */
  unsigned long leftLTR_5,      /* 5' boundary of left LTR */
                leftLTR_3,      /* 3' boundary of left LTR */
                rightLTR_5,     /* 5' boundary of right LTR */
                rightLTR_3,     /* 3' boundary of right LTR */
                lenleftTSD,
                lenrightTSD;
  bool          tsd,            /* If true, then TSDs exist. */
                motif_near_tsd, /* If true, then motif near the TSD exists. */
                motif_far_tsd,  /* If true, then motif at the inner borders
                                   of LTRs exist. */
                lengthdistconstraint; /* If true, length and distance
                                         constraints are satisfied */
  double similarity;            /* similarity value of LTRs */
  bool skipped;                 /* if skipped then because of an overlap
                                   with a higher similarity prediction or
                                   because of "overlap=no" option */
} LTRboundaries;

GT_DECLAREARRAYSTRUCT(LTRboundaries);

typedef enum {
  GT_LTRHARVEST_STREAM_STATE_START,
  GT_LTRHARVEST_STREAM_STATE_REGIONS,
  GT_LTRHARVEST_STREAM_STATE_COMMENTS,
  GT_LTRHARVEST_STREAM_STATE_FEATURES
} GtLTRharvestStreamState;

struct GtLTRharvestStream
{
  const GtNodeStream parent_instance;
  RepeatInfo repeatinfo;
  GtStr *str_indexname;
  unsigned long minseedlength;
  double similaritythreshold;
  int xdropbelowscore;
  Arbitraryscores arbitscores;
  unsigned int allowedmismatches;
  GtLTRFourCharMotif *motif;
  unsigned long offset;
  unsigned int minlengthTSD,
               maxlengthTSD;
  unsigned long numofboundaries;
  unsigned long vicinityforcorrectboundaries;
  unsigned long prevseqnum;
  const LTRboundaries **bdptrtab;
  GtArrayLTRboundaries arrayLTRboundaries;
  const GtEncseq *encseq;
  Sequentialsuffixarrayreader *ssar;
  bool verbosemode;
  bool nooverlaps;
  bool bestoverlaps;
  unsigned long cur_elem_index;
  GtLTRharvestStreamState state;
};

#define gt_ltrharvest_stream_cast(GS)\
        gt_node_stream_cast(gt_ltrharvest_stream_class(), GS);

static int bdptrcompare(const void *a, const void *b)
{
  const LTRboundaries **bda, **bdb;

  bda = (const LTRboundaries **) a;
  bdb = (const LTRboundaries **) b;
  if ((*bda)->contignumber < (*bdb)->contignumber)
  {
    return -1;
  }
  if ((*bda)->contignumber > (*bdb)->contignumber)
  {
    return 1;
  }
  if ((*bda)->leftLTR_5 < (*bdb)->leftLTR_5)
  {
    return -1;
  }
  if ((*bda)->leftLTR_5 > (*bdb)->leftLTR_5)
  {
    return 1;
  }
  return 0;
}

static int gt_simpleexactselfmatchstore(void *info,
                                         GT_UNUSED const GtEncseq *encseq,
                                         unsigned long len,
                                         unsigned long pos1,
                                         unsigned long pos2,
                                         GT_UNUSED GtError *err)
{
  unsigned long tmp;
  unsigned long contignumber = 0,
                seqnum1,
                seqnum2;
  bool samecontig = false;
  RepeatInfo *repeatinfo = (RepeatInfo *) info;

  gt_error_check(err);
  if (pos1 > pos2)
  {
    tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
  if (repeatinfo->ltrsearchseqrange.start != 0 ||
      repeatinfo->ltrsearchseqrange.end != 0)
  {
    if (pos1 < (unsigned long) repeatinfo->ltrsearchseqrange.start  ||
        pos2 + len - 1 > (unsigned long) repeatinfo->ltrsearchseqrange.end)
    {
      return 0;
    }
  }

  tmp = pos2 - pos1;
  seqnum1 = gt_encseq_seqnum(repeatinfo->encseq,pos1);
  seqnum2 = gt_encseq_seqnum(repeatinfo->encseq,pos2);
  if (seqnum1 == seqnum2)
  {
    gt_log_log("accepted:\n");
    gt_log_log("pos1: %lu\n", pos1);
    gt_log_log("pos2: %lu\n", pos2);
    samecontig = true;
    contignumber = seqnum1;
  }

  /*test maximal length of candidate pair and distance constraints*/
  if (samecontig
        && len <= (unsigned long) repeatinfo->lmax
        && (unsigned long) repeatinfo->dmin <= tmp
        && tmp <= (unsigned long) repeatinfo->dmax)
  {
    Repeat *nextfreerepeatptr;

    GT_GETNEXTFREEINARRAY(nextfreerepeatptr, &repeatinfo->repeats,
                       Repeat, 10);
    gt_log_log("maximal repeat pos1: %lu\n",
               pos1);
    gt_log_log("maximal repeat pos2: %lu\n",
               pos2);
    gt_log_log("len: %lu\n", len);
    gt_log_log("seq number: %lu\n\n", contignumber);
    nextfreerepeatptr->pos1 = pos1;
    nextfreerepeatptr->offset = tmp;
    nextfreerepeatptr->len = len;
    nextfreerepeatptr->contignumber = contignumber;
  }
  return 0;
}

static int gt_subsimpleexactselfmatchstore(void *info,
                                            GT_UNUSED const GtEncseq *encseq,
                                            const Querymatch *querymatch,
                                            GT_UNUSED GtError *err)
{
  Repeat *nextfreerepeatptr;
  SubRepeatInfo *sri = (SubRepeatInfo *) info;

  GT_GETNEXTFREEINARRAY (nextfreerepeatptr, &sri->repeats, Repeat, 10);
  nextfreerepeatptr->pos1 = sri->offset1 + gt_querymatch_dbstart(querymatch);
  nextfreerepeatptr->offset = sri->offset2 +
                             gt_querymatch_querystart(querymatch) -
                             (sri->offset1 + gt_querymatch_dbstart(querymatch));
  nextfreerepeatptr->len = gt_querymatch_len(querymatch);
  return 0;
}

/* XXX: better directly sort the ArrayLTRboundaries */
static const LTRboundaries **sortedltrboundaries(unsigned long *numofboundaries,
                                                 const GtArrayLTRboundaries
                                                 *ltr)
{
  unsigned long countboundaries = 0, nextfill = 0;
  const LTRboundaries *bd, **bdptrtab;

  for (bd = ltr->spaceLTRboundaries; bd < ltr->spaceLTRboundaries +
                                          ltr->nextfreeLTRboundaries; bd++)
  {
    if (!bd->skipped)
    {
      countboundaries++;
    }
  }
  bdptrtab = gt_malloc(sizeof (LTRboundaries *) * countboundaries);
  nextfill = 0;
  for (bd = ltr->spaceLTRboundaries; bd < ltr->spaceLTRboundaries +
                                          ltr->nextfreeLTRboundaries; bd++)
  {
    if (!bd->skipped)
    {
      bdptrtab[nextfill++] = bd;
    }
  }
  qsort(bdptrtab,(size_t) countboundaries, sizeof (LTRboundaries *),
        bdptrcompare);
  *numofboundaries = countboundaries;
  return bdptrtab;
}

/*
 The following function searches for TSDs and/or a specified palindromic
 motif at the 5'-border of left LTR and 3'-border of right LTR. Thereby,
 all maximal repeats from the vicinity are processed one after another
 to find the TSD with the minimum deviation with regard to the boundary
 position from the x-drop alignment. If also a motif is searched,
 a simple motif check at the boundaries of the TSDs is performed.
 */

static void searchforbestTSDandormotifatborders(SubRepeatInfo *info,
                                                GtLTRharvestStream *lo,
                                                const GtEncseq *encseq,
                                                LTRboundaries *boundaries,
                                                unsigned int
                                                *motifmismatchesleftLTR,
                                                unsigned int
                                                *motifmismatchesrightLTR)
{
  unsigned long i;
  unsigned long motifpos1,
         motifpos2,
         back, forward,
         oldleftLTR_5  = boundaries->leftLTR_5,
         oldrightLTR_3 = boundaries->rightLTR_3,
         difffromoldboundary1 = 0,
         difffromoldboundary2 = 0;
  unsigned int tmp_motifmismatchesleftLTR,
               tmp_motifmismatchesrightLTR,
               hitcounter = 0;

  if (info->repeats.nextfreeRepeat > 0 )
  {
    boundaries->tsd = true;
  }
  boundaries->motif_near_tsd = false;

  for (i = 0; i < info->repeats.nextfreeRepeat; i++)
  {
    /* motifpos1 is the first position after the left repeat */
    motifpos1 = info->repeats.spaceRepeat[i].pos1 +
                           info->repeats.spaceRepeat[i].len;
    /* motifpos2 is two positions before the right repeat */
    motifpos2 = info->repeats.spaceRepeat[i].pos1
                  + info->repeats.spaceRepeat[i].offset - 2;

    for (back = 0;
         back < info->repeats.spaceRepeat[i].len - info->lmin + 1;
         back++)
    {
      for (forward = 0;
           forward < info->repeats.spaceRepeat[i].len -
                     info->lmin + 1 - back;
           forward++)
      {
        tmp_motifmismatchesleftLTR = tmp_motifmismatchesrightLTR = 0;
        if (gt_encseq_get_encoded_char(/* Random access */ encseq,
                            motifpos1 - back, GT_READMODE_FORWARD)
            != lo->motif->firstleft)
        {
          tmp_motifmismatchesleftLTR++;
        }
        if (gt_encseq_get_encoded_char(/* Random access */ encseq,
                                              motifpos1 + 1 - back,
                                              GT_READMODE_FORWARD)
            != lo->motif->secondleft)
        {
          tmp_motifmismatchesleftLTR++;
        }
        if (gt_encseq_get_encoded_char(/* Random access */ encseq,
                                              motifpos2 + forward,
                                              GT_READMODE_FORWARD)
            != lo->motif->firstright)
        {
          tmp_motifmismatchesrightLTR++;
        }
        if (gt_encseq_get_encoded_char(/* Random access */ encseq,
                                              motifpos2 + 1 + forward,
                                              GT_READMODE_FORWARD)
            != lo->motif->secondright)
        {
          tmp_motifmismatchesrightLTR++;
        }

        if (tmp_motifmismatchesleftLTR <= lo->motif->allowedmismatches
            &&
            tmp_motifmismatchesrightLTR <= lo->motif->allowedmismatches)
        {
          unsigned long tsd_len;
          tsd_len = info->repeats.spaceRepeat[i].len - back - forward;

          /* TSD length not too big */
          if (tsd_len <= (unsigned long)info->lmax)
          {
            if ( !boundaries->motif_near_tsd )
            {
              unsigned long max, min;

              /* save number of mismatches */
              *motifmismatchesleftLTR  = tmp_motifmismatchesleftLTR;
              *motifmismatchesrightLTR = tmp_motifmismatchesrightLTR;

              /* adjust boundaries */
              boundaries->motif_near_tsd = true;
              boundaries->leftLTR_5  = motifpos1 - back;
              boundaries->rightLTR_3 = motifpos2 + 1 + forward;

              /* store TSD length */
              boundaries->lenleftTSD = boundaries->lenrightTSD = tsd_len;

              max = MAX(oldleftLTR_5, boundaries->leftLTR_5);
              min = MIN(oldleftLTR_5, boundaries->leftLTR_5);
              difffromoldboundary1 = max - min;

              max = MAX(oldrightLTR_3, boundaries->rightLTR_3);
              min = MIN(oldrightLTR_3, boundaries->rightLTR_3);
              difffromoldboundary2 = max - min;

              hitcounter++;
            }
            else
            {
              unsigned long max, min, difffromnewboundary1,
                   difffromnewboundary2;

              /* test if hit is nearer to old boundaries than previous hit */
              max = MAX(oldleftLTR_5, (motifpos1 - back));
              min = MIN(oldleftLTR_5, (motifpos1 - back));
              difffromnewboundary1 = max - min;
              max = MAX(oldrightLTR_3, (motifpos2 + 1 + forward));
              min = MIN(oldrightLTR_3, (motifpos2 + 1 + forward));
              difffromnewboundary2 = max - min;

              if (difffromnewboundary1 + difffromnewboundary2 <
                  difffromoldboundary1 + difffromoldboundary2)
              {
                /* save number of mismatches */
                *motifmismatchesleftLTR  = tmp_motifmismatchesleftLTR;
                *motifmismatchesrightLTR = tmp_motifmismatchesrightLTR;

                /* adjust boundaries */
                boundaries->leftLTR_5  = motifpos1 - back;
                boundaries->rightLTR_3 = motifpos2 + 1 + forward;

                /* store TSD length */
                boundaries->lenleftTSD = boundaries->lenrightTSD = tsd_len;

                difffromoldboundary1 = difffromnewboundary1;
                difffromoldboundary2 = difffromnewboundary2;
                hitcounter++;
              }
            }
          }
        }
      }
    }
  }
}

/*
 The following function searches only for a specified palindromic motif
 at the 5'-border of left LTR and 3'-border of right LTR.
 */

static void searchformotifonlyborders(GtLTRharvestStream *lo,
    LTRboundaries *boundaries,
    const GtEncseq *encseq,
    unsigned long startleftLTR,
    unsigned long endleftLTR,
    unsigned long startrightLTR,
    unsigned long endrightLTR,
    unsigned int *motifmismatchesleftLTR,
    unsigned int *motifmismatchesrightLTR
    )
{
  bool motif1 = false,
       motif2 = false;
  unsigned int tmp_motifmismatchesleftLTR,
         tmp_motifmismatchesrightLTR,
         motifmismatches_frombestmatch = 0;
  unsigned long idx,
         oldleftLTR_5  = boundaries->leftLTR_5,
         oldrightLTR_3 = boundaries->rightLTR_3,
         difffromoldboundary = 0;

  /**** search for left motif around leftLTR_5 ****/

  for (idx = startleftLTR; idx < endleftLTR; idx++)
  {
    tmp_motifmismatchesleftLTR = 0;
    if (gt_encseq_get_encoded_char(/* XXX */ encseq, idx,
                                          GT_READMODE_FORWARD)
        != lo->motif->firstleft)
    {
      tmp_motifmismatchesleftLTR++;
    }
    if (gt_encseq_get_encoded_char(/* XXX */ encseq,
                                          idx+1,
                                          GT_READMODE_FORWARD) !=
        lo->motif->secondleft)
    {
      tmp_motifmismatchesleftLTR++;
    }
    if (tmp_motifmismatchesleftLTR + (*motifmismatchesleftLTR)
                                <= lo->motif->allowedmismatches)
    {
       /* first hit */
       if ( !motif1 )
       {
         unsigned long max, min;

         motifmismatches_frombestmatch = tmp_motifmismatchesleftLTR;
         boundaries->leftLTR_5 = idx;
         motif1 = true;
         max = MAX(oldleftLTR_5, boundaries->leftLTR_5);
         min = MIN(oldleftLTR_5, boundaries->leftLTR_5);
         difffromoldboundary = max - min;
       }
       /* next hit */
       else
       {
         unsigned long maxval, minval, difffromnewboundary;

         /* test if hit is nearer to old boundaries than previous hit */
         maxval = MAX(oldleftLTR_5, idx);
         minval = MIN(oldleftLTR_5, idx);
         difffromnewboundary = maxval - minval;

         if ( difffromnewboundary < difffromoldboundary )
         {
           motifmismatches_frombestmatch = tmp_motifmismatchesleftLTR;
           boundaries->leftLTR_5 = idx;
           difffromoldboundary = difffromnewboundary;
         }
       }
    }
  }
  *motifmismatchesleftLTR += motifmismatches_frombestmatch;
  motifmismatches_frombestmatch = 0;

  for (idx = startrightLTR + 1; idx <= endrightLTR; idx++)
  {
    tmp_motifmismatchesrightLTR = 0;
    if (gt_encseq_get_encoded_char(/* XXX */ encseq, idx,
                                          GT_READMODE_FORWARD) !=
                       lo->motif->secondright)
    {
      tmp_motifmismatchesrightLTR++;
    }
    if (gt_encseq_get_encoded_char(/* XXX */ encseq,
                                          idx-1,
                                          GT_READMODE_FORWARD) !=
                       lo->motif->firstright)
    {
      tmp_motifmismatchesrightLTR++;
    }

    if (tmp_motifmismatchesrightLTR + (*motifmismatchesrightLTR)
                                   <= lo->motif->allowedmismatches)
    {
       /* first hit */
       if ( !motif2 )
       {
         unsigned long max, min;

         motifmismatches_frombestmatch = tmp_motifmismatchesrightLTR;
         boundaries->rightLTR_3 = idx;
         motif2 = true;
         max = MAX(oldrightLTR_3, boundaries->rightLTR_3);
         min = MIN(oldrightLTR_3, boundaries->rightLTR_3);
         difffromoldboundary = max - min;
       }
       /* next hit */
       else
       {
         unsigned long maxval, minval, difffromnewboundary;

         /* test if hit is nearer to old boundaries than previous hit */
         maxval = MAX(oldrightLTR_3, idx);
         minval = MIN(oldrightLTR_3, idx);
         difffromnewboundary = maxval - minval;
         if ( difffromnewboundary < difffromoldboundary )
         {
           motifmismatches_frombestmatch = tmp_motifmismatchesrightLTR;
           boundaries->rightLTR_3 = idx;
           difffromoldboundary = difffromnewboundary;
         }
       }
    }
  }
  if (idx > endrightLTR && (!motif2))
  {
    gt_log_log("no right motif found.\n");
  }
  *motifmismatchesrightLTR += motifmismatches_frombestmatch;

  if (motif1 && motif2)
  {
    boundaries->motif_near_tsd = true;
  }
  else
  {
    boundaries->motif_near_tsd = false;
  }
}

/*
 The following function searches for a specified palindromic motif at the
 3'-border of left LTR and the 5'-border of right LTR.
 */

static void searchformotifonlyinside(GtLTRharvestStream *lo,
                                     LTRboundaries *boundaries,
                                     const GtEncseq *encseq,
                                     unsigned int *motifmismatchesleftLTR,
                                     unsigned int *motifmismatchesrightLTR)
{
  bool motif1 = false,
       motif2 = false;
  unsigned long startleftLTR,
         endleftLTR,
         startrightLTR,
         endrightLTR,
         oldleftLTR_3  = boundaries->leftLTR_3,
         oldrightLTR_5 = boundaries->rightLTR_5,
         difffromoldboundary = 0,
         idx;
  unsigned int tmp_motifmismatchesleftLTR,
               tmp_motifmismatchesrightLTR,
               motifmismatches_frombestmatch = 0;

  /** vicinity of 3'-border of left LTR **/
  /* do not align over 5'border of left LTR,
     in case of need decrease alignment length */
  if ((boundaries->leftLTR_3 < lo->vicinityforcorrectboundaries)
      || (startleftLTR = boundaries->leftLTR_3 -
         lo->vicinityforcorrectboundaries + 1) <
      boundaries->leftLTR_5 + 2)
  {
    startleftLTR = boundaries->leftLTR_5 + 2;
  }
  /* do not align over 5'-border of right LTR */
  if ((endleftLTR = boundaries->leftLTR_3 +
       lo->vicinityforcorrectboundaries - 1) >
      boundaries->rightLTR_5 - 1)
  {
    endleftLTR = boundaries->rightLTR_5 - 1;
  }
  /** vicinity of 5'-border of right LTR **/
  /* do not align over 3'-border of left LTR */
  if ((startrightLTR = boundaries->rightLTR_5 -
         lo->vicinityforcorrectboundaries + 1)
       < boundaries->leftLTR_3 + 1)
  {
    startrightLTR = boundaries->leftLTR_3 + 1;
  }
  /* do not align over 3'border of right LTR */
  if ((endrightLTR = boundaries->rightLTR_5 +
       lo->vicinityforcorrectboundaries - 1) >
      boundaries->rightLTR_3 - 2)
  {
    endrightLTR = boundaries->rightLTR_3 - 2;
  }

  /**** search for right motif around leftLTR_3 ****/

  for (idx = startleftLTR + 1; idx <= endleftLTR; idx++)
  {
    tmp_motifmismatchesleftLTR = (unsigned int)0;
    if (gt_encseq_get_encoded_char(/* XXX */ encseq, idx,
                                          GT_READMODE_FORWARD)
                       != lo->motif->secondright)
    {
      tmp_motifmismatchesleftLTR++;
    }
    if (gt_encseq_get_encoded_char(/* XXX */ encseq,
                                          idx-1,
                                          GT_READMODE_FORWARD) !=
                       lo->motif->firstright)
    {
      tmp_motifmismatchesleftLTR++;
    }
    if (tmp_motifmismatchesleftLTR + (*motifmismatchesleftLTR)
                                <= lo->motif->allowedmismatches)
    {
       /* first hit */
       if ( !motif1 )
       {
         unsigned long maxval, minval;

         motifmismatches_frombestmatch = tmp_motifmismatchesleftLTR;
         boundaries->leftLTR_3 = idx;
         motif1 = true;
         maxval = MAX(oldleftLTR_3, boundaries->leftLTR_3);
         minval = MIN(oldleftLTR_3, boundaries->leftLTR_3);
         difffromoldboundary = maxval - minval;
       }
       /* next hit */
       else
       {
         unsigned long maxval, minval, difffromnewboundary;

         /* test if hit is nearer to old boundaries than previous hit */
         maxval = MAX(oldleftLTR_3, idx);
         minval = MIN(oldleftLTR_3, idx);
         difffromnewboundary = maxval - minval;

         if ( difffromnewboundary < difffromoldboundary )
         {
           motifmismatches_frombestmatch = tmp_motifmismatchesleftLTR;
           boundaries->leftLTR_3 = idx;
           difffromoldboundary = difffromnewboundary;
         }
       }
    }
  }
  *motifmismatchesleftLTR += motifmismatches_frombestmatch;
  motifmismatches_frombestmatch = 0;

  /**** search for left motif around rightLTR_5 ****/

  for (idx = startrightLTR ; idx < endrightLTR; idx++)
  {
    tmp_motifmismatchesrightLTR = 0;
    if (gt_encseq_get_encoded_char(/* XXX */ encseq, idx,
                                          GT_READMODE_FORWARD)
                       != lo->motif->firstleft)
    {
      tmp_motifmismatchesrightLTR++;
    }
    if (gt_encseq_get_encoded_char(/* XXX */ encseq, idx+1,
                                          GT_READMODE_FORWARD)
                       != lo->motif->secondleft)
    {
      tmp_motifmismatchesrightLTR++;
    }
    if (tmp_motifmismatchesrightLTR + (*motifmismatchesrightLTR)
                                   <= lo->motif->allowedmismatches)
    {
       /* first hit */
       if ( !motif2 )
       {
         unsigned long maxval, minval;

         motifmismatches_frombestmatch = tmp_motifmismatchesrightLTR;
         boundaries->rightLTR_5 = idx;
         motif2 = true;
         maxval = MAX(oldrightLTR_5, boundaries->rightLTR_5);
         minval = MIN(oldrightLTR_5, boundaries->rightLTR_5);
         difffromoldboundary = maxval - minval;
       }
       /* next hit */
       else
       {
         unsigned long maxval, minval, difffromnewboundary;

         /* test if hit is nearer to old boundaries than previous hit */
         maxval = MAX(oldrightLTR_5, idx);
         minval = MIN(oldrightLTR_5, idx);
         difffromnewboundary = maxval - minval;

         if ( difffromnewboundary < difffromoldboundary )
         {
           motifmismatches_frombestmatch = tmp_motifmismatchesrightLTR;
           boundaries->rightLTR_5 = idx;
           difffromoldboundary = difffromnewboundary;
         }
       }
    }
  }
  *motifmismatchesrightLTR += motifmismatches_frombestmatch;
  if (motif1 && motif2)
  {
    boundaries->motif_far_tsd = true;
  }
  else
  {
    boundaries->motif_far_tsd = false;
  }
}

/*
 The following function searches for TSDs and/or a specified palindromic motif
 at the 5'-border of left LTR and 3'-border of right LTR.
 */

static int searchforTSDandorMotifoutside(
  GtLTRharvestStream *lo,
  LTRboundaries *boundaries,
  const GtEncseq *encseq,
  unsigned int *motifmismatchesleftLTR,
  unsigned int *motifmismatchesrightLTR,
  GtError *err)
{
  unsigned long startleftLTR,
         endleftLTR,
         startrightLTR,
         endrightLTR,
         leftlen,
         rightlen,
         sequenceendpos,
         seqstartpos,
         seqlength;
  unsigned long contignumber = boundaries->contignumber;
  SubRepeatInfo subrepeatinfo;
  bool haserr = false;

  gt_error_check(err);

  /* check border cases */

  /* vicinity of 5'-border of left LTR */
  seqstartpos = gt_encseq_seqstartpos(encseq, contignumber);
  seqlength = gt_encseq_seqlength(encseq, contignumber);
  if (contignumber == 0)
  {
    /* do not align over left sequence boundary,
       in case of need decrease alignment length */
    if ( boundaries->leftLTR_5 < lo->vicinityforcorrectboundaries)
    {
      startleftLTR = seqstartpos;
    } else
    {
      startleftLTR = boundaries->leftLTR_5 - lo->vicinityforcorrectboundaries;
    }
  }
  else
  {
    /* do not align over left separator symbol
       in case of need decrease alignment length */
    if ( boundaries->leftLTR_5 < lo->vicinityforcorrectboundaries )
    {
      startleftLTR = seqstartpos;
    }
    else
    {
      if ( ((startleftLTR =
              boundaries->leftLTR_5 - lo->vicinityforcorrectboundaries) <
                seqstartpos)
            &&
          (boundaries->leftLTR_5 >= seqstartpos)
        )
      {
        startleftLTR = seqstartpos;
      }
    }
  }
  /* do not align over 3'-border of left LTR */
  if ( (endleftLTR = boundaries->leftLTR_5 + lo->vicinityforcorrectboundaries)
       > boundaries->leftLTR_3 - 2 /* -2 because of possible motif */
    )
  {
    endleftLTR = boundaries->leftLTR_3 - 2;
  }
  leftlen = endleftLTR - startleftLTR + 1;

  /* vicinity of 3'-border of right LTR
     do not align over 5'border of right LTR */
  if ( (startrightLTR =
         boundaries->rightLTR_3 - lo->vicinityforcorrectboundaries) <
      boundaries->rightLTR_5 + 2  /* +2 because of possible motif */
    )
  {
    startrightLTR = boundaries->rightLTR_5 + 2;
  }
  sequenceendpos = seqstartpos + seqlength - 1;
  /* do not align into next sequence in case of need decrease alignment
     length */
  endrightLTR = boundaries->rightLTR_3 + lo->vicinityforcorrectboundaries;
  if (endrightLTR > sequenceendpos &&
      boundaries->rightLTR_3 <= sequenceendpos)
  {
    endrightLTR = sequenceendpos;
  }
  rightlen = endrightLTR - startrightLTR + 1;

  /* now, search for correct boundaries */

  /* search for TSDs and/or motif */
  if (lo->minlengthTSD > 1U)
  {
    GtUchar *dbseq, *query;
    ALLOCASSIGNSPACE(dbseq,NULL,GtUchar,leftlen);
    ALLOCASSIGNSPACE(query,NULL,GtUchar,rightlen);

    gt_encseq_extract_encoded(encseq,dbseq,startleftLTR,endleftLTR);
    gt_encseq_extract_encoded(encseq,query,startrightLTR,
                                         endrightLTR);
    GT_INITARRAY(&subrepeatinfo.repeats, Repeat);
    subrepeatinfo.lmin = (unsigned long) lo->minlengthTSD;
    subrepeatinfo.lmax = (unsigned long) lo->maxlengthTSD;
    gt_assert(startleftLTR < startrightLTR);
    subrepeatinfo.offset1 = startleftLTR;
    subrepeatinfo.offset2 = startrightLTR;

    if (gt_sarrquerysubstringmatch(dbseq,
                                leftlen,
                                query,
                                (unsigned long) rightlen,
                                lo->minlengthTSD,
                                gt_encseq_alphabet(encseq),
                                gt_subsimpleexactselfmatchstore,
                                &subrepeatinfo,
                                NULL,
                                err) != 0)
    {
       haserr = true;
    }

    FREESPACE(dbseq);
    FREESPACE(query);

    if (!haserr)
    {
      searchforbestTSDandormotifatborders(&subrepeatinfo,
                                          lo,
                                          encseq,
                                          boundaries,
                                          motifmismatchesleftLTR,
                                          motifmismatchesrightLTR);
    }
    GT_FREEARRAY (&subrepeatinfo.repeats, Repeat);
  } else /* no search for TSDs, search for motif only */
  {
    searchformotifonlyborders(lo,
                              boundaries,
                              encseq,
                              startleftLTR,
                              endleftLTR,
                              startrightLTR,
                              endrightLTR,
                              motifmismatchesleftLTR,
                              motifmismatchesrightLTR);
  }
  return haserr ? -1 : 0;
}

/*
 The following function searches for TSD and/or a specified palindromic motif
 at the borders of left LTR and the right LTR, respectively.
 */
static int gt_findcorrectboundaries(GtLTRharvestStream *lo,
                                     LTRboundaries *boundaries,
                                     const GtEncseq *encseq,
                                     GtError *err)
{
  unsigned int motifmismatchesleftLTR = 0,
               motifmismatchesrightLTR = 0;

  gt_error_check(err);
  gt_log_log("searching for correct boundaries in vicinity...\n");
  /* first: 5'-border of left LTR and 3'-border of right LTR */

  if (searchforTSDandorMotifoutside(lo,
                                    boundaries,
                                    encseq,
                                    &motifmismatchesleftLTR,
                                    &motifmismatchesrightLTR,
                                    err) != 0 )
  {
    return -1;
  }

  /* second: 3'-border of left LTR and 5'-border of right LTR */
  if ( lo->motif->allowedmismatches < (unsigned int)4 )
  {
    gt_log_log("second: searching for motif only around 3'border of left LTR "
               "and 5'-border of right LTR...\n");
    searchformotifonlyinside(lo,boundaries,encseq,
                             &motifmismatchesleftLTR,
                             &motifmismatchesrightLTR);
  }
  return 0;
}

static bool checklengthanddistanceconstraints(LTRboundaries *boundaries,
                                             RepeatInfo *repeatinfo)
{
  unsigned long ulen, vlen, dist_between_LTRs;

  ulen = boundaries->leftLTR_3  - boundaries->leftLTR_5  + 1;
  vlen = boundaries->rightLTR_3 - boundaries->rightLTR_5 + 1;
  dist_between_LTRs = boundaries->rightLTR_5 - boundaries->leftLTR_5;
  if (ulen > (unsigned long)repeatinfo->lmax
        || vlen > (unsigned long)repeatinfo->lmax
        || ulen < (unsigned long)repeatinfo->lmin
        || vlen < (unsigned long)repeatinfo->lmin
        || dist_between_LTRs > (unsigned long)repeatinfo->dmax
        || dist_between_LTRs < (unsigned long)repeatinfo->dmin
        || boundaries->leftLTR_3 >= boundaries->rightLTR_5)
  {
    boundaries->lengthdistconstraint = false;
    boundaries->similarity = 0.0;
    return false;
  } else
  {
    boundaries->lengthdistconstraint = true;
  }
  return true;
}

static void adjustboundariesfromXdropextension(Myxdropbest xdropbest_left,
                                               Myxdropbest xdropbest_right,
                                               unsigned long seed1_startpos,
                                               unsigned long seed2_startpos,
                                               unsigned long seed1_endpos,
                                               unsigned long seed2_endpos,
                                               LTRboundaries *boundaries)
{
  /* left alignment */
  boundaries->leftLTR_5  = seed1_startpos - xdropbest_left.ivalue;
  boundaries->rightLTR_5 = seed2_startpos - xdropbest_left.jvalue;

  /* right alignment */
  boundaries->leftLTR_3  = seed1_endpos + xdropbest_right.ivalue;
  boundaries->rightLTR_3 = seed2_endpos + xdropbest_right.jvalue;
}

/*
 The following function applies the filter algorithms one after another
 to all candidate pairs.
*/
static int gt_searchforLTRs(GtLTRharvestStream *lo,
                             GtArrayLTRboundaries *arrayLTRboundaries,
                             const GtEncseq *encseq,
                             GtError *err)
{
  unsigned long repeatcounter;
  GtArrayMyfrontvalue fronts;
  Myxdropbest xdropbest_left;
  Myxdropbest xdropbest_right;
  unsigned long alilen = 0,
         totallength,
         ulen,
         vlen,
         maxulen = 0,
         maxvlen = 0;
  GtUchar *useq = NULL,
        *vseq = NULL;
  unsigned long edist;
  Repeat *repeatptr;
  LTRboundaries *boundaries;
  bool haserr = false;

  gt_error_check(err);

  for (repeatcounter = 0; repeatcounter < lo->repeatinfo.repeats.nextfreeRepeat;
       repeatcounter++)
  {

    repeatptr = &(lo->repeatinfo.repeats.spaceRepeat[repeatcounter]);
    alilen = ((unsigned long)lo->repeatinfo.lmax) - repeatptr->len;

    /**** left (reverse) xdrop alignment ****/
    GT_INITARRAY (&fronts, Myfrontvalue);
    if (alilen <= repeatptr->pos1)
    {
      gt_evalxdroparbitscoresleft(&lo->arbitscores,
                               &xdropbest_left,
                               &fronts,
                               encseq,
                               encseq,
                               repeatptr->pos1,
                               repeatptr->pos1 + repeatptr->offset,
                               (int) alilen,
                               (int) alilen,
                               (Xdropscore)lo->xdropbelowscore);
    }
    else /* do not align over left sequence boundary */
    {
      gt_evalxdroparbitscoresleft(&lo->arbitscores,
                               &xdropbest_left,
                               &fronts,
                               encseq,
                               encseq,
                               repeatptr->pos1,
                               repeatptr->pos1 + repeatptr->offset,
                               (int) repeatptr->pos1,
                               (int) (repeatptr->pos1 + repeatptr->offset),
                               (Xdropscore)lo->xdropbelowscore);
    }
    GT_FREEARRAY (&fronts, Myfrontvalue);

    /**** right xdrop alignment ****/
    GT_INITARRAY (&fronts, Myfrontvalue);
    totallength = gt_encseq_total_length(encseq);
    if (alilen <= totallength - (repeatptr->pos1 + repeatptr->offset +
                                repeatptr->len) )
    {
      gt_evalxdroparbitscoresright (&lo->arbitscores,
                                 &xdropbest_right,
                                 &fronts,
                                 encseq,
                                 encseq,
                                 repeatptr->pos1 + repeatptr->len,
                                 repeatptr->pos1 + repeatptr->offset +
                                 repeatptr->len,
                                 (int) alilen,
                                 (int) alilen,
                                 lo->xdropbelowscore);
    }
    else /* do not align over right sequence boundary */
    {
      gt_evalxdroparbitscoresright(&lo->arbitscores,
                                &xdropbest_right,
                                &fronts,
                                encseq,
                                encseq,
                                repeatptr->pos1 + repeatptr->len,
                                repeatptr->pos1 + repeatptr->offset +
                                repeatptr->len,
                                (int) (totallength -
                                (repeatptr->pos1 + repeatptr->len)),
                                (int) (totallength -
                                (repeatptr->pos1 + repeatptr->offset +
                                 repeatptr->len)),
                                lo->xdropbelowscore);
    }
    GT_FREEARRAY (&fronts, Myfrontvalue);

    GT_GETNEXTFREEINARRAY(boundaries,arrayLTRboundaries,LTRboundaries,5);

    boundaries->contignumber = repeatptr->contignumber;
    boundaries->leftLTR_5 = (unsigned long) 0;
    boundaries->leftLTR_3 = (unsigned long) 0;
    boundaries->rightLTR_5 = (unsigned long) 0;
    boundaries->rightLTR_3 = (unsigned long) 0;
    boundaries->lenleftTSD = (unsigned long) 0;
    boundaries->lenrightTSD = (unsigned long) 0;
    boundaries->tsd = false;
    boundaries->motif_near_tsd = false;
    boundaries->motif_far_tsd = false;
    boundaries->skipped = false;
    boundaries->similarity = 0.0;

    /* store new boundaries-positions in boundaries */
    adjustboundariesfromXdropextension(
                   xdropbest_left,
                   xdropbest_right,
                   repeatptr->pos1,  /*seed1 startpos*/
                   repeatptr->pos1 + repeatptr->offset, /*seed2 startpos*/
                   repeatptr->pos1 + repeatptr->len - 1, /*seed1 endpos*/
                   repeatptr->pos1 + repeatptr->offset + repeatptr->len - 1,
                                                         /*seed2 endpos*/
                   boundaries);

    /* if search for motif and/or TSD */
    if ( lo->motif->allowedmismatches < 4U || lo->minlengthTSD > 1U)
    {
      if ( gt_findcorrectboundaries(lo, boundaries, encseq, err) != 0 )
      {
        haserr = true;
        break;
      }

      /* if search for TSDs and (not) motif */
      if ( boundaries->tsd &&
          (lo->motif->allowedmismatches >= (unsigned int)4 ||
          (boundaries->motif_near_tsd && boundaries->motif_far_tsd)) )
      {
        /* predicted as full LTR-pair, keep it */
      }
      else
      {
        /* if search for motif only (and not TSD) */
        if ( lo->minlengthTSD <= 1U &&
            boundaries->motif_near_tsd &&
            boundaries->motif_far_tsd )
        {
          /* predicted as full LTR-pair, keep it */
        }
        else
        {
          /* delete this LTR-pair candidate */
          arrayLTRboundaries->nextfreeLTRboundaries--;
          continue;
        }
      }
    }

    /* check length and distance constraints again */
    if (!checklengthanddistanceconstraints(boundaries, &lo->repeatinfo))
    {
      /* delete this LTR-pair candidate */
      arrayLTRboundaries->nextfreeLTRboundaries--;
      continue;
    }

    /* check similarity from candidate pair
       copy LTR sequences for greedyunitedist function */
    ulen = boundaries->leftLTR_3 - boundaries->leftLTR_5 + 1;
    vlen = boundaries->rightLTR_3 - boundaries->rightLTR_5 + 1;
    if (ulen > maxulen)
    {
      maxulen = ulen;
      ALLOCASSIGNSPACE(useq, useq, GtUchar, maxulen);
    }
    if (vlen > maxvlen)
    {
      maxvlen = vlen;
      ALLOCASSIGNSPACE(vseq, vseq, GtUchar, maxvlen);
    }

    gt_encseq_extract_encoded(encseq, useq,
                                         boundaries->leftLTR_5,
                                         boundaries->leftLTR_3);
    gt_encseq_extract_encoded(encseq, vseq,
                                         boundaries->rightLTR_5,
                                         boundaries->rightLTR_3);
    edist = greedyunitedist(useq,(unsigned long) ulen, /*Implement for encseq */
                            vseq,(unsigned long) vlen);

    /* determine similarity */
    boundaries->similarity = 100.0 * (1 - (((double) edist)/(MAX(ulen,vlen))));

    if (gt_double_smaller_double(boundaries->similarity,
                                 lo->similaritythreshold))
    {
      /* delete this LTR-pair candidate */
      arrayLTRboundaries->nextfreeLTRboundaries--;
    }

  }
  FREESPACE(useq);
  FREESPACE(vseq);

  return haserr ? -1 : 0;
}

/*
 The following function removes exact duplicates from the array of
 predicted LTR elements. Exact duplicates occur when
 different seeds are extended to same boundary coordinates.
 */
static void gt_removeduplicates(GtArrayLTRboundaries *arrayLTRboundaries)
{
  unsigned long i, j;
  unsigned long startpos_i, endpos_i, startpos_j, endpos_j;
  LTRboundaries *boundaries_i,
                *boundaries_j;

  for (i = 0; i < arrayLTRboundaries->nextfreeLTRboundaries; i++)
  {
    boundaries_i = &(arrayLTRboundaries->spaceLTRboundaries[i]);
    if ( boundaries_i->skipped )
    {
      continue;
    }
    startpos_i = boundaries_i->leftLTR_5;
    endpos_i   = boundaries_i->rightLTR_3;

    for (j = i + 1; j < arrayLTRboundaries->nextfreeLTRboundaries; j++)
    {
      boundaries_j = &(arrayLTRboundaries->spaceLTRboundaries[j]);
      if ( boundaries_j->skipped )
      {
        continue;
      }
      startpos_j = boundaries_j->leftLTR_5;
      endpos_j   = boundaries_j->rightLTR_3;

      if (startpos_i == startpos_j && endpos_i == endpos_j)
      {
        boundaries_j->skipped = true;
      }
    }
  }
}

/* The following function removes overlaps and deletes the prediction with
   a lower similarity value. If "nooverlapallowed" is set, then all
   overlapping predictions are deleted completely.
 */
static void gt_removeoverlapswithlowersimilarity(
  GtArrayLTRboundaries *arrayLTRboundaries,
  bool nooverlapallowed
  )
{
  unsigned long i, j;
  unsigned long startpos_i, endpos_i, startpos_j, endpos_j;
  LTRboundaries *boundaries_i,
                *boundaries_j;

  for (i = 0; i < arrayLTRboundaries->nextfreeLTRboundaries; i++)
  {
    boundaries_i = &(arrayLTRboundaries->spaceLTRboundaries[i]);
    if ( boundaries_i->skipped )
    {
      continue;
    }
    startpos_i = boundaries_i->leftLTR_5;
    endpos_i   = boundaries_i->rightLTR_3;

    for (j = i + 1; j < arrayLTRboundaries->nextfreeLTRboundaries; j++)
    {
      boundaries_j = &(arrayLTRboundaries->spaceLTRboundaries[j]);
      if ( boundaries_j->skipped )
      {
        continue;
      }
      startpos_j = boundaries_j->leftLTR_5;
      endpos_j   = boundaries_j->rightLTR_3;

      /* if overlap */
      if ( !((endpos_i < startpos_j) || (endpos_j < startpos_i)) )
      {
        if (nooverlapallowed)
        {
          /* All predictions in a cluster will be deleted. */

          /* take min(startpos_i, startpos_j) */
          if (startpos_j < startpos_i)
          {
            startpos_i = startpos_j;
          }
          /* take max(endpos_i, endpos_j) */
          if (endpos_i < endpos_j)
          {
            endpos_i = endpos_j;
          }
          /* delete both predictions */
          boundaries_i->skipped = true;
          boundaries_j->skipped = true;
        }
        else
        {
          /* take prediction with higher similarity */
          if ( boundaries_i->similarity >= boundaries_j->similarity )
          {
            boundaries_j->skipped = true;
          }
          else
          {
            boundaries_i->skipped = true;
            break;
          }
        }
      }
    }
  }
}

static int gt_ltrharvest_stream_next(GT_UNUSED GtNodeStream *ns,
                                     GtGenomeNode **gn,
                                     GtError *err)
{
  GtLTRharvestStream *ltrh_stream;
  int had_err = 0;
  gt_error_check(err);
  ltrh_stream = gt_ltrharvest_stream_cast(ns);

  /* LTRharvest run */
  if (ltrh_stream->state == GT_LTRHARVEST_STREAM_STATE_START) {

    GT_INITARRAY(&ltrh_stream->repeatinfo.repeats, Repeat);
    ltrh_stream->repeatinfo.encseq = ltrh_stream->encseq;
    ltrh_stream->prevseqnum = GT_UNDEF_ULONG;
    if (!had_err && gt_enumeratemaxpairs(ltrh_stream->ssar,
                      ltrh_stream->encseq,
                      gt_readmodeSequentialsuffixarrayreader(ltrh_stream->ssar),
                      (unsigned int) ltrh_stream->minseedlength,
                      gt_simpleexactselfmatchstore,
                      &ltrh_stream->repeatinfo,
                      err) != 0)
    {
      had_err = -1;
    }

    /* apply the filter algorithms */
    if (!had_err && gt_searchforLTRs(ltrh_stream,
                                      &ltrh_stream->arrayLTRboundaries,
                                      ltrh_stream->encseq,
                                      err) != 0)
    {
      had_err = -1;
    }

    /* not needed any longer */
    GT_FREEARRAY(&ltrh_stream->repeatinfo.repeats, Repeat);

    /* remove exact duplicates */
    if (!had_err)
    {
      gt_removeduplicates(&ltrh_stream->arrayLTRboundaries);
    }

    /* remove overlapping predictions if desired */
    if (!had_err && (ltrh_stream->nooverlaps || ltrh_stream->bestoverlaps))
    {
      gt_removeoverlapswithlowersimilarity(&ltrh_stream->arrayLTRboundaries,
                                            ltrh_stream->nooverlaps);
    }

    /* sort elements */
    if (!had_err)
    {
      ltrh_stream->bdptrtab =
        sortedltrboundaries(&ltrh_stream->numofboundaries,
                            &ltrh_stream->arrayLTRboundaries);
    }
    ltrh_stream->state = GT_LTRHARVEST_STREAM_STATE_REGIONS;
  }

  /* first stream out the region nodes */
  if (!had_err && ltrh_stream->state == GT_LTRHARVEST_STREAM_STATE_REGIONS) {
    bool skip = false;
    if (ltrh_stream->cur_elem_index < ltrh_stream->numofboundaries) {
      unsigned long seqnum, seqlength;
      GtGenomeNode *rn;
      GtStr *seqid;
      seqnum = ltrh_stream->bdptrtab[ltrh_stream->cur_elem_index]->contignumber;
      if (ltrh_stream->prevseqnum == GT_UNDEF_ULONG) {
        ltrh_stream->prevseqnum = seqnum;
      } else {
        while (ltrh_stream->prevseqnum == seqnum) {
          ltrh_stream->cur_elem_index++;
          if (ltrh_stream->cur_elem_index >= ltrh_stream->numofboundaries) {
            skip = true;
            break;
          }
          seqnum =
               ltrh_stream->bdptrtab[ltrh_stream->cur_elem_index]->contignumber;
        }
      }
      if (!skip) {
        ltrh_stream->prevseqnum = seqnum;
        seqlength = gt_encseq_seqlength(ltrh_stream->encseq, seqnum);
        seqid = gt_str_new_cstr("seq");
        gt_str_append_ulong(seqid, seqnum);
        rn = gt_region_node_new(seqid,
                                1 + (unsigned long) ltrh_stream->offset,
                                seqlength
                                  + (unsigned long) ltrh_stream->offset);
        gt_str_delete(seqid);
        *gn = rn;
        ltrh_stream->cur_elem_index++;
      } else {
        ltrh_stream->cur_elem_index = 0;
        ltrh_stream->state = GT_LTRHARVEST_STREAM_STATE_COMMENTS;
        *gn = NULL;
      }
    } else {
      ltrh_stream->cur_elem_index = 0;
      ltrh_stream->state = GT_LTRHARVEST_STREAM_STATE_COMMENTS;
      *gn = NULL;
    }
  }

  /* then stream out the comment nodes */
  if (!had_err && ltrh_stream->state == GT_LTRHARVEST_STREAM_STATE_COMMENTS) {
    bool skip = false;
    if (ltrh_stream->cur_elem_index < ltrh_stream->numofboundaries) {
      const char *desc;
      char buf[BUFSIZ];
      unsigned long desclen, seqnum;
      GtGenomeNode *cn;
      seqnum = ltrh_stream->bdptrtab[ltrh_stream->cur_elem_index]->contignumber;
      if (ltrh_stream->prevseqnum == GT_UNDEF_ULONG) {
        ltrh_stream->prevseqnum = seqnum;
      } else {
        while (ltrh_stream->prevseqnum == seqnum) {
          ltrh_stream->cur_elem_index++;
          if (ltrh_stream->cur_elem_index >= ltrh_stream->numofboundaries) {
            skip = true;
            break;
          }
          seqnum =
               ltrh_stream->bdptrtab[ltrh_stream->cur_elem_index]->contignumber;
        }
      }
      if (!skip) {
        ltrh_stream->prevseqnum = seqnum;
        desc = gt_encseq_description(ltrh_stream->encseq, &desclen, seqnum);
        (void) strncpy(buf, desc, (size_t) (desclen * sizeof (char)));
        buf[desclen] = '\0';
        cn = gt_comment_node_new(buf);
        *gn = cn;
        ltrh_stream->cur_elem_index++;
      } else {
        ltrh_stream->cur_elem_index = 0;
        ltrh_stream->state = GT_LTRHARVEST_STREAM_STATE_FEATURES;
        *gn = NULL;
      }
    } else {
      ltrh_stream->cur_elem_index = 0;
      ltrh_stream->state = GT_LTRHARVEST_STREAM_STATE_FEATURES;
      *gn = NULL;
    }
  }

  /* finally, stream out the features */
  if (!had_err && ltrh_stream->state == GT_LTRHARVEST_STREAM_STATE_FEATURES) {
    if (ltrh_stream->cur_elem_index <  ltrh_stream->numofboundaries) {
      GtStr *seqid, *source;
      char buf[BUFSIZ];
      GtGenomeNode *node, *parent;
      const LTRboundaries *elem =
        (ltrh_stream->bdptrtab[ltrh_stream->cur_elem_index]);
      unsigned long seqstartpos;

      seqstartpos = gt_encseq_seqstartpos(ltrh_stream->encseq,
                                          elem->contignumber);
      seqid = gt_str_new_cstr("seq");
      source = gt_str_new_cstr(GT_LTRHARVEST_NAME);

      gt_str_append_ulong(seqid, elem->contignumber);

      /* repeat_region */
      node = gt_feature_node_new(seqid,
                                 "repeat_region",
                                 elem->leftLTR_5 - seqstartpos + 1
                                   - elem->lenleftTSD
                                   + (unsigned long) ltrh_stream->offset,
                                 elem->rightLTR_3 - seqstartpos + 1
                                   + elem->lenrightTSD
                                   + (unsigned long) ltrh_stream->offset,
                                 GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*) node, source);
      *gn = node;
      parent = node;
      if (ltrh_stream->motif->allowedmismatches < 4U) {
        node = gt_feature_node_new(seqid,
                                   "inverted_repeat",
                                   elem->leftLTR_5 - seqstartpos + 1
                                     + (unsigned long) ltrh_stream->offset,
                                   elem->leftLTR_5 - seqstartpos + 2
                                     + (unsigned long) ltrh_stream->offset,
                                   GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                                  (GtFeatureNode*) node);
        node = gt_feature_node_new(seqid,
                                   "inverted_repeat",
                                   elem->leftLTR_3 - seqstartpos
                                     + (unsigned long) ltrh_stream->offset,
                                   elem->leftLTR_3 - seqstartpos + 1
                                     + (unsigned long) ltrh_stream->offset,
                                   GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                                  (GtFeatureNode*) node);
        node = gt_feature_node_new(seqid,
                                   "inverted_repeat",
                                   elem->rightLTR_5 - seqstartpos + 1
                                     + (unsigned long) ltrh_stream->offset,
                                   elem->rightLTR_5 - seqstartpos + 2
                                     + (unsigned long) ltrh_stream->offset,
                                   GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                                  (GtFeatureNode*) node);
        node = gt_feature_node_new(seqid,
                                   "inverted_repeat",
                                   elem->rightLTR_3 - seqstartpos
                                     + (unsigned long) ltrh_stream->offset,
                                   elem->rightLTR_3 - seqstartpos + 1
                                     + (unsigned long) ltrh_stream->offset,
                                   GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                                  (GtFeatureNode*) node);
      }
      /* target_site_duplication */
      if (ltrh_stream->minlengthTSD > 1U) {
        node = gt_feature_node_new(seqid,
                                   "target_site_duplication",
                                   elem->leftLTR_5 - seqstartpos + 1
                                     - elem->lenleftTSD
                                     + (unsigned long) ltrh_stream->offset,
                                   elem->leftLTR_5 - seqstartpos
                                     + (unsigned long) ltrh_stream->offset,
                                   GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                                  (GtFeatureNode*) node);
        node = gt_feature_node_new(seqid,
                                   "target_site_duplication",
                                   elem->rightLTR_3 - seqstartpos + 2
                                     + (unsigned long) ltrh_stream->offset,
                                   elem->rightLTR_3 - seqstartpos + 1
                                     + elem->lenrightTSD
                                     + (unsigned long) ltrh_stream->offset,
                                   GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                                  (GtFeatureNode*) node);
      }
      /* LTR_retrotransposon */
      node = gt_feature_node_new(seqid,
                                 "LTR_retrotransposon",
                                 elem->leftLTR_5 - seqstartpos + 1
                                   + (unsigned long) ltrh_stream->offset,
                                 elem->rightLTR_3 - seqstartpos + 1
                                   + (unsigned long) ltrh_stream->offset,
                                 GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*) node, source);
      (void) snprintf(buf, BUFSIZ-1, "%.2f", elem->similarity);
      gt_feature_node_set_attribute((GtFeatureNode*) node,
                                    "LTR_similarity",
                                    buf);
      (void) snprintf(buf, BUFSIZ-1, "%lu", elem->contignumber);
      gt_feature_node_set_attribute((GtFeatureNode*) node,
                                    "seq_number",
                                    buf);
      gt_feature_node_add_child((GtFeatureNode*) parent,
                                (GtFeatureNode*) node);
      parent = node;
      /* long_terminal_repeat */
      node = gt_feature_node_new(seqid,
                                 "long_terminal_repeat",
                                 elem->leftLTR_5 - seqstartpos + 1
                                   + (unsigned long) ltrh_stream->offset,
                                 elem->leftLTR_3 - seqstartpos + 1
                                   + (unsigned long) ltrh_stream->offset,
                                 GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*) node, source);
      gt_feature_node_add_child((GtFeatureNode*) parent,
                                (GtFeatureNode*) node);
      node = gt_feature_node_new(seqid,
                                 "long_terminal_repeat",
                                 elem->rightLTR_5 - seqstartpos + 1
                                   + (unsigned long) ltrh_stream->offset,
                                 elem->rightLTR_3 - seqstartpos + 1
                                   + (unsigned long) ltrh_stream->offset,
                                 GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*) node, source);
      gt_feature_node_add_child((GtFeatureNode*) parent,
                                (GtFeatureNode*) node);
      gt_str_delete(seqid);
      gt_str_delete(source);
      ltrh_stream->cur_elem_index++;
    } else {
      *gn = NULL;
    }
  }
  return had_err;
}

static void gt_ltrharvest_stream_free(GtNodeStream *ns)
{
  GtLTRharvestStream *ltrh_stream = gt_ltrharvest_stream_cast(ns);
  GT_FREEARRAY(&ltrh_stream->arrayLTRboundaries, LTRboundaries);
  if (ltrh_stream->ssar != NULL)
    gt_freeSequentialsuffixarrayreader(&ltrh_stream->ssar);
  if (ltrh_stream->bdptrtab != NULL)
    gt_free(ltrh_stream->bdptrtab);
}

const GtNodeStreamClass* gt_ltrharvest_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtLTRharvestStream),
                                   gt_ltrharvest_stream_free,
                                   gt_ltrharvest_stream_next);
  }
  return nsc;
}

const GtEncseq* gt_ltrharvest_stream_get_encseq(GtNodeStream *ltrh_stream)
{
  GtLTRharvestStream *stream;
  gt_assert(ltrh_stream != NULL);
  stream = gt_ltrharvest_stream_cast(ltrh_stream);
  return stream->encseq;
}

GtNodeStream* gt_ltrharvest_stream_new(GtStr *str_indexname,
                                       GtRange searchrange,
                                       unsigned long minseedlength,
                                       unsigned long minltrlength,
                                       unsigned long maxltrlength,
                                       unsigned long mindistance,
                                       unsigned long maxdistance,
                                       double similaritythreshold,
                                       int xdropbelowscore,
                                       Arbitraryscores arbitscores,
                                       GtLTRFourCharMotif *motif,
                                       bool verbosemode,
                                       bool nooverlaps,
                                       bool bestoverlaps,
                                       unsigned long offset,
                                       unsigned int minlengthTSD,
                                       unsigned int maxlengthTSD,
                                       unsigned long vicinity,
                                       GtError *err)
{
  int had_err = 0;
  GtNodeStream *ns = gt_node_stream_create(gt_ltrharvest_stream_class(), false);
  GtLTRharvestStream *ltrh_stream = gt_ltrharvest_stream_cast(ns);

  ltrh_stream->str_indexname = str_indexname;
  ltrh_stream->minseedlength = minseedlength;
  ltrh_stream->similaritythreshold = similaritythreshold;
  ltrh_stream->xdropbelowscore = xdropbelowscore;
  ltrh_stream->arbitscores = arbitscores;
  ltrh_stream->motif = motif;
  ltrh_stream->verbosemode = verbosemode;
  ltrh_stream->offset = offset;
  ltrh_stream->repeatinfo.dmin = mindistance;
  ltrh_stream->repeatinfo.dmax = maxdistance;
  ltrh_stream->repeatinfo.lmin = minltrlength;
  ltrh_stream->repeatinfo.lmax = maxltrlength;
  ltrh_stream->repeatinfo.ltrsearchseqrange = searchrange;
  ltrh_stream->minlengthTSD = minlengthTSD;
  ltrh_stream->maxlengthTSD = maxlengthTSD;
  ltrh_stream->nooverlaps = nooverlaps;
  ltrh_stream->bestoverlaps = bestoverlaps;
  ltrh_stream->vicinityforcorrectboundaries = vicinity;
  ltrh_stream->cur_elem_index = 0;
  ltrh_stream->state = GT_LTRHARVEST_STREAM_STATE_START;
  /* init array for maximal repeats */
  GT_INITARRAY(&ltrh_stream->arrayLTRboundaries, LTRboundaries);

  ltrh_stream->ssar =
    gt_newSequentialsuffixarrayreaderfromfile(gt_str_get(str_indexname),
                                                SARR_LCPTAB | SARR_SUFTAB |
                                                SARR_ESQTAB | SARR_DESTAB |
                                                SARR_SSPTAB | SARR_SDSTAB,
                                                SEQ_mappedboth,
                                                NULL,
                                                err);
  if (ltrh_stream->ssar == NULL)
  {
    gt_node_stream_delete(ns);
    return NULL;
  }
  /* get encseq associated with suffix array */
  ltrh_stream->encseq = gt_encseqSequentialsuffixarrayreader(ltrh_stream->ssar);

  /* encode motif according to encseq alphabet */
  had_err = gt_ltr_four_char_motif_encode(motif, ltrh_stream->encseq, err);
  if (had_err) {
    gt_node_stream_delete(ns);
    return NULL;
  }

  return had_err ? NULL : ns;
}
