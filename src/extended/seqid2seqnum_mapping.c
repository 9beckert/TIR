/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include <string.h>
#include "core/cstr_api.h"
#include "core/hashmap.h"
#include "core/ma_api.h"
#include "core/parseutils.h"
#include "core/undef_api.h"
#include "extended/seqid2seqnum_mapping.h"

typedef GtArray SeqidInfo;

typedef struct {
  unsigned long seqnum;
  GtRange descrange;
} SeqidInfoElem;

static SeqidInfo* seqid_info_new(unsigned long seqnum, const GtRange *descrange)
{
  SeqidInfoElem seqid_info_elem;
  GtArray *seqid_info;
  gt_assert(descrange);
  seqid_info = gt_array_new(sizeof (SeqidInfoElem));
  seqid_info_elem.seqnum = seqnum;
  seqid_info_elem.descrange = *descrange;
  gt_array_add(seqid_info, seqid_info_elem);
  return seqid_info;
}

static void seqid_info_delete(SeqidInfo *seqid_info)
{
  if (!seqid_info) return;
  gt_array_delete(seqid_info);
}

static int seqid_info_add(SeqidInfo *seqid_info, unsigned long seqnum,
                          const GtRange *range, const char *filename,
                          const char *seqid, GtError *err)
{
  SeqidInfoElem *seqid_info_elem_ptr, seqid_info_elem;
  gt_error_check(err);
  gt_assert(range);
  seqid_info_elem_ptr = gt_array_get_first(seqid_info);
  if (range->end == GT_UNDEF_ULONG ||
      seqid_info_elem_ptr->descrange.end == GT_UNDEF_ULONG) {
    gt_error_set(err, "sequence file \"%s\" does contain multiple sequences "
                  "with ID \"%s\" and not all of them have description ranges",
                  filename, seqid);
    return -1;
  }
  seqid_info_elem.seqnum = seqnum;
  seqid_info_elem.descrange = *range;
  gt_array_add(seqid_info, seqid_info_elem);
  return 0;
}

static int seqid_info_get(SeqidInfo *seqid_info, unsigned long *seqnum,
                          GtRange *outrange, const GtRange *inrange,
                          const char *filename, const char *seqid, GtError *err)
{
  SeqidInfoElem *seqid_info_elem;
  unsigned long i;
  gt_error_check(err);
  gt_assert(seqid_info && seqnum && outrange && inrange);
  for (i = 0; i < gt_array_size(seqid_info); i++) {
    seqid_info_elem = gt_array_get(seqid_info, i);
    if (seqid_info_elem->descrange.end == GT_UNDEF_ULONG ||
        gt_range_contains(&seqid_info_elem->descrange, inrange)) {
      *seqnum = seqid_info_elem->seqnum;
      *outrange = seqid_info_elem->descrange;
      return 0;
    }
  }
  gt_error_set(err, "cannot find sequence ID \"%s\" (with range %lu,%lu) in "
               "sequence file \"%s\"", seqid, inrange->start, inrange->end,
               filename);
  return -1;
}

struct GtSeqid2SeqnumMapping {
  char *filename;
  GtHashmap *map;
  const char *cached_seqid;
  unsigned long cached_seqnum;
  GtRange cached_range;
};

static int fill_mapping(GtSeqid2SeqnumMapping *mapping, GtBioseq *bioseq,
                        GtError *err)
{
  SeqidInfo *seqid_info;
  unsigned long i, j;
  GtRange descrange;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(mapping && bioseq);
  for (i = 0; !had_err && i < gt_bioseq_number_of_sequences(bioseq); i++) {
    const char *desc = gt_bioseq_get_description(bioseq, i);
    if (gt_parse_description_range(desc, &descrange)) {
      /* no offset could be parsed -> store description as sequence id */
      descrange.start = 1;
      descrange.end = GT_UNDEF_ULONG;
      if ((seqid_info = gt_hashmap_get(mapping->map, desc))) {
        had_err = seqid_info_add(seqid_info, i, &descrange, mapping->filename,
                                 desc, err);
        gt_assert(had_err); /* adding a seqid without range should fail */
      }
      else {
        seqid_info = seqid_info_new(i, &descrange);
        gt_hashmap_add(mapping->map, gt_cstr_dup(desc), seqid_info);
      }
    }
    else {
      char *dup;
      /* offset could be parsed -> store description up to ':' as sequence id */
      j = 0;
      while (desc[j] != ':')
        j++;
      dup = gt_malloc((j + 1) * sizeof *dup);
      strncpy(dup, desc, j);
      dup[j] = '\0';
      if ((seqid_info = gt_hashmap_get(mapping->map, dup))) {
        had_err = seqid_info_add(seqid_info, i, &descrange, mapping->filename,
                                 dup, err);
        gt_free(dup);
      }
      else {
        seqid_info = seqid_info_new(i, &descrange);
        gt_hashmap_add(mapping->map, dup, seqid_info);
      }
    }
  }
  return had_err;
}

GtSeqid2SeqnumMapping* gt_seqid2seqnum_mapping_new(GtBioseq *bioseq,
                                                   GtError *err)
{
  GtSeqid2SeqnumMapping *mapping;
  gt_error_check(err);
  gt_assert(bioseq);
  mapping = gt_malloc(sizeof *mapping);
  mapping->filename = gt_cstr_dup(gt_bioseq_filename(bioseq));
  mapping->map = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                (GtFree) seqid_info_delete);
  if (fill_mapping(mapping, bioseq, err)) {
    gt_seqid2seqnum_mapping_delete(mapping);
    return NULL;
  }
  mapping->cached_seqid = NULL;
  return mapping;
}

void gt_seqid2seqnum_mapping_delete(GtSeqid2SeqnumMapping *mapping)
{
  if (!mapping) return;
  gt_hashmap_delete(mapping->map);
  gt_free(mapping->filename);
  gt_free(mapping);
}

int gt_seqid2seqnum_mapping_map(GtSeqid2SeqnumMapping *mapping,
                                const char *seqid, const GtRange *inrange,
                                unsigned long *seqnum, unsigned long *offset,
                                GtError *err)
{
  SeqidInfo *seqid_info;
  GtRange outrange;
  gt_error_check(err);
  gt_assert(mapping && seqid && seqnum && offset);
  /* try to answer request from cache */
  if (mapping->cached_seqid && !strcmp(seqid, mapping->cached_seqid) &&
      (mapping->cached_range.end == GT_UNDEF_ULONG ||
       gt_range_contains(&mapping->cached_range, inrange))) {
    *seqnum = mapping->cached_seqnum;
    *offset = mapping->cached_range.start;
    return 0;
  }
  /* cache miss -> regular mapping */
  if (!(seqid_info = gt_hashmap_get(mapping->map, seqid))) {
    gt_error_set(err, "sequence file \"%s\" does not contain a sequence with "
                      "ID \"%s\"", mapping->filename, seqid);
    return -1;
  }
  /* get results from seqid info */
  if (seqid_info_get(seqid_info, seqnum, &outrange, inrange, mapping->filename,
                     seqid, err)) {
    return -1;
  }
  /* report offset */
  *offset = outrange.start;
  /* store result in cache */
  mapping->cached_seqid = gt_hashmap_get_key(mapping->map, seqid);
  gt_assert(mapping->cached_seqid);
  mapping->cached_seqnum = *seqnum;
  mapping->cached_range = outrange;
  return 0;
}
