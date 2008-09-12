/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include "extended/cds_stream.h"
#include "extended/cds_visitor.h"
#include "extended/genome_stream_rep.h"

struct CDSStream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *cds_visitor;
};

#define cds_stream_cast(GS)\
        genome_stream_cast(cds_stream_class(), GS)

static int cds_stream_next_tree(GenomeStream *gs, GT_GenomeNode **gn,
                                GtError *err)
{
  CDSStream *cds_stream;
  int had_err;
  gt_error_check(err);
  cds_stream = cds_stream_cast(gs);
  had_err = genome_stream_next_tree(cds_stream->in_stream, gn, err);
  if (!had_err && *gn)
    had_err = gt_genome_node_accept(*gn, cds_stream->cds_visitor, err);
  return had_err;
}

static void cds_stream_free(GenomeStream *gs)
{
  CDSStream *cds_stream = cds_stream_cast(gs);
  genome_visitor_delete(cds_stream->cds_visitor);
  genome_stream_delete(cds_stream->in_stream);
}

const GenomeStreamClass* cds_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (CDSStream),
                                         cds_stream_next_tree,
                                         cds_stream_free };
  return &gsc;
}

GenomeStream* cds_stream_new(GenomeStream *in_stream, RegionMapping *rm,
                             const char *source)
{
  GenomeStream *gs;
  CDSStream *cds_stream;
  GtStr *source_str;
  int had_err = 0;
  gs = genome_stream_create(cds_stream_class(), true);
  cds_stream = cds_stream_cast(gs);
  source_str = gt_str_new_cstr(source);
  cds_stream->in_stream = genome_stream_ref(in_stream);
  cds_stream->cds_visitor = cds_visitor_new(rm, source_str);
  if (!cds_stream->cds_visitor)
    had_err = -1;
  gt_str_delete(source_str);
  if (had_err) {
    cds_stream_free(gs);
    return NULL;
  }
  return gs;
}
