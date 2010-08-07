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

#ifndef SEQIDS_TO_MD5_VISITOR_H
#define SEQIDS_TO_MD5_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct GtSeqidsToMD5Visitor GtSeqidsToMD5Visitor;

#include "extended/node_visitor.h"
#include "extended/region_mapping.h"

const GtNodeVisitorClass* gt_seqids_to_md5_visitor_class(void);
/* Takes ownership of <region_mapping>. */
GtNodeVisitor*            gt_seqids_to_md5_visitor_new(GtRegionMapping
                                                       *region_mapping);

#endif