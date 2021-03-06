/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "core/types_api.h"
#include "core/readmode.h"

void gt_inplace_reversecomplement(GtUchar *seq,unsigned long len)
{
  GtUchar tmp, *frontptr, *backptr;

  for (frontptr = seq, backptr = seq + len - 1;
       frontptr < backptr; frontptr++, backptr--)
  {
    tmp = *frontptr;
    *frontptr = GT_COMPLEMENTBASE(*backptr);
    *backptr = GT_COMPLEMENTBASE(tmp);
  }
}

void gt_copy_reversecomplement(GtUchar *dest,const GtUchar *src,
                               unsigned long len)
{
  GtUchar *destptr;
  const GtUchar *srcptr;

  for (destptr = dest, srcptr = src + len - 1;
       destptr < dest + len; destptr++, srcptr--)
  {
    *destptr = GT_COMPLEMENTBASE(*srcptr);
  }
}
