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

#ifndef FILEUTILS_H
#define FILEUTILS_H

#include "core/fileutils_api.h"

/* Returns true if the file with the name composed of the concatenation of
   <path> and <suffix> exists, false otherwise. */
bool           gt_file_with_suffix_exists(const char *path, const char *suffix);

/* Returns the size of the file whose name name is composed of the
  concatenation of <path> and <suffix>. */
off_t          gt_file_with_suffix_size(const char *path, const char *suffix);

/* Return the size of <file>. */
off_t          gt_file_size(const char *file);

/* Compare two files bytewise */

void gt_xfile_cmp(const char *file1,const char *file2);

#endif
