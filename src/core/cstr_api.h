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

#ifndef CSTR_API_H
#define CSTR_API_H

#include <stdio.h>

/* Cstr module */

/* Creates a duplicate of string <cstr> using the GenomeTools memory
   allocator. */
char*         gt_cstr_dup(const char *cstr);

/* Creates a duplicate of string <cstr> using the GenomeTools memory allocator.
   The string needs not be \0-terminated, instead its <length> must be given. */
char*         gt_cstr_dup_nt(const char *cstr, unsigned long length);

/* Replace each occurence of <f> in <cstr> to <t>. */
void          gt_cstr_rep(char *cstr, char f, char t);

/* Outputs the first <length> characters of the string <cstr> to file pointer
   <outfp>. */
void          gt_cstr_show(const char *cstr, unsigned long length, FILE *outfp);

/* Returns the length of the prefix of <cstr> ending just before <c>, if <cstr>
   does not contain <c>, strlen(cstr) is returned. */
unsigned long gt_cstr_length_up_to_char(const char *cstr, char c);

/* Removes all occurrences of <remove> from the right end of <cstr>. */
char*         gt_cstr_rtrim(char* cstr, char remove);

#endif
