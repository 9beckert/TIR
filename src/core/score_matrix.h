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

#ifndef SCORE_MATRIX_H
#define SCORE_MATRIX_H

#include "core/alpha.h"
#include "core/error.h"

typedef struct GT_ScoreMatrix GT_ScoreMatrix;

/* a score matrix is always defined over a given alphabet */
GT_ScoreMatrix* gt_score_matrix_new(GT_Alpha*);
/* reads in a protein scorematrix from the given <path> and returns it */
GT_ScoreMatrix* gt_score_matrix_new_read_protein(const char *path, GtError*);
unsigned int gt_score_matrix_get_dimension(const GT_ScoreMatrix*);
int          gt_score_matrix_get_score(const GT_ScoreMatrix*,
                                    unsigned int, unsigned int);
void         gt_score_matrix_set_score(GT_ScoreMatrix*,
                                    unsigned int, unsigned int, int);
const int**  gt_score_matrix_get_scores(const GT_ScoreMatrix*);
void         gt_score_matrix_show(const GT_ScoreMatrix*, FILE*);
void         gt_score_matrix_delete(GT_ScoreMatrix*);

#endif
