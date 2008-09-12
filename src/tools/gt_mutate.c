/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/bioseq_iterator.h"
#include "core/fasta.h"
#include "core/ma.h"
#include "core/option.h"
#include "extended/gtdatahelp.h"
#include "extended/mutate.h"
#include "tools/gt_mutate.h"

typedef struct {
  unsigned int rate; /* the mutate rate */
} MutateArguments;

static void* gt_mutate_arguments_new(void)
{
  return gt_calloc(1, sizeof (MutateArguments));
}

static void gt_mutate_arguments_delete(void *tool_arguments)
{
  MutateArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_free(arguments);
}

static OptionParser* gt_mutate_option_parser_new(void *tool_arguments)
{
  MutateArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *o;
  assert(arguments);
  op = option_parser_new("[option ...] [sequence_file ...]",
                         "Mutate the sequences of the given sequence_file(s) "
                         "and show them on stdout.");
  /* -rate */
  o = option_new_uint_max("rate", "set the mutation rate", &arguments->rate, 1,
                          100);
  option_parser_add_option(op, o);

  /* parse */
  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  return op;
}

static int gt_mutate_runner(int argc, const char **argv, int parsed_args,
                            void *tool_arguments, GtError *err)
{
  MutateArguments *arguments = tool_arguments;
  GT_BioseqIterator *bsi;
  unsigned long i;
  GT_Bioseq *bioseq;
  Seq *mutated_seq;
  int had_err;

  gt_error_check(err);
  assert(arguments);

  bsi = gt_bioseq_iterator_new(argc - parsed_args, argv + parsed_args);

  while (!(had_err = gt_bioseq_iterator_next(bsi, &bioseq, err)) && bioseq) {
    for (i = 0; i < gt_bioseq_number_of_sequences(bioseq); i++) {
      mutated_seq = mutate(gt_bioseq_get_description(bioseq, i),
                           gt_bioseq_get_sequence(bioseq, i),
                           gt_bioseq_get_sequence_length(bioseq, i),
                           gt_bioseq_get_alpha(bioseq), arguments->rate);
      gt_fasta_show_entry(seq_get_description(mutated_seq),
                          seq_get_orig(mutated_seq), seq_length(mutated_seq),
                          0);
      seq_delete(mutated_seq);
    }
    gt_bioseq_delete(bioseq);
  }

  gt_bioseq_iterator_delete(bsi);

  return had_err;
}

Tool* gt_mutate(void)
{
  return tool_new(gt_mutate_arguments_new,
                  gt_mutate_arguments_delete,
                  gt_mutate_option_parser_new,
                  NULL,
                  gt_mutate_runner);
}
