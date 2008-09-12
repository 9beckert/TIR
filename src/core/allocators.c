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

#include "core/allocators.h"
#include "core/cstr.h"
#include "core/cstr_array.h"
#include "core/fa.h"
#include "core/ma.h"
#include "core/option.h"
#include "core/splitter.h"
#include "core/symbol.h"
#include "core/versionfunc.h"
#include "core/warning.h"
#include "core/xansi.h"

static bool spacepeak = false;

static OPrval parse_env_options(int argc, const char **argv, GtError *err)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  op = option_parser_new("GT_ENV_OPTIONS='[option ...]' ...",
                         "Parse the options contained in the "
                         "environment variable GT_ENV_OPTIONS.");
  o = option_new_bool("spacepeak", "show space peak on stdout upon deletion",
                      &spacepeak, false);
  option_parser_add_option(op, o);
  option_parser_set_max_args(op, 0);
  oprval = option_parser_parse(op, NULL, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

static void proc_gt_env_options(void)
{
  int argc;
  char *env_options, **argv;
  Splitter *splitter;
  GtError *err;
  /* construct argument vector from $GT_ENV_OPTIONS */
  env_options = getenv("GT_ENV_OPTIONS");
  if (!env_options)
    return;
  env_options = gt_cstr_dup(env_options); /* make writeable copy */
  splitter = splitter_new();
  splitter_split(splitter, env_options, strlen(env_options), ' ');
  argc = splitter_size(splitter);
  argv = gt_cstr_array_preprend((const char**) splitter_get_tokens(splitter),
                             "env");
  argc++;
  /* parse options contained in $GT_ENV_OPTIONS */
  err = gt_error_new();
  switch (parse_env_options(argc, (const char**) argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      fprintf(stderr, "error parsing $GT_ENV_OPTIONS: %s\n", gt_error_get(err));
      gt_error_unset(err);
      break;
    case OPTIONPARSER_REQUESTS_EXIT: break;
  }
  gt_error_delete(err);
  gt_free(env_options);
  splitter_delete(splitter);
  gt_cstr_array_delete(argv);
}

void gt_allocators_init(void)
{
  const char *bookkeeping;
  bookkeeping = getenv("GT_MEM_BOOKKEEPING");
  gt_ma_init(bookkeeping && !strcmp(bookkeeping, "on"));
  proc_gt_env_options();
  if (spacepeak && !(bookkeeping && !strcmp(bookkeeping, "on")))
    warning("GT_ENV_OPTIONS=-spacepeak used without GT_MEM_BOOKKEEPING=on");
}

static void gt_allocators_atexit_func(void)
{
  (void) gt_allocators_clean();
}

void gt_allocators_reg_atexit_func(void)
{
  xatexit(gt_allocators_atexit_func);
}

int gt_allocators_clean(void)
{
  int fa_fptr_rval, fa_mmap_rval, gt_rval;
  if (spacepeak) {
    gt_ma_show_space_peak(stdout);
    gt_fa_show_space_peak(stdout);
  }
  fa_fptr_rval = gt_fa_check_fptr_leak();
  fa_mmap_rval = gt_fa_check_mmap_leak();
  gt_fa_clean();
  gt_symbol_clean();
  gt_rval = gt_ma_check_space_leak();
  gt_ma_clean();
  return fa_fptr_rval || fa_mmap_rval || gt_rval;
}
