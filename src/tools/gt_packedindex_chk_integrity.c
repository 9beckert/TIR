/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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

/**
 * \file gt_packedindex_chk_integrity
 * calls routines to validate basic index integrity
 */

#include "gt_packedindex_chk_integrity.h"
#include "core/error.h"
#include "core/option.h"
#include "core/versionfunc.h"
#include "match/eis-encidxseq.h"
#include "match/eis-encidxseq-param.h"
#include "match/eis-encidxseq-construct.h"
#include "match/verbose-def.h"
#include "tools/gt_packedindex_chk_integrity.h"

#define DEFAULT_PROGRESS_INTERVAL  100000UL

struct chkIndexOptions
{
  unsigned long skipCount;
  unsigned long progressInterval;
  int checkFlags;
  bool verboseOutput;
  enum seqBaseEncoding encType;
  int EISFeatureSet;
};

static OPrval
parseChkIndexOptions(int *parsed_args, int argc, const char *argv[],
                     struct chkIndexOptions *param, GtError *err);

extern int
gt_packedindex_chk_integrity(int argc, const char *argv[], GtError *err)
{
  struct encIdxSeq *seq;
  struct chkIndexOptions params;
  GtStr *inputProject;
  int parsedArgs;
  int had_err = 0;
  Verboseinfo *verbosity = NULL;
  gt_error_check(err);

  switch (parseChkIndexOptions(&parsedArgs, argc, argv, &params, err))
  {
    case OPTIONPARSER_OK:
      break;
    case OPTIONPARSER_ERROR:
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      return 0;
  }

  inputProject = gt_str_new_cstr(argv[parsedArgs]);

  verbosity = newverboseinfo(params.verboseOutput);

  seq = loadEncIdxSeq(inputProject, params.encType, params.EISFeatureSet,
                      verbosity, err);
  if ((had_err = seq == NULL))
  {
    gt_error_set(err, "Failed to load index: %s", gt_str_get(inputProject));
  }
  else
  {
    fprintf(stderr, "# Using index over sequence "FormatSeqpos
            " symbols long.\n", EISLength(seq));
    {
      int corrupt
        = EISVerifyIntegrity(seq, inputProject, params.skipCount,
                             params.progressInterval, stderr,
                             params.checkFlags, verbosity, err);
      if ((had_err = corrupt != 0))
      {
        fputs(gt_error_get(err), stderr); fputs("\n", stderr);
        gt_error_set(err, "Integrity check failed for index: %s",
                  EISIntegrityCheckResultStrings[corrupt]);
      }
    }
  }
  if (seq) deleteEncIdxSeq(seq);
  if (inputProject) gt_str_delete(inputProject);
  if (verbosity) freeverboseinfo(&verbosity);
  return had_err?-1:0;
}

static OPrval
parseChkIndexOptions(int *parsed_args, int argc, const char *argv[],
                     struct chkIndexOptions *params, GtError *err)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  bool extRankCheck;

  gt_error_check(err);
  op = option_parser_new("indexname",
                         "Map <indexname> block composition index"
                         "and bwt and check index integrity.");

  option = option_new_bool("v",
                           "print verbose progress information",
                           &params->verboseOutput,
                           false);
  option_parser_add_option(op, option);

  option = option_new_ulong("skip", "number of symbols to skip",
                            &params->skipCount, 0);
  option_parser_add_option(op, option);

  option = option_new_ulong("ticks", "print dot after this many symbols"
                            " tested okay", &params->progressInterval,
                            DEFAULT_PROGRESS_INTERVAL);
  option_parser_add_option(op, option);

  option = option_new_bool("ext-rank-check",
                           "do additional checks of rank query results",
                           &extRankCheck, false);
  option_parser_add_option(op, option);

  option_parser_set_min_max_args(op, 1, 1);
  oprval = option_parser_parse(op, parsed_args, argc, (const char**) argv,
                               versionfunc, err);
  option_parser_delete(op);
  params->checkFlags = EIS_VERIFY_BASIC | (extRankCheck?EIS_VERIFY_EXT_RANK:0);
  params->EISFeatureSet = EIS_FEATURE_REGION_SUMS;
  params->encType = BWT_ON_BLOCK_ENC;
  return oprval;
}
