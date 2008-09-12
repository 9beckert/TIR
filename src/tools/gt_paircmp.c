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

#include <string.h>
#include <inttypes.h>
#include "core/error.h"
#include "core/versionfunc.h"
#include "core/option.h"
#include "core/str.h"
#include "core/symboldef.h"
#include "match/spacedef.h"
#include "match/test-pairwise.h"
#include "tools/gt_paircmp.h"

typedef struct
{
  GtStr *charlist;
  unsigned long len;
} Charlistlen;

typedef struct
{
  GtStrArray *strings,
           *files;
  Charlistlen *charlistlen;
  GtStr *text;
} Cmppairwiseopt;

static void showsimpleoptions(const Cmppairwiseopt *opt)
{
  if (gt_strarray_size(opt->strings) > 0)
  {
    printf("# two strings \"%s\" \"%s\"\n",gt_strarray_get(opt->strings,0),
                                           gt_strarray_get(opt->strings,1UL));
    return;
  }
  if (gt_strarray_size(opt->files) > 0)
  {
    printf("# two files \"%s\" \"%s\"\n",gt_strarray_get(opt->files,0),
                                         gt_strarray_get(opt->files,1UL));
    return;
  }
  if (opt->charlistlen != NULL)
  {
    printf("# alphalen \"%s\" %lu\n",
             gt_str_get(opt->charlistlen->charlist),
             opt->charlistlen->len);
    return;
  }
  if (gt_str_length(opt->text) > 0)
  {
    printf("# text \"%s\"\n",gt_str_get(opt->text));
    return;
  }
}

static OPrval parse_options(int *parsed_args,
                            Cmppairwiseopt *pw,
                            int argc, const char **argv, GtError *err)
{
  OptionParser *op;
  Option *optionstrings,
         *optionfiles,
         *optioncharlistlen,
         *optiontext;
  GtStrArray *charlistlen;
  OPrval oprval;

  gt_error_check(err);
  charlistlen = gt_strarray_new();
  pw->strings = gt_strarray_new();
  pw->files = gt_strarray_new();
  pw->text = gt_str_new();
  pw->charlistlen = NULL;
  op = option_parser_new("options","Apply function to pairs of strings.");
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

  optionstrings = option_new_stringarray("ss","use two strings",
                                         pw->strings);
  option_parser_add_option(op, optionstrings);

  optionfiles = option_new_filenamearray("ff","use two files",
                                         pw->files);
  option_parser_add_option(op, optionfiles);

  optioncharlistlen = option_new_stringarray("a",
                                             "use character list and length",
                                             charlistlen);
  option_parser_add_option(op, optioncharlistlen);

  optiontext = option_new_string("t","use text",pw->text, NULL);
  option_parser_add_option(op, optiontext);

  option_exclude(optionstrings, optionfiles);
  option_exclude(optionstrings, optioncharlistlen);
  option_exclude(optionstrings, optiontext);
  option_exclude(optionfiles, optioncharlistlen);
  option_exclude(optionfiles, optiontext);
  option_exclude(optioncharlistlen, optiontext);

  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  if (oprval == OPTIONPARSER_OK)
  {
    if (option_is_set(optionstrings))
    {
      if (gt_strarray_size(pw->strings) != 2UL)
      {
        gt_error_set(err, "option -ss requires two string arguments");
        oprval = OPTIONPARSER_ERROR;
      }
    } else
    {
      if (option_is_set(optionfiles))
      {
        if (gt_strarray_size(pw->files) != 2UL)
        {
          gt_error_set(err, "option -ff requires two filename arguments");
          oprval = OPTIONPARSER_ERROR;
        }
      } else
      {
        if (option_is_set(optioncharlistlen))
        {
          long readint;

          if (gt_strarray_size(charlistlen) != 2UL)
          {
            gt_error_set(err,
                         "option -a requires charlist and length argument");
            oprval = OPTIONPARSER_ERROR;
          }
          ALLOCASSIGNSPACE(pw->charlistlen,NULL,Charlistlen,1);
          pw->charlistlen->charlist =
            gt_str_ref(gt_strarray_get_str(charlistlen,
                                                                  0));
          if (sscanf(gt_strarray_get(charlistlen,1UL),"%ld",&readint) != 1 ||
              readint < 1L)
          {
            gt_error_set(err,
                         "option -a requires charlist and length argument");
            oprval = OPTIONPARSER_ERROR;
          }
          pw->charlistlen->len = (unsigned long) readint;
        } else
        {
          if (!option_is_set(optiontext))
          {
            gt_error_set(err,
                         "use exactly one of the options -ss, -ff, -a, -t");
            oprval = OPTIONPARSER_ERROR;
          }
        }
      }
    }
  }
  option_parser_delete(op);
  if (oprval == OPTIONPARSER_OK && *parsed_args != argc)
  {
    gt_error_set(err, "superfluous program parameters");
    oprval = OPTIONPARSER_ERROR;
  }
  gt_strarray_delete(charlistlen);
  return oprval;
}

static void freesimpleoption(Cmppairwiseopt *cmppairwise)
{
  gt_strarray_delete(cmppairwise->strings);
  gt_strarray_delete(cmppairwise->files);
  gt_str_delete(cmppairwise->text);
  if (cmppairwise->charlistlen != NULL)
  {
    gt_str_delete(cmppairwise->charlistlen->charlist);
    FREESPACE(cmppairwise->charlistlen);
  }
}

static unsigned long applycheckfunctiontosimpleoptions(
                                  Checkcmppairfuntype checkfunction,
                                  const Cmppairwiseopt *opt)
{
  if (gt_strarray_size(opt->strings) > 0)
  {
    bool forward = true;
    while (true)
    {
      checkfunction(forward,
                    (const Uchar *) gt_strarray_get(opt->strings,0),
                    (unsigned long) strlen(gt_strarray_get(opt->strings,0)),
                    (const Uchar *) gt_strarray_get(opt->strings,1UL),
                    (unsigned long) strlen(gt_strarray_get(opt->strings,1UL)));
      if (!forward)
      {
        break;
      }
      forward = false;
    }
    return 2UL; /* number of testcases */
  }
  if (gt_strarray_size(opt->files) > 0)
  {
    runcheckfunctionontwofiles(checkfunction,
                               gt_strarray_get(opt->files,0),
                               gt_strarray_get(opt->files,1UL));
    return 2UL;
  }
  if (opt->charlistlen != NULL)
  {
    return runcheckfunctiononalphalen(checkfunction,
                                      gt_str_get(opt->charlistlen->charlist),
                                      opt->charlistlen->len);
  }
  if (gt_str_length(opt->text) > 0)
  {
    return runcheckfunctionontext(checkfunction,gt_str_get(opt->text));
  }
  assert(false);
  return 0;
}

int gt_paircmp(int argc, const char **argv, GtError *err)
{
  int parsed_args;
  Cmppairwiseopt cmppairwise;
  OPrval oprval;

  gt_error_check(err);

  oprval = parse_options(&parsed_args,&cmppairwise,argc, argv, err);
  if (oprval == OPTIONPARSER_OK)
  {
    unsigned long testcases;

    assert(parsed_args == argc);
    showsimpleoptions(&cmppairwise);
    testcases = applycheckfunctiontosimpleoptions(checkgreedyunitedist,
                                                  &cmppairwise);
    printf("# number of testcases: %lu\n",testcases);
  }
  freesimpleoption(&cmppairwise);
  if (oprval == OPTIONPARSER_REQUESTS_EXIT)
  {
    return 0;
  }
  if (oprval == OPTIONPARSER_ERROR)
  {
    return -1;
  }
  return 0;
}
