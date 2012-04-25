/*
  Copyright (c) 2012 Manuela Beckert <9beckert@informatik.uni-hamburg.de>
  Copyright (c) 2012 Dorle Osterode <9osterod@informatik.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

/*
 * Finds Terminal Inverted Repeats and does... something... with them.
 * To be edited when we know more.
 */

#include "core/ma.h"
#include "core/unused_api.h"
#include "tools/gt_tir.h"
#include "tools/gt_tir_stream.h"

/* struct with all arguments */
typedef struct {
  GtStr *str_indexname;         // for reading of the enhanced suffix array
  unsigned long minseedlength;
} GtTirArguments;

/*
 * Initializes the argument struct.
 */
static void* gt_tir_arguments_new(void)
{
  GtTirArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->str_indexname = gt_str_new();
  return arguments;
}

/*
 * Deletes the argument struct.
 */
static void gt_tir_arguments_delete(void *tool_arguments)
{
  GtTirArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->str_indexname);
  gt_free(arguments);
}

/*
 * Sets dependencies between arguments
 */
static GtOptionParser* gt_tir_option_parser_new(void *tool_arguments)
{
  GtTirArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *optionindex,
           *optionseed;
  
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] -index <indexname>", // syntax description
                         "Predict Terminal Inverted Repeats (TIR)."); 

  /* -index */
  optionindex = gt_option_new_string("index",         // name
                     "specify the name of the enhanced suffix " // short description
                     "array index (mandatory)",
                     arguments->str_indexname, NULL); // pointer to result, default value 
  gt_option_is_mandatory(optionindex);
  gt_option_parser_add_option(op, optionindex);       // add to parser
  
   /* -seed */
  optionseed = gt_option_new_ulong_min("seed",
                               "specify minimum seed length for"
                               " exact repeats",
                               &arguments->minseedlength,
                               20UL,
                               2UL);
  gt_option_parser_add_option(op, optionseed); 

  return op;
}

/*
 * Check for further dependencies.
 * Currently not needed.
 */
/*static int gt_tir_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtTirArguments *arguments = tool_arguments;
  int had_err = 0; // 0 no error, <0 error
  gt_error_check(err);
  gt_assert(arguments);

   XXX: do some checking after the option have been parsed (usally this is not
     necessary and this function can be removed completely). 
  if (gt_str_length(arguments->str_option_tir))
    printf("%s\n", gt_str_get(arguments->str_option_tir));

  return had_err;
}*/

/*
 * Print arguments to console.
 */
static void gt_tir_showargsline(int argc, const char **argv)
{
  int i;
  gt_assert(argv && argc >= 1);
  printf("# args=");
  for (i=1; i<argc; i++) {
    printf("%s", argv[i]);
    if (i != argc-1) printf(" ");
  }
  printf("\n");
}

/*
 * Here is the action!!! :)
 */
static int gt_tir_runner(int argc, const char **argv, 
  GT_UNUSED int parsed_args,	      // all arguments till parsed_args were processed already
  void *tool_arguments,   // argument struct
  GtError *err) // error messages
{
  /* Here we'll add nodestream later */
  
  GtTirArguments *arguments = tool_arguments;
  GtNodeStream *tir_stream;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* Create TIR stream MUHAHA */
  tir_stream = gt_tir_stream_new(arguments->str_indexname,
                                 arguments->minseedlength,
                                 err);
  
  if (tir_stream == NULL)
  {
    return -1;
  }
  
  /* output arguments line */
  gt_tir_showargsline(argc, argv);
  
  /* pull the features through the stream and free them afterwards */
  if (!had_err)
    had_err = gt_node_stream_pull(tir_stream, err);
    
  /* free */
  gt_node_stream_delete(tir_stream);

  return had_err;
}


/* 
 * Combination of everything,
 * NULL-arguments won't be added.
 */
GtTool* gt_tir(void)
{
  return gt_tool_new(gt_tir_arguments_new,
                  gt_tir_arguments_delete,
                  gt_tir_option_parser_new,
                  NULL, //gt_tir_arguments_check,
                  gt_tir_runner);
}
