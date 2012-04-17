/*
  Copyright (c) 2012 Manuela Beckert <9beckert@informatik.uni-hamburg.de>
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

#include "core/ma.h"
#include "core/unused_api.h"
#include "tools/gt_tir.h"

// Struktur mit allen Parametern
typedef struct {
  GtStr *str_indexname;         // darüber wird enhanced suffix array eingelesen
  unsigned long minseedlength;
} GtTirArguments;

// Struktur initialisieren und löschen
static void* gt_tir_arguments_new(void)
{
  GtTirArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->str_indexname = gt_str_new();
  return arguments;
}

static void gt_tir_arguments_delete(void *tool_arguments)
{
  GtTirArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->str_indexname);
  gt_free(arguments);
}

// Abhängigkeiten
static GtOptionParser* gt_tir_option_parser_new(void *tool_arguments)
{
  GtTirArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *optionindex,
           *optionseed;
  
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] -index <indexname>", // Beschreibung der Syntax
                         "Predict Terminal Inverted Repeats (TIR)."); 

  // Optionen hinzufügen
  /* -bool */
  //option = gt_option_new_bool("bool", "bool option tir",	// Name, kurze Beschreibung
  //                         &arguments->bool_option_tir, false);	// Zeiger auf Ergebnisstelle, default Wert
  //gt_option_parser_add_option(op, option);	// Einhängen in Parser

  /* -index */
  optionindex = gt_option_new_string("index",
                             "specify the name of the enhanced suffix "
                             "array index (mandatory)",
                             arguments->str_indexname, NULL);
  gt_option_is_mandatory(optionindex);
  gt_option_parser_add_option(op, optionindex);
  
   /* -seed */
  optionseed = gt_option_new_ulong_min("seed",
                               "specify minimum seed length for"
                               " exact repeats",
                               &arguments->minseedlength,
                               30UL,
                               1UL);
  gt_option_parser_add_option(op, optionseed); 

  return op;
}

// Überprüfen ob zusätzliche Bedingungen erfüllt sind 
/*static int gt_tir_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtTirArguments *arguments = tool_arguments;
  int had_err = 0; // sagt aus, ob Fehler aufgetreten ist, 0 kein Fehler, <0 Fehler
  gt_error_check(err);
  gt_assert(arguments);

  /* XXX: do some checking after the option have been parsed (usally this is not
     necessary and this function can be removed completely). 
  if (gt_str_length(arguments->str_option_tir))
    printf("%s\n", gt_str_get(arguments->str_option_tir));

  return had_err;
}*/

/* Ausgabe der Parameter an der Konsole */
static void gt_ltrharvest_showargsline(int argc, const char **argv)
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

// Tatsächliche Methode, hier passiert die Action! :)

static int gt_tir_runner(int argc, const char **argv, 
  int parsed_args,	      // alle Optionen bis parsed_args wurden bereits bearbeitet
  void *tool_arguments,   // Argumentstruktur
  GT_UNUSED GtError *err) // Fehlermeldungen
{
  // Hier kommt später mal Nodestream
  
  GtTirArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  // TIR  stream erzeugen MUHAHA
  tir_stream = gt_tir_stream_new(argruments->str_indexname,
                                 arguments->minseedlength,
                                 err);
  
  if (ltrh_stream == NULL)
  {
    return -1;
  }
  
  /* output arguments line */
  gt_ltrharvest_showargsline(argc, argv);
  
  /* pull the features through the stream and free them afterwards */
  if (!had_err)
    had_err = gt_node_stream_pull(last_stream, err);
    
  /* free */
  gt_node_stream_delete(tir_stream);

  return had_err;
}

// Hier wird alles zusammengeführt, Angabe was verknüpft wird, wo null steht wird nicht angesprungen
GtTool* gt_tir(void)
{
  return gt_tool_new(gt_tir_arguments_new,
                  gt_tir_arguments_delete,
                  gt_tir_option_parser_new,
                  NULL, //gt_tir_arguments_check,
                  gt_tir_runner);
}
