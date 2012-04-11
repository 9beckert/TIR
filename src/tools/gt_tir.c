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

// Struktur mit allen Operationen
typedef struct {
  bool bool_option_tir;
  GtStr  *str_option_tir;
} GtTirArguments;

// Struktur initialisierung und löschen
static void* gt_tir_arguments_new(void)
{
  GtTirArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->str_option_tir = gt_str_new();
  return arguments;
}

static void gt_tir_arguments_delete(void *tool_arguments)
{
  GtTirArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->str_option_tir);
  gt_free(arguments);
}

// Abhängigkeiten
static GtOptionParser* gt_tir_option_parser_new(void *tool_arguments)
{
  GtTirArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [file]", /* XXX */	// Beschreibung der Syntax
                         "DESCRIBE YOUR TOOL IN ONE LINE HERE."); /* XXX */

	// Optionen hinzufügen
  /* -bool */
  option = gt_option_new_bool("bool", "bool option tir",	// Name, kurze Beschreibung
                           &arguments->bool_option_tir, false);	// Zeiger auf Ergebnisstelle, default Wert
  gt_option_parser_add_option(op, option);	// Einhängen in Parser

  /* -str */
  option = gt_option_new_string("str", "str option tir",
                             arguments->str_option_tir, NULL);
  gt_option_parser_add_option(op, option);

  return op;
}

// Überprüfen ob zusätzliche Bedingungen erfüllt sind oder nicht
static int gt_tir_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtTirArguments *arguments = tool_arguments;
  int had_err = 0; // sagt aus, ob Fehler aufgetreten ist, 0 kein Fehler, <0 Fehler
  gt_error_check(err);
  gt_assert(arguments);

  /* XXX: do some checking after the option have been parsed (usally this is not
     necessary and this function can be removed completely). */
  if (gt_str_length(arguments->str_option_tir))
    printf("%s\n", gt_str_get(arguments->str_option_tir));

  return had_err;
}

// Tatsächliche Methode, hier passiert die Action! :)
static int gt_tir_runner(int argc, const char **argv, int parsed_args,	// alle Optionen bis parsed_args wurden bereits bearbeitet
                              void *tool_arguments, GT_UNUSED GtError *err) // Argumentstruktur, Fehlermeldungen
{
  GtTirArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* XXX */
  if (arguments->bool_option_tir)
    printf("argc=%d, parsed_args=%d\n", argc, parsed_args);
  printf("argv[0]=%s\n", argv[0]);

  return had_err;
}

// Hier wird alles zusammengeführt, Angabe was verknüpft wird, wo null steht wird nicht angesprungen
GtTool* gt_tir(void)
{
  return gt_tool_new(gt_tir_arguments_new,
                  gt_tir_arguments_delete,
                  gt_tir_option_parser_new,
                  gt_tir_arguments_check,
                  gt_tir_runner);
}
