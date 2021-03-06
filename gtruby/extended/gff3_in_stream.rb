#
# Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

require 'gtdlload'
require 'gthelper'
require 'core/str_array'
require 'extended/genome_stream'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  typealias "bool", "ibool"
  extern "GtNodeStream* gt_gff3_in_stream_new_sorted(const char *, bool)"
  extern "GtStrArray* gt_gff3_in_stream_get_used_types(GtNodeStream*)"

  class GFF3InStream < GenomeStream
    def initialize(filename)
      if not File.readable?(filename)
        GT.gterror("file '#{filename}' not readable")
      end
      @genome_stream = GT.gt_gff3_in_stream_new_sorted(filename, false)
      @genome_stream.free = GT::symbol("gt_node_stream_delete", "0P")
    end

    def get_used_types
      str_array_ptr = GT.gt_gff3_in_stream_get_used_types(@genome_stream)
      used_types = GT::StrArray.new(str_array_ptr)
      used_types.to_a
    end
  end
end
