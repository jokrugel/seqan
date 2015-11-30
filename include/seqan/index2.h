// ==========================================================================
//                                   index2
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Johannes Krugel <krugel@in.tum.de>
// ==========================================================================
// Facade header for module index2.
// ==========================================================================

#ifndef SEQAN_INDEX2_H_
#define SEQAN_INDEX2_H_

// ===========================================================================
// Prerequisites.
// ===========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/sequence.h>

// ===========================================================================
// Basics
// ===========================================================================

#include <seqan/index2/index2_base.h>

// ===========================================================================
// q-gram indexes
// ===========================================================================

#include <seqan/index2/index_qgram2l.h>
#include <seqan/index2/find_index_qgram_ext.h>

// ===========================================================================
// Compressed indexes
// ===========================================================================

#include <seqan/index2/index_compressed_bits.h>
#include <seqan/index2/index_lz.h>
#include <seqan/index2/index_sadakane.h>
#include <seqan/index2/find_index_lz.h>
#include <seqan/index2/find_index_sadakane.h>

// ===========================================================================
// Suffix trees
// ===========================================================================

#include <seqan/index2/index_digest_base.h>
#include <seqan/index2/find_index_digest_pattern.h>
#include <seqan/index2/find_index_digest_finder.h>

#include <seqan/index2/index_sttd64_base.h>
#include <seqan/index2/index_sttd64_vptree.h>
#include <seqan/index2/index_sttd64_stree.h>

#include <seqan/index2/find_index_suffix_tree.h>

#endif  // SEQAN_INDEX2_H_
