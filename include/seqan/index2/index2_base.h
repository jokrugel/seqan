// ==========================================================================
//                                   index2
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Meta functions for the module index2.
// This module contains some additional index classes for string matching.
// ==========================================================================

#ifndef SEQAN_INDEX2_INDEX2_BASE_H_
#define SEQAN_INDEX2_INDEX2_BASE_H_

namespace seqan {
   
// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// DefaultIndexPattern
// ----------------------------------------------------------------------------
// TODO(krugel) Move to find_pattern_base.h

// We use the default pattern spezialization (but e.g. for IndexDigest we need to specialize this).
template <typename TNeedle, typename TIndex>
struct DefaultIndexPattern {
    typedef typename DefaultPattern<TNeedle>::Type          Type;
};

// ----------------------------------------------------------------------------
// IsIndex<Index>
// ----------------------------------------------------------------------------
// TODO(krugel) Move to index_base.h

template <typename THaystack>
struct IsIndex {
    typedef False                                           Type;
};

template <typename TString, typename TIndexSpec>
struct IsIndex<Index<TString, TIndexSpec> > {
    typedef True                                            Type;
};

// ----------------------------------------------------------------------------
// SupportsSuffixTreeIteration<Index>
// ----------------------------------------------------------------------------
// TODO(krugel) Move to index_base.h

template <typename TIndexSpec>
struct SupportsSuffixTreeIteration {
    typedef False Type;
};

template <typename TSpec>
struct SupportsSuffixTreeIteration<IndexEsa<TSpec> > {
    typedef True Type;
};

template <typename TSpec>
struct SupportsSuffixTreeIteration<IndexWotd<TSpec> > {
    typedef True Type;
};

template <typename TOccSpec, typename TSpec>
struct SupportsSuffixTreeIteration<FMIndex<TOccSpec, TSpec> > {
    //TODO(krugel) Supports prefix tree iteration instead
    typedef False Type;
};

// ----------------------------------------------------------------------------
// GetOccurrences<Index>
// ----------------------------------------------------------------------------

// TODO(krugel) Move to index_base.h ?
// Determines the return type of getOccurrences(Iterator<VSTree>), getOccurrences(index, shape) and related
// For some indexes this is an Infix<String<SAValue> > while for others it is String<SAValue>
template <typename TIndex>
struct GetOccurrences {
    //typedef String<typename SAValue<TIndex>::Type>          Type;
    typedef typename Infix<typename Fibre<TIndex, FibreSA>::Type const>::Type Type;
};

}  // namespace seqan

#endif  // SEQAN_INDEX2_INDEX2_BASE_H_
