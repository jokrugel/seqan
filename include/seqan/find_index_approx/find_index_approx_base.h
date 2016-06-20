// ==========================================================================
//                             find_index_approx
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
// Author: Henrik Poppe <poppeh@in.tum.de>
// Author: Johannes Krugel <krugel@in.tum.de>
// ==========================================================================
// Base for module find_index_approx.
// Algorithms to find the occurences of a search pattern in a text allowing
// some errors and using an index structure. The algorithms implemented in
// this module are:
// - DPBacktracking: Classical dynamic programming algorithm performing
//   backtracking in a suffix tree. [NBYST2001]
// - IntoExactSearch: A partitioning algorithm using the Pigeonhole principle.
// - Intermediate: A mixture of the above 2 algorithms. [NBY2000]
// - Hierarchical: Partitioning with hierarchical verification. [Myers1994,
//   HN2003]
// 
// [NBYST2001] G. Navarro and R.A. Baeza-Yates and E. Sutinen and J. Tarhio
//    "Indexing Methods for Approximate String Matching", IEEE, 2001,
//    http://sites.computer.org/debull/A01DEC-CD.pdf
// [NBY2000] G. Navarro, R.A. Baeza-Yates "A Hybrid Indexing Method for
//    Approximate String Matching" Journal of Discrete Algorithms, 2000
//    http://www.dcc.uchile.cl/~gnavarro/ps/jda00.2.pdf
// [Myers1994] E.W. Myers. "A sublinear algorithm for approximate keyword
//    searching", Algorithmica 12.4-5, 1994
//    http://dx.doi.org/10.1007/BF01185432
// [HN2003] H. Hyyrö and G. Navarro "A Practical Index for Genome Searching"
//    SPIRE 2003,
//    http://dx.doi.org/10.1007/b14038.
// 
// ==========================================================================

#ifndef SEQAN_FIND_INDEX_APPROX_FIND_INDEX_APPROX_BASE_H_
#define SEQAN_FIND_INDEX_APPROX_FIND_INDEX_APPROX_BASE_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct FindState2_ {
    enum TState {
        STATE_EMPTY,           // Finder/pattern is empty.
        STATE_INITIAL,         // Finer/pattern has just been initialized.
        STATE_FOUND,           // Found the end position of a hit.
        STATE_NOTFOUND,        // No hit found, no more hits possible.
        STATE_BEGIN_FOUND,     // Found begin position.
        STATE_BEGIN_NOTFOUND,  // Found end but not begin, should not happen.
        STATE_NO_HIT           // Set manually to non-hit.
    };
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// SupportsApproximateSearch<Pattern>
// ----------------------------------------------------------------------------
// TODO(krugel) Move to find_pattern_base.h

template <typename TPatternSpec>
struct SupportsApproximateSearch {
    typedef False Type;
};

template <>
struct SupportsApproximateSearch<AbndmAlgo> {
    typedef True Type;
};

template <typename TScore, typename TSpec, typename TFindBeginPatternSpec>
struct SupportsApproximateSearch<DPSearch<TScore, TSpec, TFindBeginPatternSpec> > {
    typedef True Type;
};

template <typename TSpec, typename TFindBeginPatternSpec>
struct SupportsApproximateSearch<Myers<TSpec, TFindBeginPatternSpec> > {
    typedef True Type;
};

template <typename TSpec>
struct SupportsApproximateSearch<Swift<TSpec> > {
    typedef True Type;
};

// ----------------------------------------------------------------------------
// MatchesAreSorted<Finder>
// ----------------------------------------------------------------------------
// TODO(krugel) Move to find_base.h

template<typename TFinder>
struct MatchesAreSorted {
    typedef True Type;
};

template<typename TText, typename TIndexSpec, typename TSpec>
struct MatchesAreSorted<Finder<Index<TText, TIndexSpec>, TSpec> > {
    typedef False Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// minSupportedPatternLength
// ----------------------------------------------------------------------------

// TODO(krugel) Move to find_index.h?
// By default we assume the index does not have a minimal possible pattern length.
// Real definitions are in find_index_qgram_ext.h

template <typename TText, typename TIndexSpec, typename TFinderSpec>
inline typename Size<Finder<Index<TText, TIndexSpec>, TFinderSpec> >::Type
minSupportedPatternLength(Finder<Index<TText, TIndexSpec>, TFinderSpec> & finder) {
    typedef Index<TText, TIndexSpec>                        TIndex;
    typedef typename DefaultIndexPattern<TText, TIndex>::Type TPattern;
    TPattern pattern;
    return minSupportedPatternLength(finder, pattern);
}

template <typename THaystack, typename TFinderSpec, typename TPattern>
inline typename Size<Finder<THaystack, TFinderSpec> >::Type
minSupportedPatternLength(Finder<THaystack, TFinderSpec> & /*finder*/, TPattern & /*pattern*/) {
    return 0u;
}

// ----------------------------------------------------------------------------
// configureSave
// ----------------------------------------------------------------------------

// TODO(krugel) Move to index_base.h?
// For some indexes it is necessary to define the disk location already before
// the construction (if memory-mapped files are used). For compatibility we
// define the corresponding function here for all other indexes.

template <typename TIndex>
inline void
configureSave(TIndex & index, const char * filename) {
    ignoreUnusedVariableWarning(index);
    ignoreUnusedVariableWarning(filename);
    // do nothing
}

}  // namespace seqan

#endif  // SEQAN_FIND_INDEX_APPROX_FIND_INDEX_APPROX_BASE_H_
