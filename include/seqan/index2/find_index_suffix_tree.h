// ==========================================================================
//                                   index2
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Generic finder for any index supporting suffix tree iteration to go down
// (e.g. IndexEsa, IndexWotd, and IndexSttd64).
// ==========================================================================

#ifndef SEQAN_INDEX_FIND_INDEX_SUFFIX_TREE_H_
#define SEQAN_INDEX_FIND_INDEX_SUFFIX_TREE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template<typename TSpec = void>
struct FinderSuffixTreeGoDown {};

template <typename TText, typename TIndexSpec, typename TFinderSpec>
class Finder<Index<TText, TIndexSpec>, FinderSuffixTreeGoDown<TFinderSpec> > {
protected:
    typedef Index<TText, TIndexSpec>                        TIndex;
    typedef typename Position<TText>::Type                  TPosition;
    typedef typename Size<TIndex>::Type						TSize;

    typedef typename GetOccurrences<TIndex>::Type           TOccurrences;
    typedef typename Iterator<TOccurrences>::Type           TOccurrencesIterator;
    // TODO(krugel) Use the OccurrencesIterator class here as soon as adapted for IndexSttd64

public:
    Holder<TIndex>  index;
    TSize           data_length;
    
    TOccurrences occurrences;
    TOccurrencesIterator occurrencesIter;

    Finder() {
        SEQAN_CHECKPOINT;
        clear(*this);
    }

    Finder(TIndex &_index): index(_index) {
        SEQAN_CHECKPOINT;
        clear(*this);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TText, typename TIndexSpec>
struct DefaultFinder<Index<TText, IndexSttd64<TIndexSpec> > > {
    typedef FinderSuffixTreeGoDown<>                        Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TFinderSpec>
inline void
clear(Finder<Index<TText, TIndexSpec>, FinderSuffixTreeGoDown<TFinderSpec> > & me) {
    SEQAN_CHECKPOINT;
    me.data_length = 0;
    clear(me.occurrences);
}

// ----------------------------------------------------------------------------
// empty
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TFinderSpec>
inline bool
empty(Finder<Index<TText, TIndexSpec>, FinderSuffixTreeGoDown<TFinderSpec> > & me) {
    SEQAN_CHECKPOINT;
    return empty(me.occurrences);
}

// ----------------------------------------------------------------------------
// find
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TFinderSpec, typename TPattern>
inline bool find(Finder<Index<TText, TIndexSpec>, FinderSuffixTreeGoDown<TFinderSpec> > &finder, TPattern const &pattern) {
    typedef Index<TText, TIndexSpec>                        TIndex;
    typedef typename Iterator<TIndex, TopDown<> >::Type     TIterator;
    
    // TODO(krugel) Add extra case for finding only the first occurrence? Does not necessarily help for approximate matching
    
    // TODO(krugel) IndexSttd64: Better retrieve matches one after the other?
    if (empty(finder)) {
        TIterator it(haystack(finder));
        if (!goDown(it, needle(pattern))) return false;
        finder.occurrences = getOccurrences(it);
        if (empty(finder.occurrences)) return false;
        finder.occurrencesIter = begin(finder.occurrences);
        _setFinderLength(finder, length(needle(pattern)));
    } else {
        ++finder.occurrencesIter;
    }
    if (atEnd(finder.occurrencesIter, finder.occurrences)) return false;
    return true;
}

// ----------------------------------------------------------------------------
// beginPosition
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TFinderSpec>
inline typename Position<Finder<Index<TText, TIndexSpec>, FinderSuffixTreeGoDown<TFinderSpec> > >::Type
beginPosition(Finder<Index<TText, TIndexSpec>, FinderSuffixTreeGoDown<TFinderSpec> > & me) {
    SEQAN_CHECKPOINT;
    return value(me.occurrencesIter);
}

// ----------------------------------------------------------------------------
// setHost
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TFinderSpec>
void setHost(Finder<Index<TText, TIndexSpec>, FinderSuffixTreeGoDown<TFinderSpec> > & finder, typename Parameter_<Index<TText, TIndexSpec> >::Type & index) {
    //typedef Index<TText, TIndexSpec>  TIndex;

    SEQAN_CHECKPOINT;
    clear(finder);
    setValue(finder.index, index);
}

}  // namespace seqan

#endif  // SEQAN_INDEX_FIND_INDEX_SUFFIX_TREE_H_
