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
// Author: Johannes Krugel <krugel@in.tum.de>
// ==========================================================================
// Finder class for the dynamic programming backtracking algorithm.
// ==========================================================================

#ifndef SEQAN_FIND_INDEX_APPROX_FIND_INDEX_BACKTRACKING_FINDER_BASE_H_
#define SEQAN_FIND_INDEX_APPROX_FIND_INDEX_BACKTRACKING_FINDER_BASE_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template<typename TScore = EditDistanceScore>
struct DPBacktracking { };

template <typename TText, typename TIndexSpec, typename TScore>
struct Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > : FindState2_ {
    typedef Index<TText, TIndexSpec>                        TIndex;
    typedef typename Position<TIndex>::Type                 TPosition;
    typedef typename Iterator<TIndex, TopDown<ParentLinks<> > >::Type TIterator;

    typedef typename GetOccurrences<TIndex>::Type           TOccurrences;
    typedef typename Iterator<TOccurrences>::Type           TOccurrencesIterator;

    TState _state;                       // The finder's state.
    Holder<TIndex> _holder;              // The index we are working on.
    TPosition _beginPosition;            // The begin position of the current match.
    TPosition _endPosition;              // The end position of the current match.

    TIterator _data_iterator;            // The suffix tree top-down-iterator
    TOccurrences _occurrences;           // The found occurrences in the current suffix tree node.
    TOccurrencesIterator _occIt;         // The current occurrence.

    Finder()
        : _state(STATE_EMPTY)
        , _holder()
        , _data_iterator()
    {
    }

    Finder(TIndex & index)
        : _state(STATE_INITIAL)
        , _holder(index)
        , _data_iterator(index)
    {
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// haystack, host
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TScore>
inline Index<TText, TIndexSpec> &
haystack(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder) {
    return value(finder._holder);
}

template <typename TText, typename TIndexSpec, typename TScore>
inline Index<TText, TIndexSpec> const &
haystack(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > const & finder) {
    return value(finder._holder);
}

template <typename TText, typename TIndexSpec, typename TScore>
inline Index<TText, TIndexSpec> &
host(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder) {
    return value(finder._holder);
}

template <typename TText, typename TIndexSpec, typename TScore>
inline Index<TText, TIndexSpec> const &
host(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > const & finder) {
    return value(finder._holder);
}

// ----------------------------------------------------------------------------
// begin, end
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TScore, typename TTag>
inline typename Iterator<Index<TText, TIndexSpec> const>::Type
begin(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > const & finder, Tag<TTag> const & spec) {
    typedef Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > TFinder;
    // for compatibility with exact finders we also allow STATE_FOUND here
    SEQAN_ASSERT(finder._state == TFinder::STATE_BEGIN_FOUND || finder._state == TFinder::STATE_FOUND);
    return begin(haystack(finder), spec) + finder._beginPosition;
}

template <typename TText, typename TIndexSpec, typename TScore, typename TTag>
inline typename Iterator<Index<TText, TIndexSpec> >::Type
begin(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder, Tag<TTag> const & spec) {
    typedef Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > TFinder;
    return begin(const_cast<TFinder const &>(finder), spec);
}

template <typename TText, typename TIndexSpec, typename TScore, typename TTag>
inline typename Iterator<Index<TText, TIndexSpec> const>::Type
end(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > const & finder, Tag<TTag> const & spec) {
    typedef Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > TFinder;
    SEQAN_ASSERT(finder._state == TFinder::STATE_BEGIN_FOUND || finder._state == TFinder::STATE_FOUND);
    return begin(haystack(finder), spec) + finder._endPosition;
}

template <typename TText, typename TIndexSpec, typename TScore, typename TTag>
inline typename Iterator<Index<TText, TIndexSpec> >::Type
end(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder, Tag<TTag> const & spec) {
    typedef Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > TFinder;
    return end(const_cast<TFinder const &>(finder), spec);
}

// ----------------------------------------------------------------------------
// beginPosition, endPosition, position
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TScore>
inline typename Position<Index<TText, TIndexSpec> const>::Type
beginPosition(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > const & finder) {
    typedef Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > TFinder;
    // for compatibility with exact finders we also allow STATE_FOUND here
    SEQAN_ASSERT(finder._state == TFinder::STATE_BEGIN_FOUND || finder._state == TFinder::STATE_FOUND);
    return finder._beginPosition;
}

template <typename TText, typename TIndexSpec, typename TScore>
inline typename Position<Index<TText, TIndexSpec> >::Type
beginPosition(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder) {
    typedef Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > TFinder;
    return beginPosition(const_cast<TFinder const &>(finder));
}

template <typename TText, typename TIndexSpec, typename TScore>
inline typename Position<Index<TText, TIndexSpec> const>::Type
endPosition(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > const & finder) {
    typedef Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > TFinder;
    SEQAN_ASSERT(finder._state == TFinder::STATE_BEGIN_FOUND || finder._state == TFinder::STATE_FOUND);
    return finder._endPosition;
}

template <typename TText, typename TIndexSpec, typename TScore>
inline typename Position<Index<TText, TIndexSpec> >::Type
endPosition(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder) {
    typedef Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > TFinder;
    return endPosition(const_cast<TFinder const &>(finder));
}

template <typename TText, typename TIndexSpec, typename TScore>
inline typename Position<Index<TText, TIndexSpec> const>::Type
position(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > const & finder) {
    typedef Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > TFinder;
    // for compatibility with IndexEsa we also allow STATE_FOUND here
    SEQAN_ASSERT(finder._state == TFinder::STATE_BEGIN_FOUND || finder._state == TFinder::STATE_FOUND);
    return endPosition(finder) - 1u;
}

template <typename TText, typename TIndexSpec, typename TScore>
inline typename Position<Index<TText, TIndexSpec> >::Type
position(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder) {
    typedef Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > TFinder;
    return position(const_cast<TFinder const &>(finder));
}

// ----------------------------------------------------------------------------
// length
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TScore>
inline typename Size<Index<TText, TIndexSpec> const>::Type
length(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > const & finder) {
    return endPosition(finder) - beginPosition(finder);
}

template <typename TText, typename TIndexSpec, typename TScore>
inline typename Size<Index<TText, TIndexSpec> >::Type
length(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder) {
    typedef Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > TFinder;
    return length(const_cast<TFinder const &>(finder));
}

// ----------------------------------------------------------------------------
// empty
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TScore>
inline bool
empty(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > const & finder) {
    return endPosition(finder) == beginPosition(finder);
}

template <typename TText, typename TIndexSpec, typename TScore>
inline bool
empty(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder) {
    typedef Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > TFinder;
    return empty(const_cast<TFinder const &>(finder));
}

// ----------------------------------------------------------------------------
// setHost, setHaystack
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TScore>
inline void
setHost(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder, typename Parameter_<Index<TText, TIndexSpec> >::Type & index) {
    typedef Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > TFinder;
    typedef Index<TText, TIndexSpec>                        TIndex;
    typedef typename Iterator<TIndex, TopDown<ParentLinks<> > >::Type TIterator;
    SEQAN_ASSERT(finder._state == TFinder::STATE_EMPTY
              || finder._state == TFinder::STATE_INITIAL);
    finder._state = TFinder::STATE_INITIAL;
    setValue(finder._holder, index);
    finder._data_iterator = TIterator(index);
}

template <typename TText, typename TIndexSpec, typename TScore>
inline void
setHaystack(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder, typename Parameter_<Index<TText, TIndexSpec> >::Type & index) {
    setHost(finder, index);
}

// ----------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TScore>
inline void
clear(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder) {
    typedef Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > TFinder;
    if (finder._state != TFinder::STATE_EMPTY) finder._state = TFinder::STATE_INITIAL;
    goBegin(finder._data_iterator);
    //clear(finder._occIt);
    // TODO(krugel) clear(finder._occurrences);  // Does not compile for infixes?
}

// ----------------------------------------------------------------------------
// goBegin
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TScore>
inline void
goBegin(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder) {
    clear(finder);
}

}  // namespace seqan
          
#endif  // SEQAN_FIND_INDEX_APPROX_FIND_INDEX_BACKTRACKING_FINDER_BASE_H_
