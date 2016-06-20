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
// Finder class for the partitioning algorithms.
// ==========================================================================

#ifndef SEQAN_FIND_INDEX_APPROX_FIND_INDEX_PARTITIONING_FINDER_H_
#define SEQAN_FIND_INDEX_APPROX_FIND_INDEX_PARTITIONING_FINDER_H_

#include <set>

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template<typename TPieceFinderSpec = Default, typename TVerifyFinderSpec = Default, typename DoPreparePatterns = False>
struct FinderPartitioning { };

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
struct Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > : FindState2_ {
    typedef Index<TText, TIndexSpec>                        TIndex;
    typedef typename Position<TIndex>::Type                 TPosition;
    typedef Finder<TIndex, TPieceFinderSpec>                TPieceFinder;
    typedef typename Infix<TText>::Type                     TTextInfix;
    typedef Finder<TTextInfix, TVerifyFinderSpec>           TVerifyFinder;
    typedef Map<Pair<TText, String<TPosition> > >           TPreparedPositions;
    // TODO(krugel) Prepare: somehow use TNeedleInfix instead

    TState _state;                       // The finder's state.
    Holder<TIndex> _holder;              // The index we are working on.
    TPosition _beginPosition;            // The begin position of the current match.
    TPosition _endPosition;              // The end position of the current match.

    TPieceFinder _pieceFinder;           // The finder object to find the pieces of the needle in the text.
    TVerifyFinder _verifyFinder;         // The finder object to verify the candidate matches of the needle in the text.
    TTextInfix _verifyRegion;            // The region in the text to be verified.
    std::set<TPosition> _returnedOcc;    // Contains the already returned occurrences (necessary to return each occurrence only once).
    
    TPreparedPositions _preparedPiecePositions;

#ifndef NDEBUG
    unsigned _candidatePositions;
#endif

    Finder()
        : _state(STATE_EMPTY)
        , _holder()
        , _pieceFinder()
        , _verifyFinder()
        , _verifyRegion()
        , _returnedOcc()
#ifndef NDEBUG
        , _candidatePositions(0)
#endif
    {
    }

    Finder(TIndex & index)
        : _state(STATE_INITIAL)
        , _holder(index)
        , _pieceFinder(index)
        , _verifyFinder()
        , _verifyRegion()
        , _returnedOcc()
#ifndef NDEBUG
        , _candidatePositions(0)
#endif
    {
    }
};

// ----------------------------------------------------------------------------
// Provide default type parameters
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
struct Finder<Index<TText, TIndexSpec>, FinderPartitioning<Default, TVerifyFinderSpec, DoPreparePatterns> >
        : Finder<Index<TText, TIndexSpec>, FinderPartitioning<typename DefaultFinder<Index<TText, TIndexSpec> >::Type, TVerifyFinderSpec, DoPreparePatterns> > {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<typename DefaultFinder<Index<TText, TIndexSpec> >::Type, TVerifyFinderSpec, DoPreparePatterns> > TBase;
    typedef Index<TText, TIndexSpec>                        TIndex;
    Finder() : TBase() { }
    Finder(TIndex & index) : TBase(index) { }
};

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename DoPreparePatterns>
struct Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, Default, DoPreparePatterns> >
        : Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, typename DefaultFinder<typename Infix<TText>::Type>::Type, DoPreparePatterns> > {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, typename DefaultFinder<typename Infix<TText>::Type>::Type, DoPreparePatterns> > TBase;
    typedef Index<TText, TIndexSpec>                        TIndex;
    Finder() : TBase() { }
    Finder(TIndex & index) : TBase(index) { }
};

template <typename TText, typename TIndexSpec, typename DoPreparePatterns>
struct Finder<Index<TText, TIndexSpec>, FinderPartitioning<Default, Default, DoPreparePatterns> >
        : Finder<Index<TText, TIndexSpec>, FinderPartitioning<typename DefaultFinder<Index<TText, TIndexSpec> >::Type, typename DefaultFinder<typename Infix<TText>::Type>::Type, DoPreparePatterns> > {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<typename DefaultFinder<Index<TText, TIndexSpec> >::Type, typename DefaultFinder<typename Infix<TText>::Type>::Type, DoPreparePatterns> > TBase;
    typedef Index<TText, TIndexSpec>                        TIndex;
    Finder() : TBase() { }
    Finder(TIndex & index) : TBase(index) { }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#ifndef NDEBUG
template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline void
printDebugInfo(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder) {
    std::cerr << "candidate positions = " << finder._candidatePositions << std::endl;
}
#endif

// ----------------------------------------------------------------------------
// preparePatterns
// ----------------------------------------------------------------------------

template <typename TFinder, typename TPattern, typename TNeedle, typename TScoreValue>
inline void
preparePatterns(TFinder & finder, TPattern & pattern, StringSet<TNeedle> & needles, TScoreValue tolerance) {
    // Do nothing
    ignoreUnusedVariableWarning(finder);
    ignoreUnusedVariableWarning(pattern);
    ignoreUnusedVariableWarning(needles);
    ignoreUnusedVariableWarning(tolerance);
}

// This specialization is for DoPreparePatterns == True and at the moment only used for IndexDigest to preprocess a set of patterns
template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename TPattern, typename TNeedle, typename TScoreValue>
inline void
preparePatterns(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, True> > & finder, TPattern & pattern, StringSet<TNeedle> & needles, TScoreValue tolerance) {
    
    typedef True                                            DoPreparePatterns;
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    typedef typename TFinder::TPreparedPositions            TPreparedPositions;
    typedef typename Cargo<TPreparedPositions>::Type        TPreparedPositionsValue;
    
    typedef StringSet<TNeedle>                              TNeedles;
    typedef typename Iterator<TNeedles>::Type               TNeedlesIter;
    typedef typename Infix<TNeedle>::Type                   TNeedleInfix;
    typedef String<TNeedleInfix>                            TPieces;
    typedef typename Iterator<TPieces>::Type                TPiecesIter;
    
    for (TNeedlesIter iter = begin(needles); ! atEnd(iter, needles); ++iter) {
        TPieces tmpPieces;
        reserve(tmpPieces, tolerance + 1u, Exact());
        _split(value(iter), tolerance + 1u, tmpPieces);
        
        for (TPiecesIter pieceIter = begin(tmpPieces); ! atEnd(pieceIter, tmpPieces); ++pieceIter) {
            // TODO(krugel) Use find and insert instead
            finder._preparedPiecePositions[*pieceIter] = TPreparedPositionsValue();
        }
    }
    
    preparePiecePositions(finder._pieceFinder, pattern._piecePattern, finder._preparedPiecePositions);
}

// ----------------------------------------------------------------------------
// preparePiecePositions
// ----------------------------------------------------------------------------

template <typename TPieceFinder, typename TPattern, typename TPreparedPositions>
inline void
preparePiecePositions(TPieceFinder & pieceFinder, TPattern & pattern, TPreparedPositions & preparedPositions) {
    // Do nothing
    ignoreUnusedVariableWarning(pieceFinder);
    ignoreUnusedVariableWarning(pattern);
    ignoreUnusedVariableWarning(preparedPositions);
}

// ----------------------------------------------------------------------------
// setPreparedPiecePositions
// ----------------------------------------------------------------------------

template <typename TFinder, typename TPreparedPositions>
inline void
setPreparedPiecePositions(TFinder & finder, TPreparedPositions & preparedPositions) {
    // Do nothing
    ignoreUnusedVariableWarning(finder);
    ignoreUnusedVariableWarning(preparedPositions);
}


// ----------------------------------------------------------------------------
// haystack, host
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline Index<TText, TIndexSpec> &
haystack(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder) {
    return value(finder._holder);
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline Index<TText, TIndexSpec> const &
haystack(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > const & finder) {
    return value(finder._holder);
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline Index<TText, TIndexSpec> &
host(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder) {
    return value(finder._holder);
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline Index<TText, TIndexSpec> const &
host(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > const & finder) {
    return value(finder._holder);
}

// ----------------------------------------------------------------------------
// begin, end
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns, typename TTag>
inline typename Iterator<Index<TText, TIndexSpec> const>::Type
begin(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > const & finder, Tag<TTag> const & spec) {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    SEQAN_ASSERT_EQ(finder._state, TFinder::STATE_BEGIN_FOUND);
    return begin(haystack(finder), spec) + beginPosition(finder);
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns, typename TTag>
inline typename Iterator<Index<TText, TIndexSpec> const>::Type
begin(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder, Tag<TTag> const & spec) {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    return begin(const_cast<TFinder const &>(finder), spec);
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns, typename TTag>
inline typename Iterator<Index<TText, TIndexSpec> const>::Type
end(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > const & finder, Tag<TTag> const & spec) {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    SEQAN_ASSERT(finder._state == TFinder::STATE_BEGIN_FOUND || finder._state == TFinder::STATE_FOUND);
    return begin(haystack(finder), spec) + endPosition(finder);
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns, typename TTag>
inline typename Iterator<Index<TText, TIndexSpec> const>::Type
end(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder,
        Tag<TTag> const & spec) {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    return end(const_cast<TFinder const &>(finder), spec);
}

// ----------------------------------------------------------------------------
// beginPosition, endPosition, position
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline typename Position<Index<TText, TIndexSpec> const>::Type
beginPosition(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > const & finder) {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    SEQAN_ASSERT(finder._state == TFinder::STATE_BEGIN_FOUND || finder._state == TFinder::STATE_FOUND);
    //return beginPosition(finder._verifyRegion) - beginPosition(haystack(finder)) + beginPosition(finder._verifyFinder);
    return finder._beginPosition;
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline typename Position<Index<TText, TIndexSpec> >::Type
beginPosition(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder) {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    return beginPosition(const_cast<TFinder const &>(finder));
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline typename Position<Index<TText, TIndexSpec> const>::Type
endPosition(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > const & finder) {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    SEQAN_ASSERT(finder._state == TFinder::STATE_BEGIN_FOUND || finder._state == TFinder::STATE_FOUND);
    //return beginPosition(finder._verifyRegion) - beginPosition(haystack(finder)) + endPosition(finder._verifyFinder);
    return finder._endPosition;
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline typename Position<Index<TText, TIndexSpec> >::Type
endPosition(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder) {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    return endPosition(const_cast<TFinder const &>(finder));
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline typename Position<Index<TText, TIndexSpec> const>::Type
position(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > const & finder) {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    SEQAN_ASSERT(finder._state == TFinder::STATE_BEGIN_FOUND || finder._state == TFinder::STATE_FOUND);
    return endPosition(finder) - 1u;
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline typename Position<Index<TText, TIndexSpec> >::Type
position(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder) {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    return position(const_cast<TFinder const &>(finder));
}

// ----------------------------------------------------------------------------
// length
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline typename Size<Index<TText, TIndexSpec> const>::Type
length(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > const & finder) {
    return endPosition(finder) - beginPosition(finder);
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline typename Size<Index<TText, TIndexSpec> >::Type
length(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder) {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    return length(const_cast<TFinder const &>(finder));
}

// ----------------------------------------------------------------------------
// empty
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline bool
empty(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > const & finder) {
    return endPosition(finder) == beginPosition(finder);
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline bool
empty(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder) {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    return empty(const_cast<TFinder const &>(finder));
}

// ----------------------------------------------------------------------------
// setHost, setHaystack
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline void
setHost(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder, typename Parameter_<Index<TText, TIndexSpec> >::Type index) {
    //typedef Index<TText, TIndexSpec>                        TIndex;
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    SEQAN_ASSERT(finder._state == TFinder::STATE_EMPTY
              || finder._state == TFinder::STATE_INITIAL);
    finder._state = TFinder::STATE_INITIAL;
    //finder._holder = Holder<TIndex>(index);
    setValue(finder._holder, index);
    setHost(finder._pieceFinder, haystack(finder));
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline void
setHaystack(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder, typename Parameter_<Index<TText, TIndexSpec> >::Type & index) {
    setHost(finder, index);
}

// ----------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline void
clear(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder) {
    typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    if (finder._state != TFinder::STATE_EMPTY) finder._state = TFinder::STATE_INITIAL;
    clear(finder._pieceFinder);
    clear(finder._verifyFinder);
    clear(finder._returnedOcc);
    // TODO(krugel) Do not clear: clear(finder._preparedPiecePositions);
}

// ----------------------------------------------------------------------------
// goBegin
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns>
inline void
goBegin(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder) {
    clear(finder);
}

}  // namespace seqan
          
#endif  // SEQAN_FIND_INDEX_APPROX_FIND_INDEX_PARTITIONING_FINDER_H_
