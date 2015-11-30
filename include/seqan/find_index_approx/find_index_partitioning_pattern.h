// ==========================================================================
//                             find_index_approx
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
// Pattern class for the partitioning algorithms.
// ==========================================================================

#ifndef SEQAN_FIND_INDEX_APPROX_FIND_INDEX_PARTITIONING_PATTERN_H_
#define SEQAN_FIND_INDEX_APPROX_FIND_INDEX_PARTITIONING_PATTERN_H_

#include <cmath> // std::log
#include <algorithm> // std::min, std::max

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct _IntoExactSearch;
typedef Tag<_IntoExactSearch> IntoExactSearch;

struct _Intermediate;
typedef Tag<_Intermediate> Intermediate;

struct _PartitioningHierarchical;
typedef Tag<_PartitioningHierarchical> PartitioningHierarchical;

struct DefaultVerifyPatternSpec {
    typedef Myers<FindInfix> Type;
    //typedef DPSearch<EditDistanceScore> Type;
};

template<typename TSpec = IntoExactSearch, typename TScore = EditDistanceScore,
    typename TPiecePatternSpec = Default, typename TVerifyPatternSpec = Default>
struct Partitioning { };

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
struct Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > : public FindState2_ {
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;
    typedef typename Infix<TNeedle>::Type                   TNeedleInfix;
    typedef String<TNeedleInfix>                            TPieces;
    typedef typename Size<TPieces>::Type                    TPiecesSize;
    typedef typename Value<TScore>::Type                    TScoreValue;
    typedef typename Position<TPieces>::Type                TPositionStringSet;
    typedef Pattern<TNeedleInfix, TPiecePatternSpec>        TPiecePattern;
    typedef Pattern<TNeedle, TVerifyPatternSpec>            TVerifyPattern;

    typename TPattern::TState _state;    // The pattern's state.
    Holder<TNeedle> _host;               // The needle we work on.
    TScoreValue _scoreLimit;             // The minimal score of a match.
    TScoreValue _pieceScoreLimit;        // The minimal score of a match for one piece.
    TScore _scoringScheme;               // The scoring scheme to use.
    TScoreValue _currentScore;           // The current score of a match.

    TPiecePattern _piecePattern;         // The pattern object to search for a substring of the needle.
    TVerifyPattern _verifyPattern;       // The pattern object to verify a candidate match in the text.
    TPieces _pieces;                     // Contains substrings of the needle.
    TPiecesSize _numberOfPieces;         // The number of pieces to split the needle.
    TPositionStringSet _pieceNo;         // The current position in pieces.

    Pattern()
        : _state(TPattern::STATE_EMPTY)
    {
        SEQAN_CHECKPOINT;
    }

    // TODO(krugel) Explicit definition should normally not be necessary?
    template <typename TPiecePatternSpec2, typename TVerifyPatternSpec2>
    Pattern(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec2, TVerifyPatternSpec2> > & rhs)
        : _state(rhs._state)
        , _host(rhs._host)
        , _scoreLimit(rhs._scoreLimit)
        , _scoringScheme(rhs._scoringScheme)
        , _currentScore(rhs._currentScore)
        //, _piecePattern(rhs._piecePattern)
        //, _verifyPattern(rhs._verifyPattern)
        , _pieces(rhs._pieces)
        , _numberOfPieces(rhs._numberOfPieces)
        , _pieceNo(rhs._pieceNo)
    {
        _piecePattern = rhs._piecePattern;
        _verifyPattern = rhs._verifyPattern;
        SEQAN_CHECKPOINT;
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle)
        : _state(TPattern::STATE_INITIAL)
        , _host(needle)
        , _scoreLimit(0)
        , _scoringScheme()
        , _piecePattern()
        , _verifyPattern()
        , _pieces()
        , _numberOfPieces(0u)
        , _pieceNo()
    {
        SEQAN_CHECKPOINT;
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, TScoreValue scoreLimit)
        : _state(TPattern::STATE_INITIAL)
        , _host(needle)
        , _scoreLimit(scoreLimit)
        , _scoringScheme()
        , _piecePattern()
        , _verifyPattern()
        , _pieces()
        , _numberOfPieces(0u)
        , _pieceNo()
    {
        SEQAN_CHECKPOINT;
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, TScoreValue scoreLimit, TScore scoringScheme)
        : _state(TPattern::STATE_INITIAL)
        , _host(needle)
        , _scoreLimit(scoreLimit)
        , _scoringScheme(scoringScheme)
        , _piecePattern()
        , _verifyPattern()
        , _pieces()
        , _numberOfPieces(0u)
        , _pieceNo()
    {
        SEQAN_CHECKPOINT;
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, TScoreValue scoreLimit, TScore scoringScheme, TPiecesSize numberOfPieces)
        : _state(TPattern::STATE_INITIAL)
        , _host(needle)
        , _scoreLimit(scoreLimit)
        , _scoringScheme(scoringScheme)
        , _piecePattern()
        , _verifyPattern()
        , _pieces()
        , _numberOfPieces(numberOfPieces)
        , _pieceNo(0u)
    {
        SEQAN_CHECKPOINT;
    }
};

// ----------------------------------------------------------------------------
// Provide default type parameters
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec, typename TScore, typename TVerifyPatternSpec>
struct Pattern<TNeedle, Partitioning<TSpec, TScore, Default, TVerifyPatternSpec> >
        : public Pattern<TNeedle, Partitioning<TSpec, TScore, typename DefaultPattern<typename Infix<TNeedle>::Type>::Type, TVerifyPatternSpec> > {
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, Default, TVerifyPatternSpec> > TPattern;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, typename DefaultPattern<typename Infix<TNeedle>::Type>::Type, TVerifyPatternSpec> > TBase;
    typedef typename Value<TScore>::Type                    TScoreValue;
    typedef typename Infix<TNeedle>::Type                   TNeedleInfix;
    typedef String<TNeedleInfix>                            TPieces;
    typedef typename Size<TPieces>::Type                    TPiecesSize;

    Pattern() : TBase() { }

    //template <typename TPiecePatternSpec2, typename TVerifyPatternSpec2>
    //Pattern(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec2, TVerifyPatternSpec2> > const & rhs) : TBase(rhs) {}

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle) : TBase(needle) { }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, TScoreValue scoreLimit) : TBase(needle, scoreLimit) { }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, TScoreValue scoreLimit, TScore scoringScheme) : TBase(needle, scoreLimit, scoringScheme) { }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, TScoreValue scoreLimit, TScore scoringScheme, TPiecesSize numberOfPieces) : TBase(needle, scoreLimit, scoringScheme, numberOfPieces) { }
};

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec>
struct Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, Default> >
        : public Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, typename DefaultVerifyPatternSpec::Type> > {
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, Default> > TPattern;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, typename DefaultVerifyPatternSpec::Type> > TBase;
    typedef typename Value<TScore>::Type                    TScoreValue;
    typedef typename Infix<TNeedle>::Type                   TNeedleInfix;
    typedef String<TNeedleInfix>                            TPieces;
    typedef typename Size<TPieces>::Type                    TPiecesSize;

    Pattern() : TBase() { }

    //template <typename TNeedle2, typename TPiecePatternSpec2, typename TVerifyPatternSpec2>
    //Pattern(Pattern<TNeedle2, Partitioning<TSpec, TScore, TPiecePatternSpec2, TVerifyPatternSpec2> > const & rhs) : TBase(rhs) {}

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle) : TBase(needle) { }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, TScoreValue scoreLimit) : TBase(needle, scoreLimit) { }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, TScoreValue scoreLimit, TScore scoringScheme) : TBase(needle, scoreLimit, scoringScheme) { }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, TScoreValue scoreLimit, TScore scoringScheme, TPiecesSize numberOfPieces) : TBase(needle, scoreLimit, scoringScheme, numberOfPieces) { }
};

template <typename TNeedle, typename TSpec, typename TScore>
struct Pattern<TNeedle, Partitioning<TSpec, TScore, Default, Default> >
        : public Pattern<TNeedle, Partitioning<TSpec, TScore, typename DefaultPattern<typename Infix<TNeedle>::Type>::Type, typename DefaultVerifyPatternSpec::Type> > {
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, Default, Default> > TPattern;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, typename DefaultPattern<typename Infix<TNeedle>::Type>::Type, typename DefaultVerifyPatternSpec::Type> > TBase;
    typedef typename Value<TScore>::Type                    TScoreValue;
    typedef typename Infix<TNeedle>::Type                   TNeedleInfix;
    typedef String<TNeedleInfix>                            TPieces;
    typedef typename Size<TPieces>::Type                    TPiecesSize;

    Pattern() : TBase() { }

    //template <typename TNeedle2, typename TPiecePatternSpec2, typename TVerifyPatternSpec2>
    //Pattern(Pattern<TNeedle2, Partitioning<TSpec, TScore, TPiecePatternSpec2, TVerifyPatternSpec2> > const & rhs) : TBase(rhs) {}

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle) : TBase(needle) { }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, TScoreValue scoreLimit) : TBase(needle, scoreLimit) { }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, TScoreValue scoreLimit, TScore scoringScheme) : TBase(needle, scoreLimit, scoringScheme) { }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, TScoreValue scoreLimit, TScore scoringScheme, TPiecesSize numberOfPieces) : TBase(needle, scoreLimit, scoringScheme, numberOfPieces) { }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
struct SupportsApproximateSearch<Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > {
    typedef True Type;
};

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
struct ScoringScheme<Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > > {
    typedef TScore Type;
};

// ============================================================================
// Functions (private)
// ============================================================================

// Partitioning into exact search: split the pattern in pieces, so that one piece has to occur in the text without errors
template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns, typename TNeedle, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline void
_setDefaultNumberOfPieces(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > const & /* finder not used here */,
        Pattern<TNeedle, Partitioning<IntoExactSearch, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    typedef typename Value<TScore>::Type                    TScoreValue;
    SEQAN_ASSERT_LEQ_MSG(scoreLimit(pattern), 0, "Parameter scoreLimit has to be a non-positive value, since the scoring scheme is EditDistanceScore.");

    const TScoreValue limit = -scoreLimit(pattern);
    pattern._numberOfPieces = limit + 1u;
}

// Intermediate partitioning: Heuristic from Navarro, Baeza-Yates "A Hybrid Indexing Method for Approximate String Matching"
template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns, typename TNeedle, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline void
_setDefaultNumberOfPieces(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > const & finder,
        Pattern<TNeedle, Partitioning<Intermediate, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    typedef typename Value<TScore>::Type                    TScoreValue;
    typedef typename Value<TText>::Type                     TAlphabet;
    typedef typename ValueSize<TAlphabet>::Type             TAlphabetSize;
    typedef typename Size<TText>::Type                      TSize;
    SEQAN_ASSERT_LEQ_MSG(scoreLimit(pattern), 0, "Parameter scoreLimit has to be a non-positive value, since the scoring scheme is EditDistanceScore.");

    const TAlphabetSize sigma = ValueSize<TAlphabet>::VALUE;
    const TScoreValue limit = - scoreLimit(pattern);
    double j = (length(needle(pattern)) + limit) / (std::log((double) length(haystack(finder))) / std::log((double) sigma));
    j = std::max(j, 1.0);
    j = std::min(j, limit + 1.0);
    pattern._numberOfPieces = static_cast<TSize>(std::floor(j + 0.5));
}

//f(x) : maximal tolerance on level x when using j pieces    (level 0 = innermost)
//f(x+1) = f(x)*j+j-1 = (f(x)+1)*j-1
//2, 2, 2, 2     3, 3, 3     3, 2, 2     2, 3, 3
//1, 3, 7, 15    2, 8, 26    2, 5, 11    1, 5, 17

// PartitioningHierarchical: by default use 2 pieces
template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns, typename TNeedle, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline void
_setDefaultNumberOfPieces(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > const & /*finder*/,
        Pattern<TNeedle, Partitioning<PartitioningHierarchical, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_LEQ_MSG(scoreLimit(pattern), 0, "Parameter scoreLimit has to be a non-positive value, since the scoring scheme is EditDistanceScore.");
    pattern._numberOfPieces = std::min(- scoreLimit(pattern) + 1, 2);
}

// IntoExactSearch
template <typename TNeedle, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline void
_setPieceScoreLimit(Pattern<TNeedle, Partitioning<IntoExactSearch, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    pattern._pieceScoreLimit = 0;
    return;
}

// Intermediate, PartitioningHierarchical
template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline void
_setPieceScoreLimit(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    typedef typename Value<TScore>::Type                    TScoreValue;
    SEQAN_ASSERT_GT(pattern._numberOfPieces, 0u);
    pattern._pieceScoreLimit = scoreLimit(pattern) / (TScoreValue) pattern._numberOfPieces; // round downwards
    setScoreLimit(pattern._piecePattern, pattern._pieceScoreLimit);
    return;
}

// split the pattern in numberOfPieces pieces, all (nearly) of the same size
template <typename TNeedle, typename TPosition, typename TNeedles>
inline void
_split(TNeedle & needle, TPosition numberOfPieces, TNeedles & set) {
    typedef typename Size<TNeedle>::Type                    TSize;

    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GT(numberOfPieces, 0u);
    TPosition len  = length(needle) / numberOfPieces;
    TPosition rest = length(needle) % numberOfPieces;
    TPosition beg = 0u;
    // The first pieces are one character longer
    for (TPosition i = 0u; i < numberOfPieces; ++i) {
        TSize pieceLength = (i < rest) ? (len + 1u) : len;
        appendValue(set, infixWithLength(needle, beg, pieceLength));
        //std::cerr << beg << " ~ " << infixWithLength(needle, beg, pieceLength) << " * " << begin(infixWithLength(needle, beg, pieceLength)) << " * " << beginPosition(infixWithLength(needle, beg, pieceLength)) << " * " << beginPosition(set[i]) << "|";
        //std::cerr << pieceLength << ", " << std::endl;
        beg += pieceLength;
    }
    //std::cerr << std::endl;
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// minSupportedPatternLength
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns, typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline typename Size<Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > >::Type
minSupportedPatternLength(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder,
    Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern) {
    if (pattern._numberOfPieces == 0u) _setDefaultNumberOfPieces(finder, pattern);
    _setPieceScoreLimit(pattern);
    return minSupportedPatternLength(finder._pieceFinder, pattern._piecePattern) * pattern._numberOfPieces;
}

// ----------------------------------------------------------------------------
// findBegin, find
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns, typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec, typename TScoreValue>
inline bool
findBegin(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder,
        Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern, TScoreValue scoreLimit2) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;

    SEQAN_ASSERT(pattern._state == TPattern::STATE_FOUND
              || pattern._state == TPattern::STATE_BEGIN_FOUND);

    setScoreLimit(pattern._verifyPattern, scoreLimit2);
    bool found = findBegin(finder, pattern);
    setScoreLimit(pattern._verifyPattern, pattern._scoreLimit);
    return found;
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns, typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline bool
findBegin(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder,
        Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;
    typedef typename Position<TText>::Type                  TPosition;

    SEQAN_ASSERT(pattern._state == TPattern::STATE_FOUND
              || pattern._state == TPattern::STATE_BEGIN_FOUND);
    bool found;
    TPosition pos;
    if (pattern._scoreLimit == 0u) {
        found = true;
        pos = 0u;
    } else {
        found = findBegin(finder._verifyFinder, pattern._verifyPattern);
        if (found) pos = beginPosition(finder._verifyFinder);
    }
    if (found) {
        finder._state  = TPattern::STATE_BEGIN_FOUND;
        pattern._state = TPattern::STATE_BEGIN_FOUND;
        finder._beginPosition = beginPosition(finder._verifyRegion) - beginPosition(haystack(finder)) + pos;
    } else {
        finder._state  = TPattern::STATE_BEGIN_NOTFOUND;
        pattern._state = TPattern::STATE_BEGIN_NOTFOUND;
    }
    return found;
}

// TODO(krugel) Move to find_index.h after refactoring using some meta-function that determines if findBegin has to be called before accessing beginPosition
template <typename TFinder, typename TPattern>
static inline bool
myFindBegin(TFinder & finder, TPattern & pattern) {
    ignoreUnusedVariableWarning(finder);
    ignoreUnusedVariableWarning(pattern);
    return true;
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns, typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline bool
myFindBegin(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder,
        Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    return findBegin(finder, pattern);
}

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns, typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec, typename TScoreValue2>
inline bool
find(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder,
        Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern,
        const TScoreValue2 scoreLimit2) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;
    if(pattern._state == TPattern::STATE_EMPTY || pattern._state == TPattern::STATE_INITIAL) {
        setScoreLimit(pattern, scoreLimit2);
    } else {
        SEQAN_ASSERT_EQ_MSG(scoreLimit, scoreLimit(pattern), "The score limit cannot be changed after the initialization of the pattern.");
    }
    return find(finder, pattern);
}

// TODO(krugel) Cleanup comments

template <typename TText, typename TIndexSpec, typename TPieceFinderSpec, typename TVerifyFinderSpec, typename DoPreparePatterns, typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline bool
find(Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > & finder,
        Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    //typedef Finder<Index<TText, TIndexSpec>, FinderPartitioning<TPieceFinderSpec, TVerifyFinderSpec, DoPreparePatterns> > TFinder;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;
    typedef Index<TText, TIndexSpec>                        TIndex;
    typedef typename Position<TIndex>::Type                 TPosition;
    typedef typename Value<TScore>::Type                    TScoreValue;
    //typedef typename TFinder::TPieceFinder                  TPieceFinder;
    //typedef typename TFinder::TVerifyFinder                 TVerifyFinder;
    //typedef typename TPattern::TPiecePattern                TPiecePattern;
    //typedef typename Infix<TNeedle>::Type                   TNeedleInfix;
    //typedef typename Infix<TText>::Type                     TTextInfix;
    //typedef typename TPattern::TVerifyPattern               TVerifyPattern;

    //std::cerr << "findPiece(" << needle(pattern) << std::endl;

    if (DoPreparePatterns::VALUE) {
        setPreparedPiecePositions(finder._pieceFinder, finder._preparedPiecePositions);
    }
    
    SEQAN_ASSERT_GEQ(length(needle(pattern)), minSupportedPatternLength(finder, pattern));
    
    bool findPiece;
    if (finder._state == TPattern::STATE_EMPTY || pattern._state == TPattern::STATE_EMPTY) {
        return false;
    } else if (finder._state == TPattern::STATE_INITIAL || pattern._state == TPattern::STATE_INITIAL) {
        if (pattern._state != TPattern::STATE_INITIAL) clear(pattern);
        if (finder._state  != TPattern::STATE_INITIAL) clear(finder);
        if (pattern._numberOfPieces == 0u) _setDefaultNumberOfPieces(finder, pattern);
        _split(needle(pattern), pattern._numberOfPieces, pattern._pieces);

        //* finder._pieceFinder = TPieceFinder(haystack(finder)); // Don't overwrite the prepared positions
        //* pattern._piecePattern = TPiecePattern(pattern._pieces[pattern._pieceNo]);
        //- setHaystack(finder._pieceFinder, haystack(finder));
        setNeedle(pattern._piecePattern, pattern._pieces[pattern._pieceNo]);
        _setPieceScoreLimit(pattern);
        SEQAN_ASSERT_GEQ(length(needle(pattern._piecePattern)), minSupportedPatternLength(finder._pieceFinder, pattern._piecePattern));
        //std::cerr << pattern._pieceScoreLimit << " " << pattern._pieceScoreLimit << std::endl;
        setNeedle(pattern._verifyPattern, needle(pattern));
        setScoreLimit(pattern._verifyPattern, scoreLimit(pattern));
        findPiece = true;
         
        //std::cerr << "h = " << length(haystack(finder._pieceFinder)) << std::endl;
        //std::cerr << "n = " << needle(pattern._piecePattern) << std::endl;
        //std::cerr << "f = " << find(finder._pieceFinder, pattern._piecePattern) << std::endl;
        //std::cerr << "length(pieces) = " << length(pattern._pieces) << ", numberOfPieces = " << pattern._numberOfPieces << ", pieces[0] = " << pattern._pieces[0] << " * " << pattern._pieces[1] << std::endl;
    } else if (finder._state == TPattern::STATE_NOTFOUND) {
        pattern._state = TPattern::STATE_NOTFOUND;
        return false;
    } else if (pattern._scoreLimit == 0u) { // FOUND or BEGIN_FOUND or BEGIN_NOTFOUND
        findPiece = true;  // If we are doing exact search, there cannot be more than one possible match per verifyRegion
    } else { // FOUND or BEGIN_FOUND or BEGIN_NOTFOUND
        findPiece = false;
    }
    
    TScoreValue limit = -scoreLimit(pattern);
    while (true) {
        if (findPiece) {
            if (DoPreparePatterns::VALUE) {
                setPreparedPiecePositions(finder._pieceFinder, finder._preparedPiecePositions);
            }
        
            //std::cerr << "find()" << std::endl;
            //TPieceFinder fndr2(haystack(finder));
            //Pattern<TNeedle> ptrn(TNeedle(needle(pattern._piecePattern)));
            //while (find(fndr2, pattern._piecePattern)) {
            //    std::cerr << '+' << '[' << beginPosition(fndr2) << ',' << endPosition(fndr2) << ")\t" << infix(fndr2) << std::endl;
            //}
            //std::cerr << "findPiece(" << needle(pattern._piecePattern) << ") with " << pattern._pieceScoreLimit << std::endl;
            // In case of problems try to use "fndr" instead of "finder._pieceFinder": TPieceFinder fndr(haystack(finder));
            if (find(finder._pieceFinder, pattern._piecePattern)) {
                myFindBegin(finder._pieceFinder, pattern._piecePattern);
                #ifndef NDEBUG
                    ++finder._candidatePositions;
                #endif

                // TODO(krugel): This depends on minimal gap costs of -1
                // The end position can be #limit positions before or behind the hit.
                // If we verify an end position, we want to be sure to find all matches ending at this position.
                // Therefore we start 2*limit ahead.
                // This gives a region of total size: needleLength + 3 limit.
                // We guarantee to find all matches with: pos + needleLength - limit <= endPosition <= pos + needleLength + limit
                TPosition pieceOffset = beginPosition(pattern._pieces[pattern._pieceNo]) - beginPosition(needle(pattern));
                TPosition pos = beginPosition(finder._pieceFinder);
                TPosition begPos = (pos <= pieceOffset + 2u * limit) ? 0u : pos - pieceOffset - 2u * limit;
                TPosition endPos = pos - pieceOffset + length(needle(pattern)) + limit;
                endPos = std::min(endPos, length(haystack(finder)));
                finder._verifyRegion = infix(indexText(haystack(finder)), begPos, endPos);

	            //std::cerr << "found piece \"" << needle(pattern._piecePattern) << "\" in infix \"" << finder._verifyRegion << "\", pieceOffset = " << pieceOffset << ", limit = " << limit << std::endl;

                // Set up verification (if not possible by using "==")
                if (pattern._scoreLimit != 0u) {
                    // TO-DO(krugel): Avoid construction
                    clear(pattern._verifyPattern);
                    clear(finder._verifyFinder);
                    setHaystack(finder._verifyFinder, finder._verifyRegion);
                    //* pattern._verifyPattern = TVerifyPattern(needle(pattern));
                    //* finder._verifyFinder = TVerifyFinder(finder._verifyRegion);

                    // TODO(krugel): scoringScheme(pattern._verifyPattern), 
                    setScoreLimit(pattern._verifyPattern, pattern._scoreLimit);
                }

                findPiece = false;
                continue;
            } else {
                // When a piece is not found, continue with the next piece
                ++pattern._pieceNo;
                if (pattern._pieceNo == length(pattern._pieces)) {
                    finder._state  = TPattern::STATE_NOTFOUND;
                    pattern._state = TPattern::STATE_NOTFOUND;
                    return false;
                }
                clear(finder._pieceFinder);
                clear(pattern._piecePattern);
                setNeedle(pattern._piecePattern, pattern._pieces[pattern._pieceNo]);
                //finder._pieceFinder = TPieceFinder(haystack(finder));
                //pattern._piecePattern = TPiecePattern(pattern._pieces[pattern._pieceNo]);

                //_setPieceScoreLimit(pattern);
                // TODO(krugel): scoringScheme(pattern._piecePattern), 
                findPiece = true;
                continue;
            }
        } else {
            //std::cerr << "verify(" << haystack(finder._verifyFinder) << ")" << std::endl;
            
            if (pattern._scoreLimit == 0u) {
                // If we do exact search, there is no need to verify since we only have one piece which certainly matches
                //std::cerr << needle(pattern) << " != " << finder._verifyRegion << std::endl;
                //if (needle(pattern) != finder._verifyRegion) {
                //    std::cerr << needle(pattern) << " != " << finder._verifyRegion << std::endl;
                //    findPiece = true;
                //    continue;
                //}
                SEQAN_ASSERT_EQ(needle(pattern), finder._verifyRegion);
                finder._endPosition = beginPosition(finder._verifyRegion) - beginPosition(haystack(finder)) + length(needle(pattern));
                //std::cerr << "A finder._endPosition = " << finder._endPosition << std::endl;
                pattern._currentScore = 0u;
            } else {
                // Verify if the whole pattern matches (not only the piece)
                // TODO(krugel) Better don't verify the whole needle but everything without the piece?
                if (! find(finder._verifyFinder, pattern._verifyPattern)) {
                    findPiece = true;
                    continue;
                }
                //std::cerr << "a verified(" << haystack(finder._verifyFinder) << ") offset = " << beginPosition(pattern._pieces[pattern._pieceNo]) << " ## " << endPosition(finder._verifyFinder) << " < " << length(needle(pattern)) - limit << std::endl;

                // We must exclude matches that start too early (that do not include the found piece) because we maybe cannot find all begin positions
                //std::cerr << "skip? " << endPosition(finder._verifyFinder) << " < " << length(needle(pattern)) << " + " << limit << std::endl ;
                while (beginPosition(finder._verifyRegion) != 0 && endPosition(finder._verifyFinder) < length(needle(pattern)) + limit) {
                    if (!find(finder._verifyFinder, pattern._verifyPattern)) {
                        findPiece = true;
                        break;
                    }
                }
                if (findPiece) continue;
                finder._endPosition = beginPosition(finder._verifyRegion) - beginPosition(haystack(finder)) + endPosition(finder._verifyFinder);
                //std::cerr << "B finder._endPosition = " << finder._endPosition << std::endl;
                pattern._currentScore = getScore(pattern._verifyPattern);
            }

            //std::cerr << "b " << beginPosition(finder._verifyRegion) << " verified(" << haystack(finder._verifyFinder) << ") offset = " << beginPosition(pattern._pieces[pattern._pieceNo]) << " ## " << endPosition(finder._verifyFinder) << " < " << length(needle(pattern)) + limit << std::endl;
            finder._state  = TPattern::STATE_FOUND;
            pattern._state = TPattern::STATE_FOUND;
            // don't find the same position twice
            if (finder._returnedOcc.find(endPosition(finder)) == finder._returnedOcc.end()) {
                //std::cerr << "new" << std::endl;
                insert(finder._returnedOcc, endPosition(finder));
                return true;
            } else {
                //std::cerr << "duplicate" << std::endl;
                if (pattern._scoreLimit == 0u) findPiece = true; // Little hack to make this case work
                continue;
            }
        }
    }
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// scoringScheme, scoreLimit, getScore, beginScore
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline TScore const &
scoringScheme(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    return pattern._scoringScheme;
}

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline typename Value<TScore>::Type
scoreLimit(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    return pattern._scoreLimit;
}

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline typename Value<TScore>::Type
getScore(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;
    SEQAN_ASSERT(pattern._state == TPattern::STATE_FOUND
              || pattern._state == TPattern::STATE_BEGIN_FOUND
              || pattern._state == TPattern::STATE_BEGIN_NOTFOUND);
    //return getScore(pattern._verifyPattern);
    return pattern._currentScore;
}

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline typename Value<TScore>::Type
getBeginScore(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;
    SEQAN_ASSERT_EQ(pattern._state, TPattern::STATE_BEGIN_FOUND);
    if (pattern._scoreLimit == 0u) return 0u;
    return getBeginScore(pattern._verifyPattern);
}

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline typename Value<TScore>::Type
getBeginScore(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;
    SEQAN_ASSERT_EQ(pattern._state, TPattern::STATE_BEGIN_FOUND);
    if (pattern._scoreLimit == 0u) return 0u;
    return getBeginScore(pattern._verifyPattern);
}

// ----------------------------------------------------------------------------
// host, needle
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline TNeedle const &
host(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline TNeedle &
host(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline TNeedle const &
needle(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    return host(pattern);
}

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline TNeedle &
needle(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    return host(pattern);
}

// ----------------------------------------------------------------------------
// length
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline typename Size<TNeedle>::Type
length(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    return length(needle(pattern));
}

// ----------------------------------------------------------------------------
// begin, end, beginPosition, endPosition
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec, typename TTag>
inline typename Iterator<TNeedle const, Tag<TTag> const>::Type
begin(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;
    SEQAN_ASSERT_EQ(pattern._state, TPattern::STATE_BEGIN_FOUND);
    return begin(needle(pattern), spec);
}

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec, typename TTag>
inline typename Iterator<TNeedle const, Tag<TTag> const>::Type
end(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;
    SEQAN_ASSERT(pattern._state == TPattern::STATE_FOUND
              || pattern._state == TPattern::STATE_BEGIN_FOUND
              || pattern._state == TPattern::STATE_BEGIN_NOTFOUND);
    return end(needle(pattern), spec);
}

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline typename Position<TNeedle>::Type
beginPosition(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;
    SEQAN_ASSERT_EQ(pattern._state, TPattern::STATE_BEGIN_FOUND);
    return 0u;
}

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline typename Position<TNeedle>::Type
endPosition(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;
    SEQAN_ASSERT(pattern._state == TPattern::STATE_FOUND
              || pattern._state == TPattern::STATE_BEGIN_FOUND
              || pattern._state == TPattern::STATE_BEGIN_NOTFOUND);
    return length(needle(pattern));
}

// ----------------------------------------------------------------------------
// setScoreLimit, setScoringScheme
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec, typename TScoreValue2>
inline void
setScoreLimit(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern, TScoreValue2 scoreLimit2) {
    SEQAN_CHECKPOINT;
    clear(pattern);
    pattern._scoreLimit = scoreLimit2;
}

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec, typename TScore2>
inline void
setScoringScheme(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern, TScore2 const & scoringScheme2) {
    SEQAN_CHECKPOINT;
    clear(pattern);
    pattern._scoringScheme = scoringScheme2;
}

// ----------------------------------------------------------------------------
// setNumberOfPieces
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec, typename TPiecesSize>
inline void
setNumberOfPieces(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern, TPiecesSize numberOfPieces) {
    SEQAN_CHECKPOINT;
    clear(pattern);
    pattern._numberOfPieces = numberOfPieces;
}

// ----------------------------------------------------------------------------
// setHost
// ----------------------------------------------------------------------------

#ifdef SEQAN_CXX11_STANDARD

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec, typename TNeedle2>
inline void
setHost(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern, TNeedle2 && needle2) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;
    clear(pattern);
    setValue(pattern._host, std::forward<TNeedle2>(needle2));
    //pattern._host = Holder<TNeedle>(needle2);
    pattern._state = TPattern::STATE_INITIAL;
}

#else  // SEQAN_CXX11_STANDARD

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec, typename TNeedle2>
inline void
setHost(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern, TNeedle2 & needle2) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;
    clear(pattern);
    setValue(pattern._host, needle2);
    //pattern._host = Holder<TNeedle>(needle2);
    pattern._state = TPattern::STATE_INITIAL;
}

template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec, typename TNeedle2>
inline void
setHost(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern, TNeedle2 const & needle2) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;
    clear(pattern);
    setValue(pattern._host, needle2);
    pattern._state = TPattern::STATE_INITIAL;
}

#endif  // SEQAN_CXX11_STANDARD

// ----------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------

// Does not reset needle
template <typename TNeedle, typename TSpec, typename TScore, typename TPiecePatternSpec, typename TVerifyPatternSpec>
inline void
clear(Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> > TPattern;
    if (pattern._state != TPattern::STATE_EMPTY) pattern._state = TPattern::STATE_INITIAL;
    //pattern._currentScore = 0;
    pattern._numberOfPieces = 0u;
    pattern._pieceNo = 0u;
    clear(pattern._pieces);
    //clear(pattern._piecePattern);
    //clear(pattern._verifyPattern);
}

}  // namespace seqan

#endif  // SEQAN_FIND_INDEX_APPROX_FIND_INDEX_PARTITIONING_PATTERN_H_
