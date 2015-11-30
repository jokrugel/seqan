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
// Pattern class for the dynamic programming backtracking algorithm.
// ==========================================================================

#ifndef SEQAN_FIND_INDEX_APPROX_FIND_INDEX_BACKTRACKING_PATTERN_H_
#define SEQAN_FIND_INDEX_APPROX_FIND_INDEX_BACKTRACKING_PATTERN_H_

#include <algorithm> // std::min, std::max

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TNeedle, typename TScore>
struct Pattern<TNeedle, DPBacktracking<TScore> > : public FindState2_ {
    typedef Pattern<TNeedle, DPBacktracking<TScore> >       TPattern;
    typedef typename Value<TScore>::Type                    TScoreValue;
    typedef String<TScoreValue>                             TColumn;

    typename TPattern::TState _state;    // The pattern's state.
    Holder<TNeedle> _host;               // The needle we work on.
    TScoreValue _scoreLimit;             // The minimal score of a match.
    TScore _scoringScheme;               // The scoring scheme to use.
    TScoreValue _currentScore;           // The current score of a match.

    String<TColumn> _stack; // The matrix for the dynamic programming algorithm.

    Pattern()
        : _state(TPattern::STATE_EMPTY)
    {
        SEQAN_CHECKPOINT;
    }

    // TODO(krugel) Explicit definition should normally not be necessary?
    Pattern(TPattern const & rhs)
        : _state(rhs._state)
        , _host(rhs._host)
        , _scoreLimit(rhs._scoreLimit)
        , _scoringScheme(rhs._scoringScheme)
        , _currentScore(rhs._currentScore)
        , _stack(rhs._stack)
    {
        SEQAN_CHECKPOINT;
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 const & needle)
        : _state(TPattern::STATE_INITIAL)
        , _host(needle)
        , _scoreLimit(0)
        , _scoringScheme()
        , _stack()
    {
        SEQAN_CHECKPOINT;
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, TScoreValue scoreLimit)
        : _state(TPattern::STATE_INITIAL)
        , _host(needle)
        , _scoreLimit(scoreLimit)
        , _scoringScheme()
        , _stack()
    {
        SEQAN_CHECKPOINT;
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, TScoreValue scoreLimit, TScore scoringScheme)
        : _state(TPattern::STATE_INITIAL)
        , _host(needle)
        , _scoreLimit(scoreLimit)
        , _scoringScheme(scoringScheme)
        , _stack()
    {
        SEQAN_CHECKPOINT;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TSpec>
struct SupportsApproximateSearch<DPBacktracking<TSpec> > {
    typedef True Type;
};

template <typename TNeedle, typename TScore>
struct ScoringScheme<Pattern<TNeedle, DPBacktracking<TScore> > > {
    typedef TScore Type;
};

// ============================================================================
// Functions (private)
// ============================================================================

// returns true iff v implies a match
template <typename TScoreValue>
inline bool
_impliesMatch(String<TScoreValue> const & column, const TScoreValue scoreLimit) {
    SEQAN_CHECKPOINT;
    return back(column) >= scoreLimit;
}

// returns true iff v implies the dismissal of the subtree
template <typename TScoreValue>
inline bool
_impliesDismiss(String<TScoreValue> const & column, const TScoreValue scoreLimit) {
    SEQAN_CHECKPOINT;
    typedef String<TScoreValue>                             TColumn;
    typedef typename Position<TColumn>::Type                TPosition;

    for (TPosition i = 0u; i < length(column); ++i) {
        if (column[i] >= scoreLimit) return false;
    }
    return true;
}

// sets finder._data_iterator to the next node in the traversal and resizes the _stack respectively
template <typename TText, typename TIndexSpec, typename TNeedle, typename TScore>
inline bool
_flyNext(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder,
        Pattern<TNeedle, DPBacktracking<TScore> > & pattern) {
    SEQAN_CHECKPOINT;
    while (!goRight(finder._data_iterator)) {
        if (!goUp(finder._data_iterator)) return false;
    }
    resize(pattern._stack, parentRepLength(finder._data_iterator) + 1u);
    return true;
}

template <typename TNeedle, typename TScore>
inline void
_initColumn(Pattern<TNeedle, DPBacktracking<TScore> > & pattern) {
    SEQAN_CHECKPOINT;
    typedef typename Size<TNeedle>::Type                    TSize;
    typedef typename Value<TScore>::Type                    TScoreValue;
    typedef String<TScoreValue>                             TColumn;

    TColumn column;
    TSize lenPatt = length(needle(pattern));
    resize(column, lenPatt + 1u, Exact());
    for (TSize i = 0u; i <= lenPatt; ++i) column[i] = - (TScoreValue) i;
    appendValue(pattern._stack, column);
}

template <typename TNeedle, typename TScore, typename TTextInfix, typename TPosition>
inline void
_computeColumn(Pattern<TNeedle, DPBacktracking<TScore> > & pattern, TTextInfix const & text, TPosition pos) {
    SEQAN_CHECKPOINT;
    typedef typename Size<TNeedle>::Type                    TSize;
    typedef typename Value<TScore>::Type                    TScoreValue;
    typedef String<TScoreValue>                             TColumn;
    typedef typename SequenceEntryForScore<TScore, TNeedle>::Type TScoreEntryNeedle;
    typedef typename SequenceEntryForScore<TScore, TTextInfix>::Type TScoreEntryText;

    TNeedle & ndl = needle(pattern);
    TSize lenPatt = length(ndl);
    TSize lenPath = length(pattern._stack);
    TColumn column;
    TScoreValue a, b, c;

    resize(column, lenPatt + 1u, Exact());
    column[0u] = - (TScoreValue) lenPath;

    TScoreEntryText scoreEntryText = sequenceEntryForScore(scoringScheme(pattern), text, pos);
    for (TSize i = 1u; i <= lenPatt; ++i) {
        TScoreEntryNeedle scoreEntryNeedle = sequenceEntryForScore(scoringScheme(pattern), ndl, i - 1u);
        a = pattern._stack[lenPath - 1u][i - 1u] + score                 (scoringScheme(pattern), scoreEntryText, scoreEntryNeedle);
        b = pattern._stack[lenPath - 1u][i]      + scoreGapOpenHorizontal(scoringScheme(pattern), scoreEntryText, scoreEntryNeedle);
        c = column[i - 1u]                       + scoreGapOpenVertical  (scoringScheme(pattern), scoreEntryText, scoreEntryNeedle);
        column[i] = std::max(a, std::max(b, c));
    }

    appendValue(pattern._stack, column);
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// findBegin
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TNeedle, typename TScore, typename TScoreValue>
inline bool
findBegin(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder,
        Pattern<TNeedle, DPBacktracking<TScore> > & pattern, TScoreValue scoreLimit2) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPBacktracking<TScore> >       TPattern;
    SEQAN_ASSERT(pattern._state == TPattern::STATE_FOUND
              || pattern._state == TPattern::STATE_BEGIN_FOUND);
    
    if (pattern._state == TPattern::STATE_BEGIN_FOUND) {
        finder._state  =  TPattern::STATE_BEGIN_NOTFOUND;
        pattern._state =  TPattern::STATE_BEGIN_NOTFOUND;
        return false;
    } else {
        if (getScore(pattern) >= scoreLimit2) {
            finder._state  =  TPattern::STATE_BEGIN_FOUND;
            pattern._state =  TPattern::STATE_BEGIN_FOUND;
            return true;
        } else {
            finder._state  =  TPattern::STATE_BEGIN_NOTFOUND;
            pattern._state =  TPattern::STATE_BEGIN_NOTFOUND;
            return false;
        }
    }
}

template <typename TText, typename TIndexSpec, typename TNeedle, typename TScore>
inline bool
findBegin(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder,
        Pattern<TNeedle, DPBacktracking<TScore> > & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPBacktracking<TScore> >       TPattern;
    SEQAN_ASSERT(pattern._state == TPattern::STATE_FOUND
                   || pattern._state == TPattern::STATE_BEGIN_FOUND);
    
    if (pattern._state == TPattern::STATE_BEGIN_FOUND) {
        finder._state  =  TPattern::STATE_BEGIN_NOTFOUND;
        pattern._state =  TPattern::STATE_BEGIN_NOTFOUND;
        return false;
    } else {
        finder._state  =  TPattern::STATE_BEGIN_FOUND;
        pattern._state =  TPattern::STATE_BEGIN_FOUND;
        return true;
    }
}

// ----------------------------------------------------------------------------
// find
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TNeedle, typename TScore, typename TScoreValue2>
inline bool
find(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder,
        Pattern<TNeedle, DPBacktracking<TScore> > & pattern, const TScoreValue2 scoreLimit2) {
    SEQAN_CHECKPOINT;
    setScoreLimit(pattern, scoreLimit2);
    return find(finder, pattern);
}

template <typename TText, typename TIndexSpec, typename TNeedle, typename TScore>
inline bool
find(Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > & finder,
        Pattern<TNeedle, DPBacktracking<TScore> > & pattern) {
    SEQAN_CHECKPOINT;
    typedef Finder<Index<TText, TIndexSpec>, DPBacktracking<TScore> > TFinder;
    typedef Pattern<TNeedle, DPBacktracking<TScore> >       TPattern;
    typedef typename Position<Index<TText, TIndexSpec> >::Type TPosition;
    typedef typename Size<Index<TText, TIndexSpec> >::Type  TSize;
    typedef typename Value<TScore>::Type                    TScoreValue;
    typedef typename TFinder::TIterator                     TIterator;
    
    TIterator & it = finder._data_iterator;
    bool nextNodeExists;
    TScoreValue limit = - scoreLimit(pattern);
    TSize maxLen = length(needle(pattern)) + limit + 1u;

    if (finder._state == TPattern::STATE_EMPTY || pattern._state == TPattern::STATE_EMPTY) {
        return false;
    } else if (finder._state == TPattern::STATE_INITIAL || pattern._state == TPattern::STATE_INITIAL) {
        if (pattern._state != TPattern::STATE_INITIAL) clear(pattern);
        if (finder._state  != TPattern::STATE_INITIAL) clear(finder);

        _initColumn(pattern);
        nextNodeExists = goDown(it);
        if (! nextNodeExists) {
            finder._state  = TPattern::STATE_NOTFOUND;
            pattern._state = TPattern::STATE_NOTFOUND;
            return false;
        }
    } else if (finder._state == TPattern::STATE_NOTFOUND) {
        pattern._state = TPattern::STATE_NOTFOUND;
        return false;
    } else { // FOUND or BEGIN_FOUND or BEGIN_NOTFOUND
        ++finder._occIt; // the next occurrence of the current suffix tree node
        if (! atEnd(finder._occIt, finder._occurrences)) {
            finder._state  = TPattern::STATE_FOUND;
            pattern._state = TPattern::STATE_FOUND;
            // The length of the current occurrence is the same as the length of the previous,
            // because they come from the same suffix tree node
            finder._endPosition -= finder._beginPosition;
            finder._beginPosition = *finder._occIt;
            finder._endPosition += finder._beginPosition;
            return true;
        }
    }
    
    while (true) { // search for the next occurrence
        // TODO(krugel): This depends on minimal gap costs of -1
        TPosition currentPos = length(pattern._stack) - parentRepLength(it) - 1u;
        if (length(pattern._stack) + 1u > maxLen) { // Cutoff the tree below
            nextNodeExists = _flyNext(finder, pattern);
            if (! nextNodeExists) break;
        } else if (currentPos < parentEdgeLength(it)) {
            // Take the character at position currentPos of the current suffix tree edge and compute the dynamic programming matrix
            _computeColumn(pattern, parentEdgeLabel(it), currentPos);

            if (_impliesMatch(back(pattern._stack), scoreLimit(pattern))) {
                finder._state  = TPattern::STATE_FOUND;
                pattern._state = TPattern::STATE_FOUND;
                finder._occurrences = getOccurrences(it);
                finder._occIt = begin(finder._occurrences);
                finder._beginPosition = *finder._occIt;
                // The length of the occurrence depends on the current depth in the suffix tree
                // plus the length of the edge we are going down at the moment
                finder._endPosition = finder._beginPosition + length(pattern._stack) - 1u;
                pattern._currentScore = back(back(pattern._stack));
                return true;
            } else if (_impliesDismiss(back(pattern._stack), scoreLimit(pattern))) {
                nextNodeExists = _flyNext(finder, pattern);
                if (! nextNodeExists) break;
            }
        } else {
            nextNodeExists = goDown(it) || _flyNext(finder, pattern);
            if (! nextNodeExists) break;
        }
    }
    finder._state  = TPattern::STATE_NOTFOUND;
    pattern._state = TPattern::STATE_NOTFOUND;
    return false;
}

// ----------------------------------------------------------------------------
// scoringScheme, scoreLimit, getScore, beginScore
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TScore>
inline TScore const &
scoringScheme(Pattern<TNeedle, DPBacktracking<TScore> > const & pattern) {
    SEQAN_CHECKPOINT;
    return pattern._scoringScheme;
}

template <typename TNeedle, typename TScore>
inline typename Value<TScore>::Type
scoreLimit(Pattern<TNeedle, DPBacktracking<TScore> > const & pattern) {
    SEQAN_CHECKPOINT;
    return pattern._scoreLimit;
}

template <typename TNeedle, typename TScore>
inline typename Value<TScore>::Type
getScore(Pattern<TNeedle, DPBacktracking<TScore> > const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPBacktracking<TScore> >       TPattern;
    SEQAN_ASSERT(pattern._state == TPattern::STATE_FOUND
              || pattern._state == TPattern::STATE_BEGIN_FOUND
              || pattern._state == TPattern::STATE_BEGIN_NOTFOUND);
    return pattern._currentScore;
}

template <typename TNeedle, typename TScore>
inline typename Value<TScore>::Type
getBeginScore(Pattern<TNeedle, DPBacktracking<TScore> > const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPBacktracking<TScore> >       TPattern;
    SEQAN_ASSERT_EQ(pattern._state, TPattern::STATE_BEGIN_FOUND);
    return getScore(pattern);
}

// ----------------------------------------------------------------------------
// host, needle
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TScore>
inline TNeedle const &
host(Pattern<TNeedle, DPBacktracking<TScore> > const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}

template <typename TNeedle, typename TScore>
inline TNeedle &
host(Pattern<TNeedle, DPBacktracking<TScore> > & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}

template <typename TNeedle, typename TScore>
inline TNeedle const &
needle(Pattern<TNeedle, DPBacktracking<TScore> > const & pattern) {
    SEQAN_CHECKPOINT;
    return host(pattern);
}

template <typename TNeedle, typename TScore>
inline TNeedle &
needle(Pattern<TNeedle, DPBacktracking<TScore> > & pattern) {
    SEQAN_CHECKPOINT;
    return host(pattern);
}

// ----------------------------------------------------------------------------
// length
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TScore>
inline typename Size<TNeedle>::Type
length(Pattern<TNeedle, DPBacktracking<TScore> > const & pattern) {
    SEQAN_CHECKPOINT;
    return length(needle(pattern));
}

// ----------------------------------------------------------------------------
// begin, end, beginPosition, endPosition
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TScore, typename TTag>
inline typename Iterator<TNeedle const, Tag<TTag> const>::Type
begin(Pattern<TNeedle, DPBacktracking<TScore> > const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPBacktracking<TScore> >       TPattern;
    SEQAN_ASSERT_EQ(pattern._state, TPattern::STATE_BEGIN_FOUND);
    return begin(needle(pattern), spec);
}

template <typename TNeedle, typename TScore, typename TTag>
inline typename Iterator<TNeedle const, Tag<TTag> const>::Type
end(Pattern<TNeedle, DPBacktracking<TScore> > const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPBacktracking<TScore> >       TPattern;
    SEQAN_ASSERT(pattern._state == TPattern::STATE_FOUND
              || pattern._state == TPattern::STATE_BEGIN_FOUND);
    return end(needle(pattern), spec);
}

template <typename TNeedle, typename TScore>
inline typename Position<TNeedle>::Type
beginPosition(Pattern<TNeedle, DPBacktracking<TScore> > const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPBacktracking<TScore> >       TPattern;
    SEQAN_ASSERT_EQ(pattern._state, TPattern::STATE_BEGIN_FOUND);
    return 0u;
}

template <typename TNeedle, typename TScore>
inline typename Position<TNeedle>::Type
endPosition(Pattern<TNeedle, DPBacktracking<TScore> > const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPBacktracking<TScore> >       TPattern;
    SEQAN_ASSERT(pattern._state == TPattern::STATE_FOUND
              || pattern._state == TPattern::STATE_BEGIN_FOUND);
    return length(needle(pattern));
}

// ----------------------------------------------------------------------------
// setScoreLimit, setScoringScheme
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TScore, typename TScoreValue2>
inline void
setScoreLimit(Pattern<TNeedle, DPBacktracking<TScore> > & pattern, TScoreValue2 scoreLimit) {
    SEQAN_CHECKPOINT;
    clear(pattern);
    pattern._scoreLimit = scoreLimit;
}

template <typename TNeedle, typename TScore, typename TScore2>
inline void
setScoringScheme(Pattern<TNeedle, DPBacktracking<TScore> > & pattern, TScore2 const & scoringScheme2) {
    SEQAN_CHECKPOINT;
    clear(pattern);
    pattern._scoringScheme = scoringScheme2;
}

// ----------------------------------------------------------------------------
// setHost
// ----------------------------------------------------------------------------

#ifdef SEQAN_CXX11_STANDARD

template <typename TNeedle, typename TScore, typename TNeedle2>
inline void
setHost(Pattern<TNeedle, DPBacktracking<TScore> > & pattern, TNeedle2 && needle2) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPBacktracking<TScore> >       TPattern;
    clear(pattern);
    setValue(pattern._host, std::forward<TNeedle2>(needle2));
    pattern._state = TPattern::STATE_INITIAL;
}

#else  // SEQAN_CXX11_STANDARD

template <typename TNeedle, typename TScore, typename TNeedle2>
inline void
setHost(Pattern<TNeedle, DPBacktracking<TScore> > & pattern, TNeedle2 & needle2) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPBacktracking<TScore> >       TPattern;
    clear(pattern);
    setValue(pattern._host, needle2);
    pattern._state = TPattern::STATE_INITIAL;
}

template <typename TNeedle, typename TScore, typename TNeedle2>
inline void
setHost(Pattern<TNeedle, DPBacktracking<TScore> > & pattern, TNeedle2 const & needle2) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPBacktracking<TScore> >       TPattern;
    clear(pattern);
    setValue(pattern._host, needle2);
    pattern._state = TPattern::STATE_INITIAL;
}
#endif  // SEQAN_CXX11_STANDARD

// ----------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TScore>
inline void
clear(Pattern<TNeedle, DPBacktracking<TScore> > & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPBacktracking<TScore> >       TPattern;
    if (pattern._state != TPattern::STATE_EMPTY) pattern._state = TPattern::STATE_INITIAL;
    pattern._currentScore = 0;
    clear(pattern._stack);
}

}  // namespace seqan

#endif  // SEQAN_FIND_INDEX_APPROX_FIND_INDEX_BACKTRACKING_PATTERN_H_
