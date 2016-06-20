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
// Finder class for q-gram index.
// Works with IndexQGram and IndexQGram2L and is able to also search for
// patterns that are shorter or longer than the gram size Q.
// ==========================================================================

#ifndef SEQAN_INDEX_FIND_INDEX_QGRAM_EXT_H_
#define SEQAN_INDEX_FIND_INDEX_QGRAM_EXT_H_

// For metafunction GetOccurrences
#include <seqan/find_index_approx/find_index_approx_base.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// OccurrencesIterator
// ----------------------------------------------------------------------------

template <typename TIndex>
struct OccurrencesIterator{};

template <typename TObject, typename TShapeSpec, typename TSpec>
struct OccurrencesIterator<Index<TObject, IndexQGram<TShapeSpec, TSpec> > > {
    typedef Index<TObject, IndexQGram<TShapeSpec, TSpec> >  TIndex;
    typedef typename Position<TObject>::Type                TPosition;
    typedef typename GetOccurrences<TIndex>::Type           TOccurrences;
    typedef typename Iterator<TOccurrences>::Type           TOccurrencesIterator;
    
    //Holder<TIndex> index;
    TOccurrences occ;
    TOccurrencesIterator realOccIt;
    
    OccurrencesIterator() {}
    
    template <typename TValue, typename TShapeSpec2>
    void init(TIndex & index2, Shape<TValue, TShapeSpec2> const &shape, typename Value<Shape<TValue, TShapeSpec2> >::Type upperHash) {
        //setValue(index, index2);
        occ = getOccurrences(index2, shape, upperHash);
        realOccIt = begin(occ);
    }
    
    TPosition operator*() {
        return value(realOccIt);
    }

    void operator++() {
        goNext(realOccIt);
    }

    template <typename TNeedle, typename TNeedleSize>
    void moveForward(TNeedle & wholeNeedle, TNeedleSize gramOffsetInWholeNeedle) {
        ignoreUnusedVariableWarning(wholeNeedle);
        ignoreUnusedVariableWarning(gramOffsetInWholeNeedle);
        goNext(realOccIt);
    }
};

template <typename TObject, unsigned int Q, typename TAddressingSpec, typename TSuffixTreeSpec>
struct OccurrencesIterator<Index<TObject, IndexQGram2L<UngappedShape<Q>, TAddressingSpec, TSuffixTreeSpec> > > {
    typedef Index<TObject, IndexQGram2L<UngappedShape<Q>, TAddressingSpec, TSuffixTreeSpec> > TIndex;

    typedef typename Position<TObject>::Type                TPosition;
    typedef typename Size<TObject>::Type                    TSize;
    typedef typename Fibre<TIndex, QGram2LFrontEndIndex>::Type TFrontEndIndex;
    typedef typename Fibre<TIndex, QGram2LBackEndIndex>::Type TBackEndIndex;
    typedef typename Fibre<TFrontEndIndex, QGramDir>::Type  TFrontEndIndexDir; // String<Size>
    typedef typename Fibre<TFrontEndIndex, QGramSA>::Type   TFrontEndIndexSA;  // String<SAValue>
    //typedef typename Fibre<TBackEndIndex, QGramDir>::Type   TBackEndIndexDir;  // String<Size>
    //typedef typename Fibre<TBackEndIndex, QGramSA>::Type    TBackEndIndexSA;   // String<SAValue>
    typedef typename Position<TFrontEndIndexDir>::Type      TFrontEndIndexDirPosition;
    typedef typename Position<TFrontEndIndexSA>::Type       TFrontEndIndexSAPosition;
    typedef typename SAValue<TFrontEndIndex>::Type          TFrontEndIndexSAValue;
    typedef typename Position<TBackEndIndex>::Type          TBackEndIndexDirPosition;
    typedef typename Position<TBackEndIndex>::Type          TBackEndIndexSAPosition;

    Holder<TIndex> index;
    TFrontEndIndexSAPosition saFrontPos;
    TFrontEndIndexSAPosition saFrontPosEnd;
    TBackEndIndexSAPosition saBackPos;
    TBackEndIndexSAPosition saBackPosEnd;
    
    OccurrencesIterator() {}
    
    void nextSubsequence() {
        TFrontEndIndexSAValue saValue = saAt(saFrontPos, value(index).frontEndIndex);
        TBackEndIndexDirPosition seqNo = getSeqNo(saValue);
        saBackPos = dirAt(seqNo, value(index).backEndIndex);
        saBackPosEnd = dirAt(seqNo + 1u, value(index).backEndIndex);
    }

    template <typename TValue, typename TShapeSpec2>
    void init(TIndex & index2, Shape<TValue, TShapeSpec2> const &shape, typename Value<Shape<TValue, TShapeSpec2> >::Type upperHash) {
        setValue(index, index2);
        indexRequire(value(index), QGram2LBackEndIndex());
        indexRequire(value(index), QGram2LFrontEndIndex());
        TFrontEndIndexDirPosition hashValue = value(shape);
        hashValue = getBucket(indexBucketMap(value(index).frontEndIndex), hashValue);
        upperHash = getBucket(indexBucketMap(value(index).frontEndIndex), upperHash);
        saFrontPos = dirAt(hashValue, value(index).frontEndIndex);
        saFrontPosEnd = dirAt(upperHash, value(index).frontEndIndex);
        if (saFrontPos >= saFrontPosEnd) return;
        nextSubsequence();
    }
    
    TPosition operator*() {
        // position of subsequence in text + offest of gram within subsequence
        return saAt(saBackPos, value(index).backEndIndex) + getSeqOffset(saAt(saFrontPos, value(index).frontEndIndex));
    }

    void operator++() {
        ++saBackPos;
        if (saBackPos >= saBackPosEnd) {
            ++saFrontPos;
            if (saFrontPos >= saFrontPosEnd) return;
            nextSubsequence();
        }
    }

    // Do the same as operator++, but skip subsequences that 
    template <typename TNeedle, typename TNeedleSize>
    void moveForward(TNeedle & wholeNeedle, TNeedleSize gramOffsetInWholeNeedle) {
        typedef typename Infix<TNeedle>::Type               TNeedleInfix;
        typedef typename Infix<TObject>::Type               TTextInfix;

        ignoreUnusedVariableWarning(wholeNeedle);
        ignoreUnusedVariableWarning(gramOffsetInWholeNeedle);
        ++saBackPos;
        if (saBackPos >= saBackPosEnd) {
            while(true) {
                ++saFrontPos;
                if (saFrontPos >= saFrontPosEnd) return;
                nextSubsequence();

                // Skip currenct subsequence if overlap with needle is not equal
                TSize needleLength = length(wholeNeedle);
                if (needleLength <= Q) return;

                TSize subsequenceLength = getSubsequenceLength(value(index));
                TSize subsequenceOffset = getSeqOffset(saAt(saFrontPos, value(index).frontEndIndex));
                TSize commonOffset = std::min(gramOffsetInWholeNeedle, subsequenceOffset);
                TSize commonLength = std::min(subsequenceLength - subsequenceOffset, needleLength - gramOffsetInWholeNeedle) + commonOffset;
                if (commonLength <= Q) return;
                TSize subsequenceStart = saAt(saBackPos, value(index).backEndIndex);
                TSize textStart = subsequenceStart + subsequenceOffset - commonOffset;
                TSize needleStart = 0u + gramOffsetInWholeNeedle - commonOffset;
                
                TTextInfix textInfix = infixWithLength(indexText(value(index)), textStart, commonLength);
                TNeedleInfix ndlInfix = infixWithLength(wholeNeedle, needleStart, commonLength);
                //std::cerr << "ndle = " << wholeNeedle << std::endl;
                //std::cerr << "gram = " << infixWithLength(wholeNeedle, gramOffsetInWholeNeedle, 3)<< std::endl;
                //std::cerr << "ndleInfix = " << ndlInfix << std::endl;
                //std::cerr << "textInfix = " << textInfix << std::endl;
                if (textInfix == ndlInfix) return;
            }
        }
    }
};

// ----------------------------------------------------------------------------
// Finder
// ----------------------------------------------------------------------------

template<typename TVerification = True, typename TSpec = void>
struct FinderQGramExtended {};

// TVerification = True | False
template <typename TText, typename TIndexSpec, typename TVerification, typename TSpecFinder>
class Finder<Index<TText, TIndexSpec>, FinderQGramExtended<TVerification, TSpecFinder> > {
protected:
    typedef Index<TText, TIndexSpec>                        TIndex;
    typedef typename Position<TText>::Type                  TPosition;
    typedef typename Size<TIndex>::Type						TSize;
    //typedef typename GetOccurrences<TIndex>::Type           TOccurrences;
    //typedef typename Iterator<TOccurrences>::Type           TOccurrencesIterator;
    typedef OccurrencesIterator<TIndex>                     TOccurrencesIterator;

public:
    Holder<TIndex>  index;
    TSize           data_length;
    
    bool first;
    //TOccurrences occurrences;
    TOccurrencesIterator occurrencesIter;
    TPosition currentStep;  // Value in the range [0 .. stepSize(index) - 1]
    TPosition posChoice; // Value is multiple of stepSize(index)

    Finder() {
        clear(*this);
    }

    Finder(TIndex &_index): index(_index) {
        clear(*this);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TText, typename TShapeSpec, typename TSuffixTreeSpec, typename TSpec>
struct DefaultFinder<Index<TText, IndexQGram2L<TShapeSpec, TSuffixTreeSpec, TSpec> > > {
    typedef FinderQGramExtended<True>                       Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// atEnd(OccurrencesIter)
// ----------------------------------------------------------------------------

template <typename TObject, typename TShapeSpec, typename TSpec>
inline bool
atEnd(OccurrencesIterator<Index<TObject, IndexQGram<TShapeSpec, TSpec> > > & occIt) {
    return atEnd(occIt.realOccIt, occIt.occ);
}

template <typename TObject, unsigned int Q, typename TAddressingSpec, typename TSuffixTreeSpec>
inline bool
atEnd(OccurrencesIterator<Index<TObject, IndexQGram2L<UngappedShape<Q>, TAddressingSpec, TSuffixTreeSpec> > > & occIt) {
    return occIt.saFrontPos >= occIt.saFrontPosEnd || occIt.saBackPos >= occIt.saBackPosEnd;
}

// ----------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TVerification, typename TSpecFinder>
inline void
clear(Finder<Index<TText, TIndexSpec>, FinderQGramExtended<TVerification, TSpecFinder> > & me) {
    me.first = true;
    me.data_length = 0;
    me.currentStep = 0;
    me.posChoice = 0;
    //me.occurrencesIter = begin(me.occurrences);
    //clear(me.occurrences);
}

// ----------------------------------------------------------------------------
// empty
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TVerification, typename TSpecFinder>
inline bool
empty(Finder<Index<TText, TIndexSpec>, FinderQGramExtended<TVerification, TSpecFinder> > & me) {
    return me.first;
}

// ----------------------------------------------------------------------------
// minSupportedPatternLength
// ----------------------------------------------------------------------------

/*
    qgram*: restriction with open addressing
        numberOfPieces = tol + 1
        pieceLength = floor(patternLength / numberOfPieces)
        pieceLength >= q + stepSize - 1
        ==> floor(patternLength / (tol + 1)) >= q + stepSize + 1   (necessary to find all matches)

    qgram*:
        pieceLength >= q + stepSize - 1                            (for efficiency to find q-grams and not shorter)

    qsample: restriction with sample
        pieceLength >= q + stepSize - 1                            (necessary to find all matches)
*/

template <typename TText, typename TShapeSpec, typename TSpec, typename TVerification, typename TSpecFinder, typename TPattern>
inline typename Size<Finder<Index<TText, IndexQGram<TShapeSpec, TSpec> >, FinderQGramExtended<TVerification, TSpecFinder> > >::Type
minSupportedPatternLength(Finder<Index<TText, IndexQGram<TShapeSpec, TSpec> >, FinderQGramExtended<TVerification, TSpecFinder> > & finder, TPattern & /*pattern*/) {
    // For a regular q-gram index no restriction (we can also search patterns < Q).
    if (getStepSize(haystack(finder)) == 1u) return 0u;

    // For a q-sample index a match has to span at least one sampled q-gram.
    typedef typename Size<Finder<Index<TText, IndexQGram<TShapeSpec, OpenAddressing> >, FinderQGramExtended<TVerification, TSpecFinder> > >::Type TSize;
    const TSize Q = length(indexShape(haystack(finder)));
    return Q + getStepSize(haystack(finder)) - 1u;
}

template <typename TText, typename TShapeSpec, typename TVerification, typename TSpecFinder, typename TPattern>
inline typename Size<Finder<Index<TText, IndexQGram<TShapeSpec, OpenAddressing> >, FinderQGramExtended<TVerification, TSpecFinder> > >::Type
minSupportedPatternLength(Finder<Index<TText, IndexQGram<TShapeSpec, OpenAddressing> >, FinderQGramExtended<TVerification, TSpecFinder> > & finder, TPattern & /*pattern*/) {
    // When using open addressing, the pattern has to be of length at least Q, even with stepSize = 1.
    typedef typename Size<Finder<Index<TText, IndexQGram<TShapeSpec, OpenAddressing> >, FinderQGramExtended<TVerification, TSpecFinder> > >::Type TSize;
    const TSize Q = length(indexShape(haystack(finder)));
    return Q + getStepSize(haystack(finder)) - 1u;
}

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TVerification, typename TSpecFinder, typename TPattern>
inline typename Size<Finder<Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> >, FinderQGramExtended<TVerification, TSpecFinder> > >::Type
minSupportedPatternLength(Finder<Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> >, FinderQGramExtended<TVerification, TSpecFinder> > & /*finder*/, TPattern & /*pattern*/) {
    // q-gram/2L index does not support q-samples.
    return 0u;
}

template <typename TText, typename TShapeSpec, typename TSuffixTreeSpec, typename TVerification, typename TSpecFinder, typename TPattern>
inline typename Size<Finder<Index<TText, IndexQGram2L<TShapeSpec, OpenAddressing, TSuffixTreeSpec> >, FinderQGramExtended<TVerification, TSpecFinder> > >::Type
minSupportedPatternLength(Finder<Index<TText, IndexQGram2L<TShapeSpec, OpenAddressing, TSuffixTreeSpec> >, FinderQGramExtended<TVerification, TSpecFinder> > & finder, TPattern & /*pattern*/) {
    typedef typename Size<Finder<Index<TText, IndexQGram2L<TShapeSpec, OpenAddressing, TSuffixTreeSpec> >, FinderQGramExtended<TVerification, TSpecFinder> > >::Type TSize;
    const TSize Q = length(indexShape(haystack(finder)));
    return Q + getStepSize(haystack(finder)) - 1u;
}

// ----------------------------------------------------------------------------
// find
// ----------------------------------------------------------------------------

// TODO(krugel) Move to index_qgram.h
template <typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue>
inline typename Infix<typename Fibre<Index<TText, IndexQGram<TShapeSpec, TSpec> >, FibreSA>::Type const>::Type 
getOccurrences(Index<TText, IndexQGram<TShapeSpec, TSpec> > const &index, Shape<TValue, TShapeSpec2> const &shape, typename Value<Shape<TValue, TShapeSpec2> >::Type upperHash) {
    typedef typename Size<typename Fibre<Index<TText, IndexQGram<TShapeSpec, TSpec> >, FibreDir>::Type>::Type TDirSize;
    TDirSize bucket = getBucket(indexBucketMap(index), value(shape));
    TDirSize bucket2 = getBucket(indexBucketMap(index), upperHash);
    return infix(indexSA(index), indexDir(index)[bucket], indexDir(index)[bucket2]);
}

template <typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue>
inline typename Infix<typename Fibre<Index<TText, IndexQGram<TShapeSpec, TSpec> >, FibreSA>::Type const>::Type 
getOccurrences(Index<TText, IndexQGram<TShapeSpec, TSpec> > &index, Shape<TValue, TShapeSpec2> const &shape, typename Value<Shape<TValue, TShapeSpec2> >::Type upperHash) {
    typedef typename Size<typename Fibre<Index<TText, IndexQGram<TShapeSpec, TSpec> >, FibreDir>::Type>::Type TDirSize;
    indexRequire(index, QGramSADir());
    TDirSize bucket = getBucket(indexBucketMap(index), value(shape));
    TDirSize bucket2 = getBucket(indexBucketMap(index), upperHash);
    return infix(indexSA(index), indexDir(index)[bucket], indexDir(index)[bucket2]);
}

template <typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue>
inline typename Size<Index<TText, IndexQGram<TShapeSpec, TSpec> > const>::Type
countOccurrences(Index<TText, IndexQGram<TShapeSpec, TSpec> > const &index, Shape<TValue, TShapeSpec2> const &shape, typename Value<Shape<TValue, TShapeSpec2> >::Type upperHash) {
    typedef typename Size<typename Fibre<Index<TText, IndexQGram<TShapeSpec, TSpec> >, FibreDir>::Type>::Type TDirSize;
    TDirSize bucket = getBucket(indexBucketMap(index), value(shape));
    TDirSize bucket2 = getBucket(indexBucketMap(index), upperHash);
    return indexDir(index)[bucket2] - indexDir(index)[bucket];
}

template <typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue>
inline typename Size<Index<TText, IndexQGram<TShapeSpec, TSpec> > const>::Type 
countOccurrences(Index<TText, IndexQGram<TShapeSpec, TSpec> > &index, Shape<TValue, TShapeSpec2> const &shape, typename Value<Shape<TValue, TShapeSpec2> >::Type upperHash) {
    typedef typename Size<typename Fibre<Index<TText, IndexQGram<TShapeSpec, TSpec> >, FibreDir>::Type>::Type TDirSize;
    indexRequire(index, QGramSADir());
    TDirSize bucket = getBucket(indexBucketMap(index), value(shape));
    TDirSize bucket2 = getBucket(indexBucketMap(index), upperHash);
    return indexDir(index)[bucket2] - indexDir(index)[bucket];
}

template <typename TText, typename TIndexSpec, typename TVerification, typename TSpecFinder, typename TNeedle, typename TPatternSpec>
inline bool
find(Finder<Index<TText, TIndexSpec>, FinderQGramExtended<TVerification, TSpecFinder> > &finder, Pattern<TNeedle, TPatternSpec> const &pattern) {
    typedef Index<TText, TIndexSpec>                        TIndex;
    typedef typename Size<TText>::Type                      TSize;
    typedef typename Fibre<TIndex, FibreShape>::Type        TShape;
    typedef typename Value<TShape>::Type                    TShapeValue;
    
    TIndex & index = haystack(finder);
    
    const TSize Q = length(indexShape(index));
    
    // The finder also works for patterns shorter than Q.
    // SEQAN_ASSERT_GEQ(length(needle(pattern)), Q);

    // If the assertion fails, not the whole text will be indexed and we might miss some occurrences.
    // But the finder should still work in this case.
    // SEQAN_ASSERT_LEQ(getStepSize(index), Q);

    // If the assertion fails, we might miss some occurrences.
    // But the finder should still work in this case.
    // SEQAN_ASSERT_GEQ(length(needle(pattern)), Q + getStepSize(index) - 1u);
    SEQAN_ASSERT_GEQ(length(needle(pattern)), minSupportedPatternLength(finder, pattern));

    while(true) {
        bool nextStep = false;
        if (finder.first) {
            nextStep = true;
            finder.first = false;
            _setFinderLength(finder, length(needle(pattern)));
        } else if (atEnd(finder.occurrencesIter)) {
            ++finder.currentStep;
            if (finder.currentStep >= getStepSize(index) || finder.currentStep >= length(needle(pattern)) - Q) return false;
            nextStep = true;
        }
        
        if (nextStep) {
            TShapeValue upperHash;
            TSize occCount = 0;
            TSize best = MaxValue<TSize>::VALUE;
            finder.posChoice = 0;
            
            // If needle is longer than Q, try different starting positions and choose the gram with the shortest positing list
            // Limit choice to first few possibilities? Is not faster  // && i < 2 * Q * getStepSize(index)
            if (length(needle(pattern)) > Q) {
                for (TSize i = 0u; i + finder.currentStep + Q <= length(needle(pattern)); i += getStepSize(index)) {
                    hashUpper(indexShape(index), begin(needle(pattern)) + i + finder.currentStep, length(needle(pattern)) - i - finder.currentStep);
                    upperHash = value(indexShape(index));
                    hash(indexShape(index), begin(needle(pattern)) + i + finder.currentStep, length(needle(pattern)) - i - finder.currentStep);
                    occCount = countOccurrences(index, indexShape(index), upperHash);
                    if (occCount < best) {
                        best = occCount;
                        finder.posChoice = i;
                    }
                }
            }

            hashUpper(indexShape(index), begin(needle(pattern)) + finder.posChoice + finder.currentStep, length(needle(pattern)) - finder.posChoice - finder.currentStep);
            upperHash = value(indexShape(index));
            hash(indexShape(index), begin(needle(pattern)) + finder.posChoice + finder.currentStep, length(needle(pattern)) - finder.posChoice - finder.currentStep);

            //hashUpper(indexShape(index), begin(needle(pattern)) + finder.currentStep, length(needle(pattern)));
            //TShapeValue upperHash = value(indexShape(index));
            //hash(indexShape(index), begin(needle(pattern)) + finder.currentStep, length(needle(pattern)));
            //std::cerr << "count = " << countOccurrences(index, indexShape(index), upperHash) << std::endl;

            //finder.occurrences = getOccurrences(index, indexShape(index), upperHash);
            finder.occurrencesIter.init(index, indexShape(index), upperHash);
        } else {
            //++finder.occurrencesIter;
            finder.occurrencesIter.moveForward(needle(pattern), finder.posChoice + finder.currentStep);
        }
        if (atEnd(finder.occurrencesIter)) continue;

        // The current position is too far to the left in the text (so the needle would have to start negative position).
        if (*(finder.occurrencesIter) < finder.currentStep + finder.posChoice) continue;
        
        //std::cerr << "# " << infix(indexText(index), beginPosition(finder), length(needle(pattern))) << std::endl;
        
        // If the caller told us not to verify we simply return the candidate match (which might be a false positive).
        if (! TVerification::VALUE) return true;
        // In this case we also don't have to verify.
        if (length(needle(pattern)) <= Q) return true;
        // Verify
        if (needle(pattern) == infixWithLength(indexText(index), beginPosition(finder), length(needle(pattern)))) return true;
    }
}

// This function is to allow calls with a String instead of a Pattern object
template <typename TText, typename TIndexSpec, typename TVerification, typename TSpecFinder, typename TNeedle>
inline bool
find(Finder<Index<TText, TIndexSpec>, FinderQGramExtended<TVerification, TSpecFinder> > &finder, TNeedle const &ndl) {
    Pattern<TNeedle> ptn(ndl);
    return find(finder, ptn);
}

// ----------------------------------------------------------------------------
// beginPosition
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TVerification, typename TSpecFinder>
inline typename Position<Finder<Index<TText, TIndexSpec>, FinderQGramExtended<TVerification, TSpecFinder> > >::Type
beginPosition(Finder<Index<TText, TIndexSpec>, FinderQGramExtended<TVerification, TSpecFinder> > & me) {
    return *(me.occurrencesIter) - me.currentStep - me.posChoice;
}

// ----------------------------------------------------------------------------
// setHost
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TVerification, typename TSpecFinder>
inline void
setHost(Finder<Index<TText, TIndexSpec>, FinderQGramExtended<TVerification, TSpecFinder> > & finder, typename Parameter_<Index<TText, TIndexSpec> >::Type & index) {
    //typedef Index<TText, TIndexSpec>  TIndex;

    clear(finder);
    setValue(finder.index, index);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INDEX_FIND_INDEX_QGRAM_EXT_H_
