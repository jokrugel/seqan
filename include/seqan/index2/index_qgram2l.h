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
// A q-gram index structure using two levels to reduce the space consumption.
//
// [KWLL2005] M.-S. Kim and K.-Y. Whang and J.-G. Lee and M.-J. Lee
//   "n-Gram/2L: A Space and Time Efficient Two-Level n-Gram Inverted Index
//   Structure", VLDB 2005, http://portal.acm.org/citation.cfm?id=1083592.1083632
// 
// Some parts of the n-Gram/2L index are based on an experimental implementation
// by Johannes Merkle.
// ==========================================================================

#ifndef SEQAN_INDEX_INDEX_QGRAM2L_H_
#define SEQAN_INDEX_INDEX_QGRAM2L_H_

#include <functional> // binary_function
#include <iostream>
#include <set>
#include <vector>

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/sequence.h>

#include <seqan/find_index_approx/find_index_approx_base.h>  // For GetOccurrences<...>::Type

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Fibres
// ----------------------------------------------------------------------------

struct FibreFrontEndIndex_;
struct FibreBackEndIndex_;
typedef Tag<FibreFrontEndIndex_> const FibreFrontEndIndex;
typedef Tag<FibreBackEndIndex_> const FibreBackEndIndex;
typedef FibreText           QGram2LText;
typedef FibreRawText        QGram2LRawText;
typedef FibreShape          QGram2LShape;
typedef FibreFrontEndIndex  QGram2LFrontEndIndex;
typedef FibreBackEndIndex   QGram2LBackEndIndex;

// ----------------------------------------------------------------------------
// Construction algorithms
// ----------------------------------------------------------------------------
struct QGram2LFixed{};
struct QGram2LVariable{};

// ----------------------------------------------------------------------------
// Index
// ----------------------------------------------------------------------------

template <typename TShapeSpec, typename TAddressingSpec = void, typename TSuffixTreeSpec = IndexEsa<> >
struct IndexQGram2L{};

// TODO(krugel) StringSet has to be StringSet<TText, Dependent>, StringSet<TInfix, Dependent>, or StringSet<TInfix, Owner> ?    

template <typename TObject, unsigned int Q, typename TAddressingSpec, typename TSuffixTreeSpec>
class Index<TObject, IndexQGram2L<UngappedShape<Q>, TAddressingSpec, TSuffixTreeSpec> > {
public:
    typedef typename Fibre<Index, QGram2LText>::Type        TText;
    typedef          Holder<TText>                          THolder;
    typedef typename Fibre<Index, QGram2LShape>::Type       TShape;
    typedef typename Cargo<Index>::Type                     TCargo;
    typedef typename Size<Index>::Type                      TSize;
    typedef typename Position<Index>::Type                  TPosition;

    typedef typename Fibre<Index, QGram2LFrontEndIndex>::Type TFrontEndIndex;
    typedef typename Fibre<Index, QGram2LBackEndIndex>::Type TBackEndIndex;
    typedef typename Infix<const TText>::Type               TInfix;
    typedef StringSet<TInfix>                               TStringSet;
    typedef typename Size<TStringSet>::Type                 TStringSetSize;

    // TO-DO(krugel) Define better default values
    static const unsigned int defaultSubsequenceLengthMin = 2 * Q;
    static const unsigned int defaultSubsequenceLengthMax = 2 * Q;

    THolder         text;    // underlying text
    TShape          shape;   // underlying shape
    TCargo          cargo;   // user-defined cargo
    bool            clearRedundantFibres;
    TSize           subsequenceLengthMin;
    TSize           subsequenceLengthMax;
    TStringSet      sequences;
    TFrontEndIndex  frontEndIndex;
    TBackEndIndex   backEndIndex;

    // frontDir: q-gram hash --> index in SA of sequences
    // frontSA:  --> numbers in sequences (with offset)
    // backDir:  [m-gram hash / number in sequences] --> index in SA of text
    // backSA    --> positions of m-sequence in the text
    
    // frontDir: not sorted
    // frontSA:  sorted like frontDir

    Index():
        text(),
        clearRedundantFibres(true),
        subsequenceLengthMin(defaultSubsequenceLengthMin),
        subsequenceLengthMax(defaultSubsequenceLengthMax),
        frontEndIndex(sequences),
        backEndIndex(text) {}

    Index(Index &other):
        text(other.text),
        shape(other.shape),
        cargo(other.cargo),
        clearRedundantFibres(other.clearRedundantFibres),
        subsequenceLengthMin(other.subsequenceLengthMin),
        subsequenceLengthMax(other.subsequenceLengthMax),
        sequences(other.sequences),
        frontEndIndex(other.frontEndIndex),
        backEndIndex(other.backEndIndex) {}

    Index(Index const &other):
        text(other.text),
        shape(other.shape),
        cargo(other.cargo),
        clearRedundantFibres(other.clearRedundantFibres),
        subsequenceLengthMin(other.subsequenceLengthMin),
        subsequenceLengthMax(other.subsequenceLengthMax),
        sequences(other.sequences),
        frontEndIndex(other.frontEndIndex),
        backEndIndex(other.backEndIndex) {}

    Index(TSize _subsequenceLengthMin, TSize _subsequenceLengthMax):
        text(),
        clearRedundantFibres(true),
        subsequenceLengthMin(_subsequenceLengthMin),
        subsequenceLengthMax(_subsequenceLengthMax),
        frontEndIndex(sequences),
        backEndIndex(text) {}

    template <typename TText_>
    Index(TText_ &_text):
        text(_text),
        clearRedundantFibres(true),
        subsequenceLengthMin(defaultSubsequenceLengthMin),
        subsequenceLengthMax(defaultSubsequenceLengthMax),
        frontEndIndex(sequences),
        backEndIndex(text) {}

    // For fixed subsequence length
    template <typename TText_>
    Index(TText_ &_text, TSize _subsequenceLength):
        text(_text),
        clearRedundantFibres(true),
        subsequenceLengthMin(_subsequenceLength),
        subsequenceLengthMax(_subsequenceLength),
        frontEndIndex(sequences),
        backEndIndex(text) {}

    template <typename TText_>
    Index(TText_ &_text, TSize _subsequenceLengthMin, TSize _subsequenceLengthMax):
        text(_text),
        clearRedundantFibres(true),
        subsequenceLengthMin(_subsequenceLengthMin),
        subsequenceLengthMax(_subsequenceLengthMax),
        frontEndIndex(sequences),
        backEndIndex(text) {}

    template <typename TText_>
    Index(TText_ const &_text, TSize _subsequenceLengthMin, TSize _subsequenceLengthMax):
        text(_text),
        clearRedundantFibres(true),
        subsequenceLengthMin(_subsequenceLengthMin),
        subsequenceLengthMax(_subsequenceLengthMax),
        frontEndIndex(sequences),
        backEndIndex(text) {}

};

// ----------------------------------------------------------------------------
// Compare functor for suffix tree nodes (based on the space savings)
// ----------------------------------------------------------------------------

template <typename TSize, typename TNodes, typename TOccCounter>
class CompareNodeQGram2L : public std::binary_function<typename Position<TNodes>::Type, typename Position<TNodes>::Type, bool> {
private:
    typedef typename Position<TNodes>::Type TPosition;
    
    TNodes & nodes;
    TOccCounter & occCounter;
    const TSize OVERLAP;

public:
    CompareNodeQGram2L(TNodes & _nodes, TOccCounter & _occCounter, TSize _overlap)
        : nodes(_nodes), occCounter(_occCounter), OVERLAP(_overlap) {
    }
    bool operator()(const TPosition & lhs, const TPosition & rhs) const {
        //return repLength(nodes[lhs]) > repLength(nodes[rhs]);
        return _spaceSavings(nodes, occCounter, lhs, OVERLAP) > _spaceSavings(nodes, occCounter, rhs, OVERLAP);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
struct Fibre<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> >, QGram2LFrontEndIndex> {
    //typedef Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > TIndex;
    //typedef typename TIndex::TSize                          TSize;
    //typedef typename TIndex::TPosition                      TPosition;
    //typedef typename TIndex::TInfix                         TInfix1;
    //typedef typename TIndex::TStringSet                     TStringSet2;

    typedef typename Infix<const TObject>::Type             TInfix;
    typedef StringSet<TInfix>                               TStringSet;
    typedef Index<TStringSet, IndexQGram<TShapeSpec, TAddressingSpec> > Type;
};

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
struct Fibre<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> >, QGram2LBackEndIndex> {
    typedef Index<TObject, IndexQGram<SimpleShape, OpenAddressing> > Type;
};

// Use the index value type as shape value type
template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
struct Fibre<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> >, FibreShape> {
    typedef Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > TIndex;
    typedef Shape<typename Value<TIndex>::Type, TShapeSpec> Type;
};

// Allow different value types for the shape
template <typename TObject, typename TShapeValue, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
struct Fibre<Index<TObject, IndexQGram2L<Shape<TShapeValue, TShapeSpec>, TAddressingSpec, TSuffixTreeSpec> >, FibreShape> {
    typedef Shape<TShapeValue, TShapeSpec>                  Type;
};

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
struct GetOccurrences<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > > {
    typedef Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > TIndex;
    typedef String<typename SAValue<TIndex>::Type>          Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Private functions for index creation
// ----------------------------------------------------------------------------

// Go to the parent node of a suffix tree (using the additional data structures)
template <typename TIndex, typename TSpec, typename TNodes, typename TMapping>
inline bool
_goUp(Iter<TIndex, VSTree<TopDown<TSpec> > > &it, TNodes & nodes, TMapping & mapping) {
    if (isRoot(it)) return false;
    //std::cout << "^" << representative(it) << " ~ " << _getId(nodeUp(it)) << " ~ " << mapping[_getId(nodeUp(it))] << std::endl;
    it = nodes[mapping[_getId(nodeUp(it))]];
    //std::cout << "'" << representative(it) << std::endl;
    return true;
}

// Go to next node of a suffix tree (using the additional data structures)
// Copied and modified from /core/include/seqan/index/index_esa_stree.h
template <typename TIndex, typename TSpec, typename TSpec2, typename TNodes, typename TMapping>
inline void
_goNextInSubtree(Iter<TIndex, VSTree<TopDown<TSpec> > > &it, const Iter<TIndex, VSTree<TSpec2> > & rootIt, TNodes & nodes, TMapping & mapping) {
    // preorder dfs
    do {
        if (!goDown(it) && !goRight(it))
            while (_goUp(it, nodes, mapping) && _getId(value(it)) != _getId(value(rootIt)) && !goRight(it) && _getId(value(it)) != _getId(value(rootIt))) {}
        if (_getId(value(it)) == _getId(value(rootIt))) {
            clear(it);
            return;
        }
    } while (!nodePredicate(it)); // always true by default
}

// Follow the suffix link of a node of a suffix tree (using the additional data structures)
template <typename TIndex, typename TSpec, typename TNodes, typename TLeaves>
inline bool
_goSuffixLink(Iter<TIndex, VSTree<TopDown<TSpec> > > & it, TNodes & nodes, TLeaves & leaves) {
    if (isRoot(it)) return false;
    it = nodes[leaves[getOccurrence(it) + 1u]];
    return true;
}

// Calculate the space savings that result, if the representative of the node at position key
// is added to the index structure (can be positive or negative).
template <typename TTreeIter, typename TOccCounter, typename TSize>
inline typename MakeSigned<TSize>::Type
_spaceSavings(String<TTreeIter> & nodes, TOccCounter & occCounter, typename Position<String<TTreeIter> >::Type key, const TSize OVERLAP) {
    typedef typename MakeSigned<TSize>::Type TSizeSigned;
    
    // TO-DO(krugel) Better estimation of space savings?
    TTreeIter & it = nodes[key];
    TSize traditionalSpace = occCounter[key];
    TSize improvedSpace = repLength(it) + countOccurrences(it) + 1u + 2u * OVERLAP * countOccurrences(it);
    //                           10            + 42                   + 1 + 2 * 1 * 42
    //                           22            + 19                   + 1 + 2 * 2 * 19  = 118
    // _spaceSavings = 344 - 118
    return (TSizeSigned) traditionalSpace - (TSizeSigned) improvedSpace;
}

// ----------------------------------------------------------------------------
// indexCreate
// ----------------------------------------------------------------------------

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline bool
indexCreate(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, QGram2LFrontEndIndex, Default const) {
    if (! indexRequire(index, QGram2LBackEndIndex())) return false;

    //std::cout << "text length = " << length(value(index.text)) << std::endl << std::endl;
    //std::cout << "----------------------------------------------------------------------" << std::endl;
    //std::cout << "Resulting sequences..." << std::endl << std::endl;
    //std::cout << length(index.sequences) << std::endl << std::endl;

    if (! indexRequire(index.frontEndIndex, QGramSADir())) return false;

    //std::cout << "----------------------------------------------------------------------" << std::endl;
    //std::cout << "Resulting sequences..." << std::endl << std::endl;
    //std::cout << length(index.sequences) << std::endl << std::endl;
    //for (TStringSetIter sequenceIt(index.sequences); !atEnd(sequenceIt); goNext(sequenceIt)) {
    //    std::cout << "- " << *sequenceIt << std::endl;
    //}

    // Destroy sequences, we don't need them anymore
    if (index.clearRedundantFibres) {
        clear(index.sequences);
        //shrinkToFit(index.sequences);
    }

    return true;
}

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline bool
indexCreate(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, QGram2LBackEndIndex, Default const) {
    if (getSubsequenceLengthMin(index) == getSubsequenceLengthMax(index)) {
        return indexCreate(index, QGram2LBackEndIndex(), QGram2LFixed());
    } else {
        SEQAN_ASSERT_MSG(false, "QGram2LVariable is not supported any more");
        return false; //indexCreate(index, QGram2LBackEndIndex(), QGram2LVariable());
    }
}

// Creation with fixed subsequence length (corresponds to the standard two-level index).
template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline bool
indexCreate(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, QGram2LBackEndIndex, QGram2LFixed const) {
    typedef Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> >  TIndex;
    typedef typename TIndex::TSize                          TSize;
    typedef typename Infix<TText>::Type                     TInfix;
    //typedef typename Iterator<TText>::Type                  TTextIter;
    //typedef typename Position<TText>::Type                  TTextPosition;
    //typedef typename TIndex::TStringSet                     TStringSet;
    //typedef typename Iterator<TStringSet>::Type             TStringSetIter;
    typedef typename TIndex::TBackEndIndex                  TBackEndIndex;
    typedef typename Fibre<TBackEndIndex, QGramDir>::Type   TBackEndIndexDir; // String<Size>
    typedef typename Fibre<TBackEndIndex, QGramSA>::Type    TBackEndIndexSA;  // String<SAValue>
    typedef typename Fibre<TBackEndIndex, QGramBucketMap>::Type TBackEndIndexBucketMap;
    typedef typename Iterator<TBackEndIndexSA>::Type        TBackEndIndexSAIter;
    typedef typename Value<TText>::Type                     TValue;
    typedef Shape<TValue, SimpleShape>                      TShapeBackEnd;

    SEQAN_ASSERT_EQ_MSG(getSubsequenceLengthMin(index), getSubsequenceLengthMax(index),
        "For using this construction algorithm, subsequenceLengthMin and subsequenceLengthMax have to be equal.");
    TSize subsequenceLength = getSubsequenceLengthMin(index);
    
    // Some shortcuts
    TText & text = indexText(index);
    TBackEndIndex & backEndIndex = index.backEndIndex;
    TShapeBackEnd & shape = indexShape(backEndIndex);
    TBackEndIndexDir & dir = indexDir(backEndIndex);
    TBackEndIndexDir & sa = indexSA(backEndIndex);

    TBackEndIndexDir dirNew;

    const TSize ALREADY_ADDED = MaxValue<TSize>::VALUE;
    TSize Q = length(indexShape(index));
    indexShape(backEndIndex) = TShapeBackEnd(subsequenceLength);
    TSize overlap = Q - 1u;
    TSize stepSize = subsequenceLength - overlap;

    //std::cerr << "Q = " << Q << std::endl;
    //std::cerr << "subsequenceLength = " << subsequenceLength << std::endl;
    //std::cerr << "overlap = " << overlap << std::endl;
    //std::cerr << "stepSize = " << stepSize << std::endl;

    setStepSize(backEndIndex, stepSize);

    // TODO(krugel) Avoid building the backEndIndex directory like this and try to build it manually
    if (! indexRequire(backEndIndex, QGramSADir())) return false;

    // TODO(krugel) Iterate instead over dir?
    for (TBackEndIndexSAIter it = begin(sa); !atEnd(it, sa); ++it) {
        TInfix sequence = infixWithLength(text, value(it), subsequenceLength);

        hash(shape, begin(sequence));
        // std::cerr << value(shape) << ", ";
        TSize dirPos = getBucket(indexBucketMap(backEndIndex), value(shape));

        //std::cerr << sequence << std::endl;
        if (dir[dirPos] != ALREADY_ADDED) { // This sequence has not been added before
            appendValue(index.sequences, sequence);
            appendValue(dirNew, dir[dirPos]);
            dir[dirPos] = ALREADY_ADDED;
        }
    }
    appendValue(dirNew, length(sa));

    TSize lastSequenceLength = length(text) - length(sa) * stepSize;
    if (lastSequenceLength >= Q) {
        TInfix sequence = infix(text, length(text) - lastSequenceLength, endPosition(text));
        //std::cerr << sequence << std::endl;
        appendValue(sa, length(text) - lastSequenceLength);
        appendValue(dirNew, length(sa));
        appendValue(index.sequences, sequence);
    }
    
    // Copy the just created new directory into the index
    // TODO(krugel) first clear(dir)?
    dir = dirNew;
    
    // Clear doesn't work so we overwrite with empty
    //clear(indexBucketMap(backEndIndex));
    //clear(indexBucketMap(backEndIndex).qgramCode);
    // shrinkToFit?
    indexBucketMap(backEndIndex) = TBackEndIndexBucketMap();
    
    return true;
}

// ----------------------------------------------------------------------------
// getOccurrences
// ----------------------------------------------------------------------------

// Find a limited number of occurrences (especially useful, if you only want to find one occurrence).
template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TShapeSpec2, typename TValue, typename TSize>
inline void
_getOccurrencesWithLimit(
        Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &index,
        Shape<TValue, TShapeSpec2> const &shape,
        typename Value<Shape<TValue, TShapeSpec2> >::Type upperHash,
        TSize limit,
        typename GetOccurrences<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type & results
        ) {
    typedef Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > TIndex;
    typedef typename Position<TObject>::Type                TPosition;
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

    TFrontEndIndexDirPosition hashValue = value(shape);
    hashValue = getBucket(indexBucketMap(index.frontEndIndex), hashValue);
    upperHash = getBucket(indexBucketMap(index.frontEndIndex), upperHash);
    
    TFrontEndIndexSAPosition saFrontPosEnd = dirAt(upperHash, index.frontEndIndex);

    //std::cerr << __LINE__ << ": " << hashValue << ".." << upperHash << std::endl;
    //std::cerr << __LINE__ << ": " << dirAt(hashValue, index.frontEndIndex) << ".." << saFrontPosEnd << std::endl;

    for (TFrontEndIndexSAPosition saFrontPos = dirAt(hashValue, index.frontEndIndex); saFrontPos < saFrontPosEnd; ++saFrontPos) {
        TFrontEndIndexSAValue saValue = saAt(saFrontPos, index.frontEndIndex);
        TBackEndIndexDirPosition seqNo = getSeqNo(saValue);
        //std::cerr << seqNo << " * " << getSeqOffset(saValue) << " (" << dirAt(seqNo, index.backEndIndex) << ".." << dirAt(seqNo + 1u, index.backEndIndex) << ")" << std::endl;
        for (TBackEndIndexSAPosition saBackPos = dirAt(seqNo, index.backEndIndex); saBackPos < dirAt(seqNo + 1u, index.backEndIndex); ++saBackPos) {
            TPosition seqPos = saAt(saBackPos, index.backEndIndex);
            //std::cerr << "- textPos = " << seqPos + getSeqOffset(saValue) << " \t(in V-sequence " << seqNo << " at position " << getSeqOffset(saValue) << ")" << std::endl;
            appendValue(results, seqPos + getSeqOffset(saValue));
            if (limit != 0u && length(results) >= limit) return;
        }
    }
}

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TShapeSpec2, typename TValue>
inline typename SAValue<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type
getOccurrence(Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &index, Shape<TValue, TShapeSpec2> const &shape) {
    typedef typename GetOccurrences<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type TOccurrences;
    TOccurrences occ;
    _getOccurrencesWithLimit(index, shape, value(shape) + 1u, 1u, occ);
    return empty(occ) ? 0u : occ[0u];
    // TODO(krugel) Better throw an exception if not found? or use MaxValue
}

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TShapeSpec2, typename TValue>
inline typename SAValue<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type
getOccurrence(Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, Shape<TValue, TShapeSpec2 > const &shape) {
    indexRequire(index, QGram2LBackEndIndex());
    indexRequire(index, QGram2LFrontEndIndex());
    return getOccurrence(const_cast<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &>(index), shape);
}

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TShapeSpec2, typename TValue>
inline typename GetOccurrences<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type
getOccurrences(Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &index, Shape<TValue, TShapeSpec2> const &shape) {
    typedef typename GetOccurrences<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type TOccurrences;
    TOccurrences occ;
    _getOccurrencesWithLimit(index, shape, value(shape) + 1u, 0u, occ);
    return occ;
}

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TShapeSpec2, typename TValue>
inline typename GetOccurrences<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type
getOccurrences(Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, Shape<TValue, TShapeSpec2> const &shape) {
    indexRequire(index, QGram2LBackEndIndex());
    indexRequire(index, QGram2LFrontEndIndex());
    return getOccurrences(const_cast<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &>(index), shape);
}

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TShapeSpec2, typename TValue>
inline typename GetOccurrences<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type
getOccurrences(Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &index, Shape<TValue, TShapeSpec2> const &shape, typename Value<Shape<TValue, TShapeSpec2> >::Type upperHash) {
    typedef typename GetOccurrences<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type TOccurrences;
    TOccurrences occ;
    _getOccurrencesWithLimit(index, shape, upperHash, 0u, occ);
    return occ;
}

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TShapeSpec2, typename TValue>
inline typename GetOccurrences<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type
getOccurrences(Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, Shape<TValue, TShapeSpec2> const &shape, typename Value<Shape<TValue, TShapeSpec2> >::Type upperHash) {
    indexRequire(index, QGram2LBackEndIndex());
    indexRequire(index, QGram2LFrontEndIndex());
    return getOccurrences(const_cast<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &>(index), shape, upperHash);
}

// ----------------------------------------------------------------------------
// countOccurrences
// ----------------------------------------------------------------------------

//template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TShapeSpec2, typename TValue>
//inline typename Size<typename GetOccurrences<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type>::Type
//countOccurrences(Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &index, Shape<TValue, TShapeSpec2> const &shape) {
//    return length(getOccurrences(index, shape));
//}

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TShapeSpec2, typename TValue>
inline typename Size<typename GetOccurrences<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type>::Type
countOccurrences(Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &index, Shape<TValue, TShapeSpec2> const &shape, typename Value<Shape<TValue, TShapeSpec2> >::Type upperHash) {
    typedef Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > TIndex;
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
    typedef typename Size<typename GetOccurrences<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type>::Type TSize;

    TFrontEndIndexDirPosition hashValue = value(shape);
    hashValue = getBucket(indexBucketMap(index.frontEndIndex), hashValue);
    upperHash = getBucket(indexBucketMap(index.frontEndIndex), upperHash);
    
    TFrontEndIndexSAPosition saFrontPosEnd = dirAt(upperHash, index.frontEndIndex);
    TSize result = 0u;
    for (TFrontEndIndexSAPosition saFrontPos = dirAt(hashValue, index.frontEndIndex); saFrontPos < saFrontPosEnd; ++saFrontPos) {
        TFrontEndIndexSAValue saValue = saAt(saFrontPos, index.frontEndIndex);
        TBackEndIndexDirPosition seqNo = getSeqNo(saValue);
        result += dirAt(seqNo + 1u, index.backEndIndex) - dirAt(seqNo, index.backEndIndex);
    }
    return result;
}

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TShapeSpec2, typename TValue>
inline typename Size<typename GetOccurrences<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type>::Type
countOccurrences(Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &index, Shape<TValue, TShapeSpec2> const &shape) {
    return countOccurrences(index, shape, value(shape) + 1u);
}

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TShapeSpec2, typename TValue>
inline typename Size<typename GetOccurrences<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type>::Type 
countOccurrences(Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, Shape<TValue, TShapeSpec2> const &shape, typename Value<Shape<TValue, TShapeSpec2> >::Type upperHash) {
    indexRequire(index, QGram2LBackEndIndex());
    indexRequire(index, QGram2LFrontEndIndex());
    return countOccurrences(const_cast<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &>(index), shape, upperHash);
}

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TShapeSpec2, typename TValue>
inline typename Size<typename GetOccurrences<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > >::Type>::Type 
countOccurrences(Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, Shape<TValue, TShapeSpec2> const &shape) {
    indexRequire(index, QGram2LBackEndIndex());
    indexRequire(index, QGram2LFrontEndIndex());
    return countOccurrences(const_cast<Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &>(index), shape);
}

// ----------------------------------------------------------------------------
// get/set for stepSize and subsequenceLength
// ----------------------------------------------------------------------------

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline typename Size<TText>::Type 
getStepSize(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &index) {
    ignoreUnusedVariableWarning(index);
    return 1u;
}

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TSize>
inline void
setStepSize(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, TSize stepSize) {
    ignoreUnusedVariableWarning(index);
    SEQAN_ASSERT_EQ_MSG(stepSize, 1u, "For IndexQGram2L the stepSize has to be 1.");
}

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TSize>
inline void
setSubsequenceLength(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, TSize subsequenceLength) {
    index.subsequenceLengthMin = subsequenceLength;
    index.subsequenceLengthMax = subsequenceLength;
}

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TSize>
inline void
setSubsequenceLengthMin(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, TSize subsequenceLengthMin) {
    index.subsequenceLengthMin = subsequenceLengthMin;
}

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec, typename TSize>
inline void
setSubsequenceLengthMax(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, TSize subsequenceLengthMax) {
    index.subsequenceLengthMax = subsequenceLengthMax;
}

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline typename Size<TText>::Type
getSubsequenceLength(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &index) {
    return index.subsequenceLengthMin;
}

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline typename Size<TText>::Type
getSubsequenceLengthMin(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &index) {
    return index.subsequenceLengthMin;
}

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline typename Size<TText>::Type
getSubsequenceLengthMax(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &index) {
    return index.subsequenceLengthMax;
}

// ----------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline void
clear(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index) {
    clear(getFibre(index, QGram2LFrontEndIndex()));
    clear(getFibre(index, QGram2LBackEndIndex()));
    clear(index.sequences);
}

// ----------------------------------------------------------------------------
// empty
// ----------------------------------------------------------------------------

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline bool
empty(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index) {
    return empty(index.frontEndIndex);
}

// ----------------------------------------------------------------------------
// getFibre
// ----------------------------------------------------------------------------

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline typename Fibre<Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> >, QGram2LBackEndIndex>::Type & 
getFibre(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, QGram2LBackEndIndex) {
    return index.backEndIndex;
}

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline typename Fibre<Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const, QGram2LBackEndIndex>::Type & 
getFibre(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &index, QGram2LBackEndIndex) {
    return index.backEndIndex;
}

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline typename Fibre<Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> >, QGram2LFrontEndIndex>::Type & 
getFibre(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, QGram2LFrontEndIndex) {
    return index.frontEndIndex;
}

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline typename Fibre<Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const, QGram2LFrontEndIndex>::Type & 
getFibre(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > const &index, QGram2LFrontEndIndex) {
    return index.frontEndIndex;
}

// ----------------------------------------------------------------------------
// indexSupplied
// ----------------------------------------------------------------------------

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline bool
indexSupplied(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, QGram2LFrontEndIndex) {
    return ! empty(indexDir(getFibre(index, QGram2LFrontEndIndex())));
}

template <typename TText, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline bool
indexSupplied(Index<TText, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, QGram2LBackEndIndex) {
    return ! empty(indexDir(getFibre(index, QGram2LBackEndIndex())));
}

// ----------------------------------------------------------------------------
// open
// ----------------------------------------------------------------------------

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline bool
open(Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, const char *fileName, int openMode) {
    String<char> name;

    name = fileName; append(name, ".txt");
    if (! open(getFibre(index, QGram2LText()), toCString(name), openMode)) return false;

    name = fileName; append(name, ".front.dir");
    if (! open(getFibre(getFibre(index, QGram2LFrontEndIndex()), QGramDir()), toCString(name), openMode)) return false;

    name = fileName; append(name, ".front.sa");
    if (! open(getFibre(getFibre(index, QGram2LFrontEndIndex()), QGramSA()), toCString(name), openMode)) return false;

    name = fileName; append(name, ".back.dir");
    if (! open(getFibre(getFibre(index, QGram2LBackEndIndex()), QGramDir()), toCString(name), openMode)) return false;

    name = fileName; append(name, ".back.sa");
    if (! open(getFibre(getFibre(index, QGram2LBackEndIndex()), QGramSA()), toCString(name), openMode)) return false;

    return true;
}

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline bool
open(Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, const char *fileName) {
    return open(index, fileName, OPEN_RDONLY);
}

// ----------------------------------------------------------------------------
// save
// ----------------------------------------------------------------------------

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline bool
save(Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, const char *fileName, int openMode) {
    String<char> name;

    name = fileName; append(name, ".txt");
    if (! save(getFibre(index, QGram2LText()), toCString(name), openMode)) return false;

    name = fileName; append(name, ".front.dir");
    if (! save(getFibre(getFibre(index, QGram2LFrontEndIndex()), QGramDir()), toCString(name), openMode)) return false;

    name = fileName; append(name, ".front.sa");
    if (! save(getFibre(getFibre(index, QGram2LFrontEndIndex()), QGramSA()), toCString(name), openMode)) return false;

    name = fileName; append(name, ".back.dir");
    if (! save(getFibre(getFibre(index, QGram2LBackEndIndex()), QGramDir()), toCString(name), openMode)) return false;

    name = fileName; append(name, ".back.sa");
    if (! save(getFibre(getFibre(index, QGram2LBackEndIndex()), QGramSA()), toCString(name), openMode)) return false;

    return true;
}

template <typename TObject, typename TShapeSpec, typename TAddressingSpec, typename TSuffixTreeSpec>
inline bool
save(Index<TObject, IndexQGram2L<TShapeSpec, TAddressingSpec, TSuffixTreeSpec> > &index, const char *fileName) {
    return save(index, fileName, OPEN_WRONLY | OPEN_CREATE);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INDEX_INDEX_QGRAM2L_H_
