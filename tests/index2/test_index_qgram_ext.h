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
// Tests for the SeqAn module index2.
// ==========================================================================

#ifndef TESTS_INDEX_TEST_INDEX_QGRAM_EXT_H_
#define TESTS_INDEX_TEST_INDEX_QGRAM_EXT_H_

#include <algorithm>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/index2.h>

namespace seqan {

SEQAN_DEFINE_TEST(test_index_qgram_construct_find)
{
    static const unsigned int Q = 2u;
    typedef String<char>                                    TText;
    typedef String<char>                                    TNeedle;
    typedef Index<TText, IndexQGram<UngappedShape<Q> > >    TIndex;
    //typedef typename Fibre<TIndex, FibreShape>::Type        TShape;
    typedef String<unsigned int>                            TOccurrences;
    typedef Finder<TIndex, FinderQGramExtended<> >          TFinder;
    typedef Pattern<TNeedle>                                TPattern;
    
    TText text("issiZmataharissiXmatahariYissi");
    TIndex index(text);

    // getOccurrences
    {
        TNeedle ndl("is");
        hash(indexShape(index), begin(ndl));
        SEQAN_ASSERT_EQ(countOccurrences(index, indexShape(index)), 3u);
        unsigned int occ = getOccurrence(index, indexShape(index));
        SEQAN_ASSERT(occ == 0u || occ == 12u || occ == 26u);

        TOccurrences pos = getOccurrences(index, indexShape(index));
        // Because index structures can return the positions in arbitrary order (not necessarily from left to right)
        std::sort(begin(pos), end(pos));
        SEQAN_ASSERT_EQ(length(pos), 3u);
        SEQAN_ASSERT_EQ(pos[0], 0u);
        SEQAN_ASSERT_EQ(pos[1], 12u);
        SEQAN_ASSERT_EQ(pos[2], 26u);
    }

    // OccurrencesIterator
    {
        TNeedle ndl("is");
        hash(indexShape(index), begin(ndl));
        OccurrencesIterator<TIndex> occIt;
        occIt.init(index, indexShape(index), value(indexShape(index)) + 1);
        TOccurrences pos;
        //for (; !atEnd(occIt); ++occIt) {
        for (; !atEnd(occIt); occIt.moveForward(ndl, 0u)) {
            appendValue(pos, *occIt);
        }

        // Because index structures can return the positions in arbitrary order (not necessarily from left to right)
        std::sort(begin(pos), end(pos));
        SEQAN_ASSERT_EQ(length(pos), 3u);
        SEQAN_ASSERT_EQ(pos[0], 0u);
        SEQAN_ASSERT_EQ(pos[1], 12u);
        SEQAN_ASSERT_EQ(pos[2], 26u);
    }

    // Finder (needle is longer than Q)
    {
        TFinder finder(index);
        TPattern pattern("ssi");
        TOccurrences pos;
        while (find(finder, pattern)) {
            appendValue(pos, position(finder));
        }
        std::sort(begin(pos), end(pos));
        SEQAN_ASSERT_EQ(length(pos), 3u);
        SEQAN_ASSERT_EQ(pos[0], 1u);
        SEQAN_ASSERT_EQ(pos[1], 13u);
        SEQAN_ASSERT_EQ(pos[2], 27u);
    }

    // Finder (needle is shorter than Q)
    {
        TFinder finder(index);
        TPattern pattern("s");
        TOccurrences pos;
        while (find(finder, pattern)) {
            appendValue(pos, position(finder));
        }
        std::sort(begin(pos), end(pos));
        SEQAN_ASSERT_EQ(length(pos), 6u);
        SEQAN_ASSERT_EQ(pos[0], 1u);
        SEQAN_ASSERT_EQ(pos[1], 2u);
        SEQAN_ASSERT_EQ(pos[2], 13u);
        SEQAN_ASSERT_EQ(pos[3], 14u);
        SEQAN_ASSERT_EQ(pos[4], 27u);
        SEQAN_ASSERT_EQ(pos[5], 28u);
    }
}

SEQAN_DEFINE_TEST(test_index_qsample_construct_find)
{
    static const unsigned int Q = 3u;
    typedef String<char>                                    TText;
    typedef String<char>                                    TNeedle;
    typedef Index<TText, IndexQGram<UngappedShape<Q> > >    TIndex;
    //typedef typename Fibre<TIndex, FibreShape>::Type        TShape;
    typedef String<unsigned int>                            TOccurrences;
    typedef Finder<TIndex, FinderQGramExtended<> >          TFinder;
    typedef Pattern<TNeedle>                                TPattern;
    
    TText text("issiZmataharissiXmatahariYissi");
    TIndex index(text);
    setStepSize(index, 2u);

    // Finder (needle is longer than Q)
    {
        TFinder finder(index);
        TPattern pattern("matah");
        TOccurrences pos;
        while (find(finder, pattern)) {
            appendValue(pos, position(finder));
        }
        std::sort(begin(pos), end(pos));
        SEQAN_ASSERT_EQ(length(pos), 2u);
        SEQAN_ASSERT_EQ(pos[0], 5u);
        SEQAN_ASSERT_EQ(pos[1], 17u);
    }
}

SEQAN_DEFINE_TEST(test_index_qgram2l_construct_find)
{
    static const unsigned int Q = 2u;
    typedef String<char>                                    TText;
    typedef String<char>                                    TNeedle;
    typedef Index<TText, IndexQGram2L<UngappedShape<Q> > >  TIndex;
    //typedef typename Fibre<TIndex, FibreShape>::Type        TShape;
    typedef String<unsigned int>                            TOccurrences;
    typedef Finder<TIndex>                                  TFinder;
    typedef Pattern<TNeedle>                                TPattern;
    
    TText text("issiZmataharissiXmatahariYissi");
    TIndex index(text);
    //setSubsequenceLengthMin(index, 3u);
    //setSubsequenceLengthMax(index, 10u);
    setSubsequenceLength(index, 4u);
    
    // Prevent deletion of temporary sequences
    index.clearRedundantFibres = false;

    // Construct
    SEQAN_ASSERT(indexRequire(index, QGram2LBackEndIndex()));
    SEQAN_ASSERT(indexRequire(index, QGram2LFrontEndIndex()));
    
    // It is not possible to predict how exactly the internal data structures look,
    // because the algorithm is free to choose the subsequences.
    // Therefore we check for some general criteria which always have to be hold.

    // No data structure is empty
    SEQAN_ASSERT_GEQ(length(index.sequences), 1u);
    SEQAN_ASSERT_GEQ(length(indexDir(index.frontEndIndex)), 1u);
    SEQAN_ASSERT_GEQ(length(indexSA(index.frontEndIndex)), 1u);
    SEQAN_ASSERT_GEQ(length(indexDir(index.backEndIndex)), 1u);
    SEQAN_ASSERT_GEQ(length(indexSA(index.backEndIndex)), 1u);

    // There is one sentinal entry in the backend index
    SEQAN_ASSERT_EQ(length(index.sequences) + 1u, length(indexDir(index.backEndIndex)));

    // Test getter
    //SEQAN_ASSERT_EQ(getSubsequenceLengthMin(index), 3u);
    //SEQAN_ASSERT_EQ(getSubsequenceLengthMax(index), 10u);
    SEQAN_ASSERT_EQ(getSubsequenceLength(index), 4u);
    SEQAN_ASSERT_EQ(getStepSize(index), 1u);

    // getOccurrences
    {
        TNeedle ndl("is");
        hash(indexShape(index), begin(ndl));

        //SEQAN_ASSERT_EQ(countOccurrences(index, indexShape(index)), 3u);
        //unsigned int occ = getOccurrence(index, indexShape(index));
        //SEQAN_ASSERT(occ == 0u || occ == 12u || occ == 26u);

        TOccurrences pos = getOccurrences(index, indexShape(index));
        // Because index structures can return the positions in arbitrary order (not necessarily from left to right)
        std::sort(begin(pos), end(pos));
        SEQAN_ASSERT_EQ(length(pos), 3u);
        SEQAN_ASSERT_EQ(pos[0], 0u);
        SEQAN_ASSERT_EQ(pos[1], 12u);
        SEQAN_ASSERT_EQ(pos[2], 26u);
    }

    // OccurrencesIterator
    {
        TNeedle ndl("is");
        hash(indexShape(index), begin(ndl));
        OccurrencesIterator<TIndex> occIt;
        occIt.init(index, indexShape(index), value(indexShape(index)) + 1);
        TOccurrences pos;
        //for (; !atEnd(occIt); ++occIt) {
        for (; !atEnd(occIt); occIt.moveForward(ndl, 0u)) {
            appendValue(pos, *occIt);
        }

        // Because index structures can return the positions in arbitrary order (not necessarily from left to right)
        std::sort(begin(pos), end(pos));
        SEQAN_ASSERT_EQ(length(pos), 3u);
        SEQAN_ASSERT_EQ(pos[0], 0u);
        SEQAN_ASSERT_EQ(pos[1], 12u);
        SEQAN_ASSERT_EQ(pos[2], 26u);
    }

    // Finder (needle is longer than Q)
    {
        TFinder finder(index);
        TPattern pattern("ssi");
        TOccurrences pos;
        while (find(finder, pattern)) {
            appendValue(pos, position(finder));
        }
        std::sort(begin(pos), end(pos));
        SEQAN_ASSERT_EQ(length(pos), 3u);
        SEQAN_ASSERT_EQ(pos[0], 1u);
        SEQAN_ASSERT_EQ(pos[1], 13u);
        SEQAN_ASSERT_EQ(pos[2], 27u);
    }

    // Finder (needle is shorter than Q)
    {
        TFinder finder(index);
        TPattern pattern("s");
        TOccurrences pos;
        while (find(finder, pattern)) {
            appendValue(pos, position(finder));
        }
        std::sort(begin(pos), end(pos));
        SEQAN_ASSERT_EQ(length(pos), 6u);
        SEQAN_ASSERT_EQ(pos[0], 1u);
        SEQAN_ASSERT_EQ(pos[1], 2u);
        SEQAN_ASSERT_EQ(pos[2], 13u);
        SEQAN_ASSERT_EQ(pos[3], 14u);
        SEQAN_ASSERT_EQ(pos[4], 27u);
        SEQAN_ASSERT_EQ(pos[5], 28u);
    }
}

}  // namespace seqan

#endif  // TESTS_INDEX_TEST_INDEX_QGRAM_EXT_H_
