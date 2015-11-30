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

// File adapted from /core/tests/find/test_find.cpp

#ifndef TESTS_INDEX2_TEST_FIND_INDEX_H_
#define TESTS_INDEX2_TEST_FIND_INDEX_H_

#include <algorithm>
#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <seqan/index2.h>

namespace seqan {

// TODO(krugel) Don't overwrite DefaultFinder
// Here we want to use the extended q-gram finder for also longer patterns
template <typename TText, unsigned Q, typename TSpec>
struct DefaultFinder<Index<TText, IndexQGram<UngappedShape<Q>, TSpec> > > {
    typedef FinderQGramExtended<True>                       Type;
};

template <typename TIndex>
inline void
configureIndex(TIndex & index) {
    ignoreUnusedVariableWarning(index);
    // do nothing
}

// The 3-gram index is sampled with step length 2 (the 2-gram index is full with step size 1).
template <typename TText>
inline void
configureIndex(Index<TText, IndexQGram<UngappedShape<3> > > & index) {
    setStepSize(index, 2u);
}

template <typename TIndexSpec>
void testFindIndex() {
    typedef String<unsigned int>                            TPositions;
    
/*
    // Save/open works with all indexes exept for IndexSttd64

    //____________________________________________________________________________
    // Test0 - save/open

    {
        typedef String<char>                                THaystack;
        typedef String<char>                                TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;    
        typedef Finder<TIndex>                              TFinder;
        typedef typename DefaultIndexPattern<TNeedle, TIndex>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;

        CharString filename = "myindex";
        
        TPositions pos;
        THaystack hstck("Dies ist ein Haystack. Ja, das ist wirklich einer!");
        
        {
            TIndex index0(hstck);
            configureIndex(index0);
            configureSave(index0, toCString(filename));
            
            TFinder finder0(index0);
            TNeedle ndl0("ist");	
            TPattern pattern0(ndl0);
            find(finder0, pattern0);
            
            // Save the index to disk.
            save(index0, toCString(filename));
        }
        
        // Reopen the index and check if it still works.
        TIndex index;
        SEQAN_ASSERT_EQ(length(index), length(hstck));
        
        // TODO(krugel) SEQAN_ASSERT(open(...), save(...))
        open(index, toCString(filename));
        
        TFinder finder(index);
        TNeedle ndl("ist");	
        TPattern pattern(ndl);

        SEQAN_ASSERT_EQ(length(index), length(hstck));
        SEQAN_ASSERT_EQ(indexText(index), hstck);
        SEQAN_ASSERT_EQ(indexText(haystack(finder)), hstck);

        while (find(finder, pattern)) {
            //std::cerr << "Position = " << position(finder) << std::endl;
            appendValue(pos, position(finder));
            SEQAN_ASSERT_EQ(position(finder), beginPosition(finder));
            SEQAN_ASSERT_EQ(endPosition(finder), beginPosition(finder) + length(finder));
            SEQAN_ASSERT_EQ(length(finder), length(ndl));
            //SEQAN_ASSERT_EQ(begin(finder), begin(hstck) + beginPosition(finder)); // doesn't compile with IndexEsa and IndexWotd
            //SEQAN_ASSERT_EQ(end(finder), begin(hstck) + endPosition(finder));     // doesn't compile with IndexEsa and IndexWotd
            SEQAN_ASSERT_EQ(infix(finder), ndl);
        }
        
        SEQAN_ASSERT_EQ(needle(pattern), ndl);
        SEQAN_ASSERT_EQ(host(pattern), ndl);
        SEQAN_ASSERT_EQ(host(reinterpret_cast<Pattern<String<char> > const &>(pattern)), ndl);

        // Because index structures can return the positions in arbitrary order (not necessarily from left to right)
        if (! MatchesAreSorted<TFinder>::Type::VALUE) {
            std::sort(begin(pos), end(pos));
        }
        // q-sample cannot work correctly because needles are smaller than Q + stepSize - 1
        SEQAN_ASSERT_EQ(length(pos), 2u);
        SEQAN_ASSERT_EQ(pos[0], 5u);
        SEQAN_ASSERT_EQ(pos[1], 31u);
    }
return;
*/

    //____________________________________________________________________________
    // Test1 - small ndl

    {
        typedef String<char>                                THaystack;
        typedef String<char>                                TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;    
        typedef Finder<TIndex>                              TFinder;
        typedef typename DefaultIndexPattern<TNeedle, TIndex>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;

        TPositions pos;
        THaystack hstck("Dies ist ein Haystack. Ja, das ist wirklich einer!");
        TIndex index(hstck);
        configureIndex(index);
        TFinder finder(index);
        TNeedle ndl("ist");	
        TPattern pattern(ndl);
        
        SEQAN_ASSERT_EQ(length(index), length(hstck));
        SEQAN_ASSERT_EQ(indexText(index), hstck);
        SEQAN_ASSERT_EQ(indexText(haystack(finder)), hstck);

        while (find(finder, pattern)) {
            //std::cerr << "Position = " << position(finder) << std::endl;
            appendValue(pos, position(finder));
            SEQAN_ASSERT_EQ(position(finder), beginPosition(finder));
            SEQAN_ASSERT_EQ(endPosition(finder), beginPosition(finder) + length(finder));
            SEQAN_ASSERT_EQ(length(finder), length(ndl));
            //SEQAN_ASSERT_EQ(begin(finder), begin(hstck) + beginPosition(finder)); // doesn't compile with IndexEsa and IndexWotd
            //SEQAN_ASSERT_EQ(end(finder), begin(hstck) + endPosition(finder));     // doesn't compile with IndexEsa and IndexWotd
            SEQAN_ASSERT_EQ(infix(finder), ndl);
        }
        
        SEQAN_ASSERT_EQ(needle(pattern), ndl);
        SEQAN_ASSERT_EQ(host(pattern), ndl);
        SEQAN_ASSERT_EQ(host(reinterpret_cast<Pattern<String<char> > const &>(pattern)), ndl);

        // Because index structures can return the positions in arbitrary order (not necessarily from left to right)
        if (! MatchesAreSorted<TFinder>::Type::VALUE) {
            std::sort(begin(pos), end(pos));
        }
        // q-sample cannot work correctly because needles are smaller than Q + stepSize - 1
        SEQAN_ASSERT_EQ(length(pos), 2u);
        SEQAN_ASSERT_EQ(pos[0], 5u);
        SEQAN_ASSERT_EQ(pos[1], 31u);
    }

    //____________________________________________________________________________
    // Test2 - large ndl

    {
        typedef String<char>                                THaystack;
        typedef String<char>                                TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;    
        typedef Finder<TIndex>                              TFinder;
        typedef typename DefaultIndexPattern<TNeedle, TIndex>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;

        TPositions pos;
        THaystack hstck("abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef");
        TIndex index(hstck);
        configureIndex(index);
        TFinder finder(index);
        TNeedle ndl("abcdefghijklmnopqrstuvwxyzabcdefg");
        TPattern pattern(ndl);

        while (find(finder, pattern)) {
            append(pos, position(finder));
            SEQAN_ASSERT_EQ(position(finder), beginPosition(finder));
            SEQAN_ASSERT_EQ(endPosition(finder), beginPosition(finder) + length(finder));
            SEQAN_ASSERT_EQ(length(finder), length(ndl));
            SEQAN_ASSERT_EQ(infix(finder), ndl);
        }

        if (! MatchesAreSorted<TFinder>::Type::VALUE) {
            std::sort(begin(pos), end(pos));
        }
        SEQAN_ASSERT_EQ(length(pos), 2u);
        SEQAN_ASSERT_EQ(pos[0], 0u);
        SEQAN_ASSERT_EQ(pos[1], 26u);
    }

    //____________________________________________________________________________
    // Test3 - different alphabet, small ndl

    {
        typedef String<Dna>                                 THaystack;
        typedef String<Dna>                                 TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;    
        typedef Finder<TIndex>                              TFinder;
        typedef typename DefaultIndexPattern<TNeedle, TIndex>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;

        // TODO(krugel) Find a solution for IndexQGram at the end of the text
        TPositions pos;
        THaystack hstck("aaaaaaacaaTT");
        TIndex index(hstck);
        configureIndex(index);
        TFinder finder(index);
        TNeedle ndl("aa");
        TPattern pattern(ndl);

        while (find(finder, pattern)) {
            append(pos, position(finder));
            SEQAN_ASSERT_EQ(position(finder), beginPosition(finder));
            SEQAN_ASSERT_EQ(endPosition(finder), beginPosition(finder) + length(finder));
            SEQAN_ASSERT_EQ(length(finder), length(ndl));
            SEQAN_ASSERT_EQ(infix(finder), ndl);
        }

        if (! MatchesAreSorted<TFinder>::Type::VALUE) {
            std::sort(begin(pos), end(pos));
        }
        // q-sample cannot work correctly because needles are smaller than Q + stepSize - 1
        SEQAN_ASSERT_EQ(length(pos), 7u);
        SEQAN_ASSERT_EQ(pos[0], 0u);
        SEQAN_ASSERT_EQ(pos[1], 1u);
        SEQAN_ASSERT_EQ(pos[2], 2u);
        SEQAN_ASSERT_EQ(pos[3], 3u);
        SEQAN_ASSERT_EQ(pos[4], 4u);
        SEQAN_ASSERT_EQ(pos[5], 5u);
        SEQAN_ASSERT_EQ(pos[6], 8u);
    }
    
    //____________________________________________________________________________
    // Test4 - different alphabet, large ndl

    {
        typedef String<Dna>                                 THaystack;
        typedef String<Dna>                                 TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;    
        typedef Finder<TIndex>                              TFinder;
        typedef typename DefaultIndexPattern<TNeedle, TIndex>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;

        TPositions pos;
        THaystack hstck("taaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaat");
        TIndex index(hstck);
        configureIndex(index);
        TFinder finder(index);
        TNeedle ndl("taaaataaaataaaataaaataaaataaaataaaataaaataaaat");
        TPattern pattern(ndl);

        while (find(finder, pattern)) {
            append(pos, position(finder));
            SEQAN_ASSERT_EQ(position(finder), beginPosition(finder));
            SEQAN_ASSERT_EQ(endPosition(finder), beginPosition(finder) + length(finder));
            SEQAN_ASSERT_EQ(length(finder), length(ndl));
            SEQAN_ASSERT_EQ(infix(finder), ndl);
        }

        if (! MatchesAreSorted<TFinder>::Type::VALUE) {
            std::sort(begin(pos), end(pos));
        }
        SEQAN_ASSERT_EQ(length(pos), 6u);
        SEQAN_ASSERT_EQ(pos[0], 0u);
        SEQAN_ASSERT_EQ(pos[1], 5u);
        SEQAN_ASSERT_EQ(pos[2], 10u);
        SEQAN_ASSERT_EQ(pos[3], 15u);
        SEQAN_ASSERT_EQ(pos[4], 20u);
        SEQAN_ASSERT_EQ(pos[5], 25u);
    }
    
    //____________________________________________________________________________
    // Test5 - fasta file

    {
        typedef String<Dna>                                 THaystack;
        typedef String<Dna>                                 TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;    
        typedef Finder<TIndex>                              TFinder;
        typedef typename DefaultIndexPattern<TNeedle, TIndex>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;

        String<char> testFile = getAbsolutePath("/tests/find_index_approx/fasta-small.fasta");

        SeqFileIn seqIO(toCString(testFile));
        String<char> id;
        THaystack hstck;
        //int res = 
        readRecord(id, hstck, seqIO);
        //SEQAN_ASSERT_EQ_MSG(res, 0, "Could not read file.");
        
        TPositions pos;
        TIndex index(hstck);
        configureIndex(index);
        TFinder finder(index);
        TNeedle ndl("TTACTTT");
        TPattern pattern(ndl);

        while (find(finder, pattern)) {
            append(pos, position(finder));
            SEQAN_ASSERT_EQ(position(finder), beginPosition(finder));
            SEQAN_ASSERT_EQ(endPosition(finder), beginPosition(finder) + length(finder));
            SEQAN_ASSERT_EQ(length(finder), length(ndl));
            SEQAN_ASSERT_EQ(infix(finder), ndl);
        }

        if (! MatchesAreSorted<TFinder>::Type::VALUE) {
            std::sort(begin(pos), end(pos));
        }
        SEQAN_ASSERT_EQ(length(pos), 8u);
        SEQAN_ASSERT_EQ(pos[0], 665u);
        SEQAN_ASSERT_EQ(pos[1], 721u);
        SEQAN_ASSERT_EQ(pos[2], 741u);
        SEQAN_ASSERT_EQ(pos[3], 763u);
        SEQAN_ASSERT_EQ(pos[4], 847u);
        SEQAN_ASSERT_EQ(pos[5], 1021u);
        SEQAN_ASSERT_EQ(pos[6], 1263u);
        SEQAN_ASSERT_EQ(pos[7], 1537u);
    }
}

SEQAN_DEFINE_TEST(test_find_index_esa) {
    testFindIndex<IndexEsa<> >();
}

SEQAN_DEFINE_TEST(test_find_index_wotd) {
    testFindIndex<IndexWotd<> >();
}

SEQAN_DEFINE_TEST(test_find_index_fm) {
    testFindIndex<FMIndex<> >();
}

SEQAN_DEFINE_TEST(test_find_index_qgram) {
    testFindIndex<IndexQGram<UngappedShape<2> > >();
}

SEQAN_DEFINE_TEST(test_find_index_qsample) {
    // q-sample index with step size 2 (see configureIndex above)
    testFindIndex<IndexQGram<UngappedShape<3> > >();
}

SEQAN_DEFINE_TEST(test_find_index_lz) {
    testFindIndex<IndexLZ<> >();
}

SEQAN_DEFINE_TEST(test_find_index_sadakane) {
    testFindIndex<IndexSadakane<> >();
}

SEQAN_DEFINE_TEST(test_find_index_sttd64) {
    testFindIndex<IndexSttd64<> >();
}

SEQAN_DEFINE_TEST(test_find_index_qgram2l) {
    testFindIndex<IndexQGram2L<UngappedShape<2> > >();
}

SEQAN_DEFINE_TEST(test_find_index_digest) {
    // TODO(krugel) DiGeST Config
    const unsigned PARTITION_SIZE = 5;
    const unsigned OUTBUF_SIZE = 5;
    const unsigned INBUF_SIZE = 10;
    const unsigned TAIL_LENGTH = 3;
    const unsigned PREFIX_LENGTH = 2;

    typedef DigestConfig<PARTITION_SIZE, OUTBUF_SIZE, INBUF_SIZE, TAIL_LENGTH, PREFIX_LENGTH> TConfig;

    testFindIndex<IndexDigest<TConfig> >();
}

}  // namespace seqan

#endif //#ifndef TESTS_INDEX2_TEST_FIND_INDEX_H_
