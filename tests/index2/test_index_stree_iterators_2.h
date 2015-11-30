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

// File adapted from /tests/index/test_stree_iterators.h

/*
+ goDown()
+ goDown(char)
+ goDown(text)
+ isRoot
+ repLength
+ representative

+ isLeaf
+ goRoot
+ countOccurrences
+ countChildren
+ getOccurrences
+ getOccurrence

/ nodeDepth
/ emptyParentEdge
/ parentEdgeFirstChar
/ parentEdgeLabel
/ parentEdgeLength
/ parentRepLength
*/

#ifndef TESTS_INDEX_TEST_STREE_ITERATORS_2_H_
#define TESTS_INDEX_TEST_STREE_ITERATORS_2_H_

#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/file.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/index2.h>

namespace seqan {

template <typename TIndexSpec>
void testSTreeIterators()
{
    typedef Index<String<char>, TIndexSpec> TIndex;
/*
    typedef typename Iterator<TIndex, TopDown<> >::Type TIterator;

    // test empty trees
    {
        TIndex index("");
        TIterator iter(index);
        TIterator piter(index);
        SEQAN_ASSERT_NOT(goDown(iter, "a"));
    }
    
    {
        String<char> text("acaaacatatz");
        TIndex index(text);
        Iter<TIndex, VSTree<TopDown<> > > it(index);
        
        while (goDown(it, "a")) {
            //std::cout << countOccurrences(it) << "\t";
            std::cerr << representative(it) << "\t";
            //std::cerr << nodeDepth(it) << "\t";
            std::cerr << std::endl;
        }
        std::cerr << std::endl;
    }
*/
    {
        String<char> text("acaaacatatz");
        TIndex index(text);
        Iter<TIndex, VSTree<TopDown<> > > it(index);
        
        SEQAN_ASSERT_EQ(isRoot(it), true);
        SEQAN_ASSERT_EQ(isLeaf(it), false);
        SEQAN_ASSERT_EQ(representative(it), "");
        SEQAN_ASSERT_EQ(repLength(it), 0u);
        //SEQAN_ASSERT_EQ(countChildren(it), 4u); // TODO(krugel): Counting the children of the root doesn't work for IndexEsa?
        SEQAN_ASSERT_EQ(countOccurrences(it), 11u);
        SEQAN_ASSERT_EQ(length(getOccurrences(it)), 11u);

        SEQAN_ASSERT(goDown(it));
        SEQAN_ASSERT_EQ(representative(it), "a");
        SEQAN_ASSERT_EQ(repLength(it), 1u);
        SEQAN_ASSERT_EQ(isRoot(it), false);
        SEQAN_ASSERT_EQ(isLeaf(it), false);
        SEQAN_ASSERT_EQ(countChildren(it), 3u);
        SEQAN_ASSERT_EQ(countOccurrences(it), 6u);
        SEQAN_ASSERT_EQ(length(getOccurrences(it)), 6u);

        SEQAN_ASSERT(goDown(it, "a"));
        SEQAN_ASSERT_EQ(representative(it), "aa");
        SEQAN_ASSERT_EQ(repLength(it), 2u);
        SEQAN_ASSERT_EQ(isRoot(it), false);
        SEQAN_ASSERT_EQ(isLeaf(it), false);
        SEQAN_ASSERT_EQ(countChildren(it), 2u);
        SEQAN_ASSERT_EQ(countOccurrences(it), 2u);
        SEQAN_ASSERT_EQ(length(getOccurrences(it)), 2u);
        SEQAN_ASSERT(getOccurrence(it) == 2u || getOccurrence(it) == 3u);

        SEQAN_ASSERT(goDown(it, 'a'));
        SEQAN_ASSERT_EQ(representative(it), "aaacatatz");
        SEQAN_ASSERT_EQ(repLength(it), 9u);
        SEQAN_ASSERT_EQ(isRoot(it), false);
        SEQAN_ASSERT_EQ(isLeaf(it), true);
        SEQAN_ASSERT_EQ(countChildren(it), 0u);
        SEQAN_ASSERT_EQ(countOccurrences(it), 1u);
        SEQAN_ASSERT_EQ(length(getOccurrences(it)), 1u);
        SEQAN_ASSERT_EQ(getOccurrence(it), 2u);
        
        SEQAN_ASSERT_NOT(goDown(it));
        
        goRoot(it);
        SEQAN_ASSERT_EQ(isRoot(it), true);

        SEQAN_ASSERT_EQ(length(index), length(text));
        SEQAN_ASSERT_EQ(indexText(index), text);
        //_dump(index);
    }
}

template <typename TIndex1, typename TIndex2>
void compareTreeIteratorsParentLinks(TIndex1 &index1, TIndex2 &index2) {
    Iter<TIndex1, VSTree<TopDown<ParentLinks<Preorder> > > > it1(index1);
    Iter<TIndex2, VSTree<TopDown<ParentLinks<Preorder> > > > it2(index2);

    while (!atEnd(it1) && !atEnd(it2)) 
    {
        SEQAN_ASSERT_EQ(representative(it1), representative(it2));
        SEQAN_ASSERT_EQ(parentEdgeLabel(it1), parentEdgeLabel(it2));
        SEQAN_ASSERT_EQ(countOccurrences(it1), countOccurrences(it2));
        SEQAN_ASSERT_EQ(isRoot(it1), isRoot(it2));
        //SEQAN_ASSERT_EQ(isLeaf(it1), isLeaf(it2));
        goNext(it1);
        goNext(it2);
    }

    SEQAN_ASSERT_EQ(atEnd(it1), atEnd(it1));
}

template <typename TIndex1, typename TIndex2, typename TString>
void compareTreeIteratorsGoDown(TIndex1 &index1, TIndex2 &index2, TString needle) {
    Iter<TIndex1, VSTree<TopDown<> > > it1(index1);
    Iter<TIndex2, VSTree<TopDown<> > > it2(index2);

    for (unsigned int i = 0u; i < length(needle); ++i) {
        goDown(it1, needle[i]);
        goDown(it2, needle[i]);

        SEQAN_ASSERT_EQ(representative(it1), representative(it2));
        SEQAN_ASSERT_EQ(parentEdgeLabel(it1), parentEdgeLabel(it2));
        SEQAN_ASSERT_EQ(countOccurrences(it1), countOccurrences(it2));
        SEQAN_ASSERT_EQ(isRoot(it1), isRoot(it2));
        //SEQAN_ASSERT_EQ(isLeaf(it1), isLeaf(it2));
    }
    SEQAN_ASSERT_EQ(atEnd(it1), atEnd(it1));
}

template <typename TIndexSpec1, typename TIndexSpec2>
void compareIndices() {
    {
        CharString text("mississippi");
        Index<CharString, TIndexSpec1> index1(text);
        Index<CharString, TIndexSpec2> index2(text);
        compareTreeIteratorsGoDown(index1, index2, "issi");
        //compareTreeIteratorsParentLinks(index1, index2);
    }
    {
        DnaString text("acaaacatat");
        Index<DnaString, TIndexSpec1> index1(text);
        Index<DnaString, TIndexSpec2> index2(text);
        compareTreeIteratorsGoDown(index1, index2, "aca");
        //compareTreeIteratorsParentLinks(index1, index2);
    }
    {
        String<char> testFile = getAbsolutePath("/tests/index2/fasta-small.fasta");
        SeqFileIn seqIO(toCString(testFile));
        String<char> id;
        DnaString text;
        //int res =
        readRecord(id, text, seqIO);
        //SEQAN_ASSERT_EQ_MSG(res, 0, "Could not read file.");

        Index<DnaString, TIndexSpec1> index1(text);
        Index<DnaString, TIndexSpec2> index2(text);
        compareTreeIteratorsGoDown(index1, index2, "TCTTAAAGGGTAAGAATTTAAGGTTCTTATTGAATTATTGTTACTATTTTAGAACTATTGTTGCTGTTTA");
        //compareTreeIteratorsParentLinks(index1, index2);
    }
}

SEQAN_DEFINE_TEST(test_index_stree_iterators_2_esa)
{
    testSTreeIterators<IndexEsa<> >();
}

SEQAN_DEFINE_TEST(test_index_stree_iterators_2_wotd)
{
    testSTreeIterators<IndexWotd<> >();
}

SEQAN_DEFINE_TEST(test_index_stree_iterators_2_fm)
{
    testSTreeIterators<FMIndex<> >();
}

SEQAN_DEFINE_TEST(test_index_stree_iterators_2_sttd64)
{
    testSTreeIterators<IndexSttd64<> >();
}

SEQAN_DEFINE_TEST(test_index_stree_iterators_2_compare_esa_wotd)
{
    compareIndices<IndexEsa<>, IndexWotd<> >();
}

SEQAN_DEFINE_TEST(test_index_stree_iterators_2_compare_wotd_sttd64)
{
    compareIndices<IndexWotd<>, IndexSttd64<> >();
}

}  // namespace seqan

#endif //TESTS_INDEX_TEST_STREE_ITERATORS_2_H_
