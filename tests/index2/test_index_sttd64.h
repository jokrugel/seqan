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
// Author: Alexander Aumann <a.aumann@web.de>
// ==========================================================================
// Tests for the SeqAn module index2.
// ==========================================================================

#ifndef TESTS_INDEX_TEST_INDEX_STTD64_H_
#define TESTS_INDEX_TEST_INDEX_STTD64_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/index2.h>

namespace seqan {

template <typename TNodeTable>
bool areNodeTablesEqual(TNodeTable const & tableOne,
                            TNodeTable const & tableTwo) {
    typename Iterator<TNodeTable const>::Type tOneIt;
    typename Iterator<TNodeTable const>::Type tTwoIt;
    Sttd64Node_ tOneNode, tTwoNode;

    if (length(tableOne) != length(tableTwo)) {
        std::cout << "Node tables are not equally long!\n";
        std::cout << "T1: " << length(tableOne) << ". T2: " << length(tableTwo) << ".\n";
        return false;
    }
    int i = 0;
    for (tOneIt = begin(tableOne), tTwoIt = begin(tableTwo);
            tOneIt != end(tableOne);
            tOneIt++, tTwoIt++, i++) {

            tOneNode = *tOneIt;
            tTwoNode = *tTwoIt;

            if (tOneNode.lp != tTwoNode.lp || tOneNode.pointerOrDepth != tTwoNode.pointerOrDepth) {
                std::cout << "Node table inequality at position " << i << "!" << std::endl;
                std::cout << "Left node:\t" << tOneNode.toString() << std::endl;
                std::cout << "Right node:\t" << tTwoNode.toString() << std::endl;
                return false;
            }

    }
    return true;
}

template<typename TNodeTable>
void parseTreeTable(std::string const & table, TNodeTable & nodeTable) {
    std::string::const_iterator strIt;

    Sttd64Node_ node;
    node.lp = 0;
    node.pointerOrDepth = 0;
    enum STATES {PARSE_INT_ONE, PARSE_INT_TWO, PARSE_TYPE, PARSE_RIGHTMOST, PARSE_IN_BETWEEN,
                    PARSE_OPEN_BRACKET_INTS};
    STATES currentState = PARSE_IN_BETWEEN;
    std::stringstream intStream("");
    clear(nodeTable);
    int nodeInt;
    for (strIt = table.begin(); strIt != table.end(); strIt++) {
        if (*strIt == '[') {
            SEQAN_ASSERT(currentState == PARSE_IN_BETWEEN);
            currentState = PARSE_TYPE;
        } else if (*strIt == ',') {
            SEQAN_ASSERT(currentState == PARSE_INT_ONE);
            SEQAN_ASSERT(intStream.good());
            intStream << '\0';
            intStream >> nodeInt;
            SEQAN_ASSERT_MSG(intStream.good(), "Error while parsing int one!");
            intStream.str("");
            node.lp = nodeInt;
            currentState = PARSE_INT_TWO;
        } else if (*strIt == '(') {
            SEQAN_ASSERT(currentState == PARSE_OPEN_BRACKET_INTS);
            currentState = PARSE_INT_ONE;
        } else if (*strIt == ']') {
            SEQAN_ASSERT(currentState == PARSE_IN_BETWEEN || currentState == PARSE_RIGHTMOST);
            appendValue(nodeTable, node);
            node.lp = 0;
            node.pointerOrDepth = 0;
            currentState = PARSE_IN_BETWEEN;
        } else if (*strIt == 'N'|| *strIt == 'L') {
            SEQAN_ASSERT(currentState == PARSE_TYPE);
            if (*strIt == 'L') {
                node.setLeafBit(true);
            }
            currentState = PARSE_OPEN_BRACKET_INTS;
        } else if (*strIt == ')') {
            SEQAN_ASSERT(currentState == PARSE_INT_TWO);
            intStream << '\0';
            SEQAN_ASSERT(intStream.good());
            intStream >> nodeInt;
            SEQAN_ASSERT_MSG(intStream.good(), "Error while parsing int two!");
            intStream.str("");
            node.pointerOrDepth += nodeInt;
            currentState = PARSE_RIGHTMOST;
        } else if (*strIt == 'R') {
            SEQAN_ASSERT(currentState == PARSE_RIGHTMOST);
            node.setRightMostBit(true);
            currentState = PARSE_IN_BETWEEN;
        } else {
            SEQAN_ASSERT(currentState == PARSE_INT_ONE || currentState == PARSE_INT_TWO);
            SEQAN_ASSERT(intStream.good());
            intStream << *strIt;
        }
    }
}

template <typename TConfig, typename TText>
bool compareIndexTablesToRefs(Index<TText, IndexSttd64<TConfig> > const & index,
                                std::string const & compTableA,
                                std::string const & compTableC,
                                std::string const & compTableG,
                                std::string const & compTableT) {
    typedef Index<TText, IndexSttd64<TConfig> > TIndex;
    typename Iterator<typename TIndex::TTrees const>::Type treesIt;

    treesIt = begin(index.trees);
    bool treesComparedEqual;
    typename TIndex::TPartTree currentIndexTree;
    typename TIndex::TPartTree currentCompTree;
    //currentIndexTree = *treesIt;
    if (empty(*treesIt)) {
        if (!compTableA.empty()) {
            std::cout << "Tables for A were not equal, tree was empty, comparison string not!" << std::endl;
            return false;
        }
    } else {
        parseTreeTable(compTableA, currentCompTree);
        treesComparedEqual = areNodeTablesEqual(*treesIt, currentCompTree);
        if (!treesComparedEqual) {
            std::cout << "Tables for A were not equal, print tables (tree first, then comparison)."
                    << std::endl;
            _printRawTreeTable(*treesIt);
            _printRawTreeTable(currentCompTree);
            return false;
        }
    }
    treesIt++;

    //table C
    if (empty(*treesIt)) {
        if (!compTableC.empty()) {
            std::cout << "Tables for C were not equal, tree was empty, comparison string not!" << std::endl;
            return false;
        }
    } else {
        parseTreeTable(compTableC, currentCompTree);
        treesComparedEqual = areNodeTablesEqual(*treesIt, currentCompTree);
        if (!treesComparedEqual) {
            std::cout << "Tables for C were not equal, print tables (tree first, then comparison)."
                    << std::endl;
            _printRawTreeTable(*treesIt);
            _printRawTreeTable(currentCompTree);
            return false;
        }
    }
    treesIt++;

    //table G
    if (empty(*treesIt)) {
        if (!compTableG.empty()) {
            std::cout << "Tables for G were not equal, tree was empty, comparison string not!" << std::endl;
            return false;
        }
    } else {
        parseTreeTable(compTableG, currentCompTree);
        treesComparedEqual = areNodeTablesEqual(*treesIt, currentCompTree);
        if (!treesComparedEqual) {
            std::cout << "Tables for G were not equal, print tables (tree first, then comparison)."
                    << std::endl;
            _printRawTreeTable(*treesIt);
            _printRawTreeTable(currentCompTree);
            return false;
        }
    }
    treesIt++;

    //table T
    if (empty(*treesIt)) {
        if (!compTableT.empty()) {
            std::cout << "Tables for T were not equal, tree was empty, comparison string not!" << std::endl;
            return false;
        }
    } else {
        parseTreeTable(compTableT, currentCompTree);
        treesComparedEqual = areNodeTablesEqual(*treesIt, currentCompTree);
        if (!treesComparedEqual) {
            std::cout << "Tables for T were not equal, print tables (tree first, then comparison)."
                    << std::endl;
            _printRawTreeTable(*treesIt);
            _printRawTreeTable(currentCompTree);
            return false;
        }
    }

    return true;
}

SEQAN_DEFINE_TEST(test_index_sttd64_basicConstruction1)
{
    using namespace seqan;
    typedef Index<String<Dna5>, IndexSttd64<> > TDna5Index;

    String<Dna5> str1 = "AGATAG";

    TDna5Index index(str1);
    createPartitions(index);
    createSuffixTrees(index);
    //a partition length of 1 is used for these tests
    std::string compTableA("[N(0,1)R][N(1,3)][L(3,0)R][L(2,0)][L(6,0)R]");
    std::string compTableC("");
    std::string compTableG("[N(1,1)R][L(2,0)][L(6,0)R]");
    std::string compTableT("[L(3,0)R]");

    SEQAN_ASSERT(compareIndexTablesToRefs(index, compTableA, compTableC, compTableG, compTableT));
}

SEQAN_DEFINE_TEST(test_index_sttd64_basicConstruction2)
{
    using namespace seqan;
    typedef Index<String<Dna5>, IndexSttd64<> > TDna5Index;

    String<Dna5> str = "AGAGAGCTT";

    TDna5Index index(str);

    createPartitions(index);
    createSuffixTrees(index);
    std::string compTableA("[N(0,1)R][N(2,3)][L(6,0)R][L(4,0)][L(6,0)R]");
    std::string compTableC("[L(6,0)R]");
    std::string compTableG("[N(1,1)R][N(2,3)][L(6,0)R][L(4,0)][L(6,0)R]");
    std::string compTableT("[N(7,1)R][L(8,0)][L(9,0)R]");
    SEQAN_ASSERT(compareIndexTablesToRefs(index, compTableA, compTableC, compTableG, compTableT));
}

SEQAN_DEFINE_TEST(test_index_sttd64_basicConstruction3)
{
    using namespace seqan;
    typedef Index<String<Dna5>, IndexSttd64<> > TDna5Index;

    String<Dna5> str3 = "AGAGACAGGCAGGG";
    String<Dna5> str4 = "TGTGTCTGGCTGGG";
    TDna5Index index(str3);

    createPartitions(index);
    createSuffixTrees(index);
    std::string compTableA("[N(0,1)R][N(1,3)][L(5,0)R][N(2,7)][N(8,5)R][L(9,0)][L(13,0)R][L(3,0)][L(5,0)R]");
    std::string compTableC("[N(5,1)R][L(9,0)][L(13,0)R]");
    std::string compTableG("[N(1,1)R][N(2,8)][N(8,5)][L(9,0)][L(14,0)R][L(9,0)][L(13,0)][L(14,0)R][L(3,0)][L(5,0)R]");
    std::string compTableT("");
    SEQAN_ASSERT(compareIndexTablesToRefs(index, compTableA, compTableC, compTableG, compTableT));

    //the node tables should be completely identical for str3 and str4, after exchanging A and T
    //clear(index);
    TDna5Index index2(str4);
    createPartitions(index2);
    createSuffixTrees(index2);
    SEQAN_ASSERT(compareIndexTablesToRefs(index2, compTableT, compTableC, compTableG, compTableA));
}

SEQAN_DEFINE_TEST(test_index_sttd64_basicConstruction4)
{
    using namespace seqan;
    typedef Index<String<Dna5>, IndexSttd64<> > TDna5Index;

    String<Dna5> str = "GATATA";

    TDna5Index index(str);

    createPartitions(index);
    createSuffixTrees(index);
    std::string compTableA("[N(1,1)R][N(2,3)][L(6,0)R][L(4,0)][L(6,0)R]");
    std::string compTableC("");
    std::string compTableG("[L(0,0)R]");
    std::string compTableT = ("[N(2,1)R][L(4,0)][L(6,0)R]");
    SEQAN_ASSERT(compareIndexTablesToRefs(index, compTableA, compTableC, compTableG, compTableT));
}

SEQAN_DEFINE_TEST(test_index_sttd64_printTree)
{
    using namespace seqan;
    typedef Index<String<Dna5>, IndexSttd64<> > TDna5Index;

    //String<Dna5> str = "AGAGAGCTT"; //this is the origin string for these tables
    std::string compTableA("[N(0,1)R][N(2,3)][L(6,0)R][L(4,0)][L(6,0)R]");
    std::string compTableC("[L(6,0)R]");
    std::string compTableG("[N(1,1)R][N(2,3)][L(6,0)R][L(4,0)][L(6,0)R]");
    std::string compTableT("[N(7,1)R][L(8,0)][L(9,0)R]");

    std::stringstream strStr;
    std::string printedTable;

    TDna5Index::TPartTree currentPrintedTree;
    TDna5Index::TPartTree currentCompTree;
    parseTreeTable(compTableA, currentCompTree);
    _printRawTreeTable(currentCompTree, strStr);
    strStr >> printedTable;

    parseTreeTable(printedTable, currentPrintedTree);
    SEQAN_ASSERT(areNodeTablesEqual(currentCompTree, currentPrintedTree));

    strStr.str("");
    parseTreeTable(compTableC, currentCompTree);
    _printRawTreeTable(currentCompTree, strStr);
    strStr >> printedTable;

    parseTreeTable(printedTable, currentPrintedTree);
    SEQAN_ASSERT(areNodeTablesEqual(currentCompTree, currentPrintedTree));

    strStr.str("");
    parseTreeTable(compTableG, currentCompTree);
    _printRawTreeTable(currentCompTree, strStr);
    strStr >> printedTable;

    parseTreeTable(printedTable, currentPrintedTree);
    SEQAN_ASSERT(areNodeTablesEqual(currentCompTree, currentPrintedTree));

    strStr.str("");
    parseTreeTable(compTableT, currentCompTree);
    _printRawTreeTable(currentCompTree, strStr);
    strStr >> printedTable;

    parseTreeTable(printedTable, currentPrintedTree);
    SEQAN_ASSERT(areNodeTablesEqual(currentCompTree, currentPrintedTree));

    std::string compNicelyPrintedTreeA("(((0)(2))(4))");
    std::string compNicelyPrintedTreeC("(6)");
    std::string compNicelyPrintedTreeG("(((1)(3))(5))");
    std::string compNicelyPrintedTreeT("((7)(8))");

    strStr.str("");
    parseTreeTable(compTableA, currentPrintedTree);
    _printPartTree(currentPrintedTree, false, strStr);
    strStr >> printedTable;
    SEQAN_ASSERT_EQ(compNicelyPrintedTreeA , printedTable);

    // TODO(krugel) Check the following tests
    // Fails in parseTreeTable: appendValue(nodeTable, node);

    //strStr.str("");
    //parseTreeTable(compTableC, currentPrintedTree);
    //_printPartTree(currentPrintedTree, false, strStr);
    //strStr >> printedTable;
    //SEQAN_ASSERT_EQ(compNicelyPrintedTreeC , printedTable);

    //strStr.str("");
    //parseTreeTable(compTableG, currentPrintedTree);
    //_printPartTree(currentPrintedTree, false, strStr);
    //strStr >> printedTable;
    //SEQAN_ASSERT_EQ(compNicelyPrintedTreeG , printedTable);

    //strStr.str("");
    //parseTreeTable(compTableT, currentPrintedTree);
    //_printPartTree(currentPrintedTree, false, strStr);
    //strStr >> printedTable;
    //SEQAN_ASSERT_EQ(compNicelyPrintedTreeT , printedTable);
}

// TODO(krugel) Use different file
SEQAN_DEFINE_TEST(test_index_sttd64_readFile)
{
    using namespace seqan;
    typedef Index<String<Dna5>, IndexSttd64<> > TDna5Index;
    TDna5Index memIndex("AGAGAGCTT");
    createPartitions(memIndex);
    createSuffixTrees(memIndex);

    typedef Index<String<Dna5>, IndexSttd64<Sttd64Config<1, 8192, 2> > > TFileIndex;
    String<Dna5> txt;
    char const* const sufTreePath = "/home/aleaum/DNA/fasta_t";
    //char const* outDirName = "/home/aleaum/DNA2/";

    String<char, MMap<> > mmapString;
    //try to load the externally stored input text
    if (!open(mmapString, sufTreePath, OPEN_RDONLY))
    {
        std::cerr << "Could not open string file!\n";
    }

    AutoSeqFormat format;
    guessFormat(mmapString, format);
    resize(txt, length(mmapString), Exact());
    assignSeq(txt, mmapString, format);
    close(mmapString);

    TFileIndex extIndex(txt);
    extIndex.inMem = false;
    
    createPartitions(extIndex);
    createSuffixTrees(extIndex);
    std::cout << "Creation complete!";

    _reopenTreeFile(extIndex, 0);
    SEQAN_ASSERT(areNodeTablesEqual(getValue(memIndex.trees, 0), getValue(extIndex.trees, 0)));
    _reopenTreeFile(extIndex, 1);
    SEQAN_ASSERT(areNodeTablesEqual(getValue(memIndex.trees, 1), getValue(extIndex.trees, 0)));
    _reopenTreeFile(extIndex, 2);
    SEQAN_ASSERT(areNodeTablesEqual(getValue(memIndex.trees, 2), getValue(extIndex.trees, 0)));
    _reopenTreeFile(extIndex, 3);
    SEQAN_ASSERT(areNodeTablesEqual(getValue(memIndex.trees, 3), getValue(extIndex.trees, 0)));
//    compareIndexTablesToRefs(extIndex, getValue(memIndex.trees, 0),
//                             getValue(memIndex.trees, 1),
//                             getValue(memIndex.trees, 2),
//                             getValue(memIndex.trees, 3));
}

SEQAN_DEFINE_TEST(test_index_sttd64_iterateTopDown1)
{
    using namespace seqan;
    typedef Index<String<Dna5>, IndexSttd64<> > TIndex;
    TIndex index("AGAGAGCTT");
    Iterator< TIndex, TopDown<> >::Type it(index);

    _printRawTreeTable(getValue(index.trees, 0), std::cerr);
    _printRawTreeTable(getValue(index.trees, 1), std::cerr);
    _printRawTreeTable(getValue(index.trees, 2), std::cerr);
    _printRawTreeTable(getValue(index.trees, 3), std::cerr);
//    _printPartTree(getValue(index.trees, 0), true, std::cerr);
//    _printPartTree(getValue(index.trees, 1), true, std::cerr);
//    _printPartTree(getValue(index.trees, 2), true, std::cerr);
//    _printPartTree(getValue(index.trees, 3), true, std::cerr);

    String<Dna5> pattern = "GAG";
    bool notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        //std::cout << representative(it) << std::endl;
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    // output match positions
    String<typename SAValue<TIndex>::Type > occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 2u);
    SEQAN_ASSERT(contains(occurrences, 1u));
    SEQAN_ASSERT(contains(occurrences, 3u));

}

SEQAN_DEFINE_TEST(test_index_sttd64_iterateTopDown2)
{
    using namespace seqan;
    typedef Index<CharString, IndexSttd64<> > TIndex;
    TIndex index("How many wood would a woodchuck chuck.");
    Iterator< TIndex, TopDown<> >::Type it(index);

    CharString pattern = "wood";
    bool notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    // output match positions
    String<typename SAValue<TIndex>::Type > occurrences = getOccurrences(it);
//    for (int i = 0; i < length(occurrences); i++)
//    {
//        std::cout << occurrences[i] << ",";
//    }
//    std::cout << "\n";
    SEQAN_ASSERT_EQ(length(occurrences) , 2u);
    SEQAN_ASSERT(contains(occurrences, 9u));
    SEQAN_ASSERT(contains(occurrences, 22u));

    pattern = "How many wood would a woodchuck chuck.";
    notFound = false;
    clear(it);
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 1u);
    SEQAN_ASSERT(contains(occurrences, 0u));

    pattern = "woodchuck chuck.";
    notFound = false;
    clear(it);
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 1u);
    SEQAN_ASSERT(contains(occurrences, 22u));
}

SEQAN_DEFINE_TEST(test_index_sttd64_iterateTopDown3)
{
    using namespace seqan;
    typedef Index<CharString, IndexSttd64<> > TIndex;
    TIndex index("How many wood would a woodchuck chuck.");
    Iterator< TIndex, TopDown<> >::Type it(index);

    CharString pattern = ".";
    bool notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    // output match positions
    String<typename SAValue<TIndex>::Type > occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences), 1u);
    SEQAN_ASSERT(contains(occurrences, 37u));
}

SEQAN_DEFINE_TEST(test_index_sttd64_iterateTopDown4)
{
    using namespace seqan;
    typedef Index<String<Dna5>, IndexSttd64<> > TIndex;
    TIndex index("AGAGAGCTT");
    Iterator< TIndex, TopDown<> >::Type it(index);

    String<Dna5> pattern = "T";
    bool notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        std::cout << representative(it) << std::endl;
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    // output match positions
    String<typename SAValue<TIndex>::Type > occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 2u);
    SEQAN_ASSERT(contains(occurrences, 7u));
    SEQAN_ASSERT(contains(occurrences, 8u));
}

SEQAN_DEFINE_TEST(test_index_sttd64_iterateTopDown5)
{
    using namespace seqan;
    typedef Index<String<Dna5>, IndexSttd64<> > TIndex;
    TIndex index("AGAGACAGGCAGGG");
    Iterator< TIndex, TopDown<> >::Type it(index);

    String<Dna5> pattern = "AG";
    bool notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    String<typename SAValue<TIndex>::Type > occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 4u);
    SEQAN_ASSERT(contains(occurrences, 0u));
    SEQAN_ASSERT(contains(occurrences, 2u));
    SEQAN_ASSERT(contains(occurrences, 6u));
    SEQAN_ASSERT(contains(occurrences, 10u));


    clear(it);
    pattern = "G";
    notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);

    occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 7u);
    SEQAN_ASSERT(contains(occurrences, 1u));
    SEQAN_ASSERT(contains(occurrences, 3u));
    SEQAN_ASSERT(contains(occurrences, 7u));
    SEQAN_ASSERT(contains(occurrences, 8u));
    SEQAN_ASSERT(contains(occurrences, 11u));
    SEQAN_ASSERT(contains(occurrences, 12u));
    SEQAN_ASSERT(contains(occurrences, 13u));

    clear(it);
    pattern = "CAGG";
    notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        std::cout << representative(it) << std::endl;
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 2u);
    SEQAN_ASSERT(contains(occurrences, 5u));
    SEQAN_ASSERT(contains(occurrences, 9u));
}

SEQAN_DEFINE_TEST(test_index_sttd64_iterateTopDownPrefLen2)
{
    using namespace seqan;
    typedef Index<String<Dna5>, IndexSttd64<Sttd64Config<2> > > TIndex;
    TIndex index("AGAGACAGGCAGGG");
    Iterator< TIndex, TopDown<> >::Type it(index);

    String<Dna5> pattern = "AG";
    bool notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    String<typename SAValue<TIndex>::Type > occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 4u);
    SEQAN_ASSERT(contains(occurrences, 0u));
    SEQAN_ASSERT(contains(occurrences, 2u));
    SEQAN_ASSERT(contains(occurrences, 6u));
    SEQAN_ASSERT(contains(occurrences, 10u));


    clear(it);
    pattern = "G";
    notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);

    occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 7u);
    SEQAN_ASSERT(contains(occurrences, 1u));
    SEQAN_ASSERT(contains(occurrences, 3u));
    SEQAN_ASSERT(contains(occurrences, 7u));
    SEQAN_ASSERT(contains(occurrences, 8u));
    SEQAN_ASSERT(contains(occurrences, 11u));
    SEQAN_ASSERT(contains(occurrences, 12u));
    SEQAN_ASSERT(contains(occurrences, 13u));

    clear(it);
    pattern = "CAGG";
    notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 2u);
    SEQAN_ASSERT(contains(occurrences, 5u));
    SEQAN_ASSERT(contains(occurrences, 9u));
}

SEQAN_DEFINE_TEST(test_index_sttd64_iterateTopDownPrefLen3)
{
    using namespace seqan;
    typedef Index<String<Dna5>, IndexSttd64<Sttd64Config<3> > > TIndex;
    TIndex index("AGAGACAGGCAGGG");
    Iterator< TIndex, TopDown<> >::Type it(index);

    String<Dna5> pattern = "AG";
    bool notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    String<typename SAValue<TIndex>::Type > occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 4u);
    SEQAN_ASSERT(contains(occurrences, 0u));
    SEQAN_ASSERT(contains(occurrences, 2u));
    SEQAN_ASSERT(contains(occurrences, 6u));
    SEQAN_ASSERT(contains(occurrences, 10u));


    clear(it);
    pattern = "G";
    notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);

    occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 7u);
    SEQAN_ASSERT(contains(occurrences, 1u));
    SEQAN_ASSERT(contains(occurrences, 3u));
    SEQAN_ASSERT(contains(occurrences, 7u));
    SEQAN_ASSERT(contains(occurrences, 8u));
    SEQAN_ASSERT(contains(occurrences, 11u));
    SEQAN_ASSERT(contains(occurrences, 12u));
    SEQAN_ASSERT(contains(occurrences, 13u));

    clear(it);
    pattern = "CAGG";
    notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 2u);
    SEQAN_ASSERT(contains(occurrences, 5u));
    SEQAN_ASSERT(contains(occurrences, 9u));
}

SEQAN_DEFINE_TEST(test_index_sttd64_iterateTopDownPrefLen2Files)
{
    using namespace seqan;
    typedef Index<String<Dna>, IndexSttd64<Sttd64Config<2> > > TIndex;

    String<Dna> txt;
    char const* const sufTreePath = "/home/aleaum/DNA/fasta_t_2";
    //char const* outDirName = "/home/aleaum/DNA/";

    String<char, MMap<> > mmapString;
    //try to load the externally stored input text
    if (!open(mmapString, sufTreePath, OPEN_RDONLY))
    {
        std::cerr << "Could not open string file!\n";
    }

    AutoSeqFormat format;
    guessFormat(mmapString, format);
    resize(txt, length(mmapString), Exact());
    assignSeq(txt, mmapString, format);
    close(mmapString);

    TIndex index(txt);
    index.inMem = false;

    Iterator< TIndex, TopDown<> >::Type it(index);

    String<Dna> pattern = "AG";
    bool notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    String<typename SAValue<TIndex>::Type > occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 4u);
    SEQAN_ASSERT(contains(occurrences, 0u));
    SEQAN_ASSERT(contains(occurrences, 2u));
    SEQAN_ASSERT(contains(occurrences, 6u));
    SEQAN_ASSERT(contains(occurrences, 10u));


    clear(it);
    pattern = "G";
    notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);

    occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 7u);
    SEQAN_ASSERT(contains(occurrences, 1u));
    SEQAN_ASSERT(contains(occurrences, 3u));
    SEQAN_ASSERT(contains(occurrences, 7u));
    SEQAN_ASSERT(contains(occurrences, 8u));
    SEQAN_ASSERT(contains(occurrences, 11u));
    SEQAN_ASSERT(contains(occurrences, 12u));
    SEQAN_ASSERT(contains(occurrences, 13u));

    clear(it);
    pattern = "CAGG";
    notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 2u);
    SEQAN_ASSERT(contains(occurrences, 5u));
    SEQAN_ASSERT(contains(occurrences, 9u));
}

SEQAN_DEFINE_TEST(test_index_sttd64_iterateTopDownPrefLen0Files)
{
    using namespace seqan;
    typedef Index<String<Dna>, IndexSttd64<Sttd64Config<0> > > TIndex;
    
    String<Dna> txt;
    char const* const sufTreePath = "/home/aleaum/DNA/fasta_t_2";
    //char const* outDirName = "/home/aleaum/DNA/";

    String<char, MMap<> > mmapString;
    //try to load the externally stored input text
    if (!open(mmapString, sufTreePath, OPEN_RDONLY))
    {
        std::cerr << "Could not open string file!\n";
    }

    AutoSeqFormat format;
    guessFormat(mmapString, format);
    resize(txt, length(mmapString), Exact());
    assignSeq(txt, mmapString, format);
    close(mmapString);

    TIndex index(txt);
    index.inMem = false;

    //TIndex index("/home/aleaum/DNA/fasta_t_2", "/home/aleaum/DNA/");
    Iterator< TIndex, TopDown<> >::Type it(index);
    _printPartTree(getValue(container(it).trees, 0), true);

    String<Dna> pattern = "AG";
    bool notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    String<typename SAValue<TIndex>::Type > occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 4u);
    SEQAN_ASSERT(contains(occurrences, 0u));
    SEQAN_ASSERT(contains(occurrences, 2u));
    SEQAN_ASSERT(contains(occurrences, 6u));
    SEQAN_ASSERT(contains(occurrences, 10u));


    clear(it);
    pattern = "G";
    notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);

    occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 7u);
    SEQAN_ASSERT(contains(occurrences, 1u));
    SEQAN_ASSERT(contains(occurrences, 3u));
    SEQAN_ASSERT(contains(occurrences, 7u));
    SEQAN_ASSERT(contains(occurrences, 8u));
    SEQAN_ASSERT(contains(occurrences, 11u));
    SEQAN_ASSERT(contains(occurrences, 12u));
    SEQAN_ASSERT(contains(occurrences, 13u));

    clear(it);
    pattern = "CAGG";
    notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 2u);
    SEQAN_ASSERT(contains(occurrences, 5u));
    SEQAN_ASSERT(contains(occurrences, 9u));
}

SEQAN_DEFINE_TEST(test_index_sttd64_iterateTopDownLongInput)
{
    using namespace seqan;
    typedef Index<String<Dna5>, IndexSttd64<Sttd64Config<0> > > TIndex;
    TIndex index("AGGCTAGTCAGTATTAACATAGACTAGCTAGCATGACATTGCCADTGAGCTACCTATTATACTCCTCTATCTACTATATTCATCTTACTTGGACATCCTAGTCCTATAGTACTATCAGTTCATTATCCAGCATGTCAGCCCTGGACTAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGT");
    Iterator< TIndex, TopDown<> >::Type it(index);

    String<Dna> pattern = "CTACTATATTCATCTTACTTGGACATCCT";
    bool notFound = false;
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
        {
            notFound = true;
            break;
        }
        unsigned endPos = _min(repLength(it), length(pattern));
        // compare remaining edge characters with pattern
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
        {
            notFound = true;
            break;
        }
    }
    SEQAN_ASSERT_NOT(notFound);
    String<typename SAValue<TIndex>::Type > occurrences = getOccurrences(it);
    SEQAN_ASSERT_EQ(length(occurrences) , 1u);
    SEQAN_ASSERT(contains(occurrences, 70u));

//
//    clear(it);
//    pattern = "G";
//    notFound = false;
//    while (repLength(it) < length(pattern))
//    {
//        // go down edge starting with the next pattern character
//        if (!goDown(it, pattern[repLength(it)]))
//        {
//            notFound = true;
//            break;
//        }
//        unsigned endPos = _min(repLength(it), length(pattern));
//        // compare remaining edge characters with pattern
//        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
//            infix(pattern, parentRepLength(it) + 1, endPos))
//        {
//            notFound = true;
//            break;
//        }
//    }
//    SEQAN_ASSERT_NOT(notFound);
//
//    occurrences = getOccurrences(it);
//    SEQAN_ASSERT_EQ(length(occurrences) , 7u);
//    SEQAN_ASSERT(contains(occurrences, 1u));
//    SEQAN_ASSERT(contains(occurrences, 3u));
//    SEQAN_ASSERT(contains(occurrences, 7u));
//    SEQAN_ASSERT(contains(occurrences, 8u));
//    SEQAN_ASSERT(contains(occurrences, 11u));
//    SEQAN_ASSERT(contains(occurrences, 12u));
//    SEQAN_ASSERT(contains(occurrences, 13u));
//
//    clear(it);
//    pattern = "CAGG";
//    notFound = false;
//    while (repLength(it) < length(pattern))
//    {
//        // go down edge starting with the next pattern character
//        if (!goDown(it, pattern[repLength(it)]))
//        {
//            notFound = true;
//            break;
//        }
//        unsigned endPos = _min(repLength(it), length(pattern));
//        // compare remaining edge characters with pattern
//        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
//            infix(pattern, parentRepLength(it) + 1, endPos))
//        {
//            notFound = true;
//            break;
//        }
//    }
//    SEQAN_ASSERT_NOT(notFound);
//    occurrences = getOccurrences(it);
//    SEQAN_ASSERT_EQ(length(occurrences) , 2u);
//    SEQAN_ASSERT(contains(occurrences, 5u));
//    SEQAN_ASSERT(contains(occurrences, 9u));

}

}  // namespace seqan

#endif //TESTS_INDEX_TEST_INDEX_STTD64_H_
