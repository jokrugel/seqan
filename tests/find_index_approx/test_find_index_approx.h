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
// Tests for the SeqAn module find_index_approx.
// ==========================================================================

#ifndef TESTS_FIND_INDEX_APPROX_TEST_FIND_INDEX_APPROX_H_
#define TESTS_FIND_INDEX_APPROX_TEST_FIND_INDEX_APPROX_H_

#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/index2.h>
#include <seqan/find_index_approx.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// This class will store the correct matches and will be used in a set<Match>.
// These "truth tables" have been populated by using the DPSearch algorithm.
// All correctly found occurrences will be deleted, so the set should be empty in the end.
struct Match {
    typedef Position<String<char> >::Type TPosition;
    typedef Value<EditDistanceScore>::Type TScoreValue;

    TPosition position, beginPosition, endPosition;
    Value<EditDistanceScore>::Type score, beginScore;
    String<char> infix;

    Match(TPosition position2, TPosition beginPosition2, TPosition endPosition2, TScoreValue score2, TScoreValue beginScore2, String<char> const & infix2)
        : position(position2), beginPosition(beginPosition2), endPosition(endPosition2), score(score2), beginScore(beginScore2), infix(infix2) {}

    bool operator<(Match const & other) const {
        return this->beginPosition < other.beginPosition || (this->beginPosition == other.beginPosition && this->endPosition < other.endPosition);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TNeedle, typename TIndex, typename TPatternSpec>
struct IndexPattern {
    typedef TPatternSpec                                    Type;
};

template <typename TNeedle, typename TIndex, typename TSpec, typename TScore, typename TVerifyPatternSpec>
struct IndexPattern<TNeedle, TIndex, Partitioning<TSpec, TScore, Default, TVerifyPatternSpec> > {
    typedef typename DefaultIndexPattern<TNeedle, TIndex>::Type TPiecePatternSpec;
    typedef Partitioning<TSpec, TScore, TPiecePatternSpec, TVerifyPatternSpec> Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TIndex>
inline void
configureIndex(TIndex & index) {
    ignoreUnusedVariableWarning(index);
    // do nothing
}

template <typename TText>
inline void
configureIndex(Index<TText, IndexQGram<UngappedShape<3> > > & index) {
    setStepSize(index, 2u);
}

// The actual test function. It will get instantiated with different combinations of Pattern, Finder and Index specializations.
template <typename TIndexSpec, typename TFinderSpec, typename TPatternSpecTmp>
void testIndexApprox() {
    typedef std::set<Match>                                 TMatches;
    typedef TMatches::iterator                              TMatchesIter;

    // First test the basic functionality: find(), position()
    // The correct positions of the occurrences are: 8, 9, 9, 9, 10
    // The index based approximate search algorithms can return the occurrences in arbitrary order.
    // The same position can be reported multiple times, if there are several beginPositions for one endPosition.
    {
        typedef String<char>                                THaystack;
        typedef String<char>                                TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;
        typedef Finder<TIndex, TFinderSpec>                 TFinder;
        typedef typename IndexPattern<TNeedle, TIndex, TPatternSpecTmp>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;

        THaystack txt("any_annealing");
        TNeedle nl("annual");
        TIndex hstk(txt);
        configureIndex(hstk);
        
        TFinder fd(hstk);
        TPattern pt(nl, -2);
        unsigned int pos1, pos2, pos3;
        
        StringSet<TNeedle> needles;
        appendValue(needles, nl);
        TPattern tempPattern;
        preparePatterns(fd, tempPattern, needles, 2);

        SEQAN_ASSERT(find(fd, pt));
        pos1 = position(fd);
        SEQAN_ASSERT(pos1 == 8u || pos1 == 9u || pos1 == 10u);

        SEQAN_ASSERT(find(fd, pt));
        pos2 = position(fd);
        SEQAN_ASSERT(pos2 == 8u || pos2 == 9u || pos2 == 10u);

        SEQAN_ASSERT(find(fd, pt));
        pos3 = position(fd);
        SEQAN_ASSERT(pos3 == 8u || pos3 == 9u || pos3 == 10u);
        
        // find() is allowed to return true 2 more times for the remaining positions.
        SEQAN_ASSERT(!find(fd, pt) || position(fd) == 8u || position(fd) == 9u || position(fd) == 10u);
        SEQAN_ASSERT(!find(fd, pt) || position(fd) == 8u || position(fd) == 9u || position(fd) == 10u);
        SEQAN_ASSERT(!find(fd, pt));
    }

    // Needle doesn't match anywhere.
    {
        typedef String<Dna>                                 THaystack;
        typedef String<Dna>                                 TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;
        typedef Finder<TIndex, TFinderSpec>                 TFinder;
        typedef typename IndexPattern<TNeedle, TIndex, TPatternSpecTmp>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;

        THaystack txt("ACGTACGTACGTACGTACGT");
        TNeedle nl("GTGTGTGTGTGTG");
        TIndex hstk(txt);
        TFinder fd(hstk);
        TPattern pt(nl, -2);

        StringSet<TNeedle> needles;
        appendValue(needles, nl);
        TPattern tempPattern;
        preparePatterns(fd, tempPattern, needles, 2);

        SEQAN_ASSERT_NOT(find(fd, pt));
    }

    // Exact search.
    // The correct positions of the occurrences are: 2, 5, 8.
    {
        typedef String<Dna>                                 THaystack;
        typedef String<Dna>                                 TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;
        typedef Finder<TIndex, TFinderSpec>                 TFinder;
        typedef typename IndexPattern<TNeedle, TIndex, TPatternSpecTmp>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;

        THaystack txt("ACGACGACG");
        TNeedle nl("ACG");
        TIndex hstk(txt);
        configureIndex(hstk);
        TFinder fd(hstk);
        TPattern pt(nl, 0);
        unsigned int pos1, pos2, pos3;
        
        StringSet<TNeedle> needles;
        appendValue(needles, nl);
        TPattern tempPattern;
        preparePatterns(fd, tempPattern, needles, 0);

        SEQAN_ASSERT(find(fd, pt));
        pos1 = position(fd);
        SEQAN_ASSERT(pos1 == 2u || pos1 == 5u || pos1 == 8u);
        
        SEQAN_ASSERT(find(fd, pt));
        pos2 = position(fd);
        SEQAN_ASSERT(pos2 == 2u || pos2 == 5u || pos2 == 8u);
        SEQAN_ASSERT_NEQ(pos2, pos1);

        SEQAN_ASSERT(find(fd, pt));
        pos3 = position(fd);
        SEQAN_ASSERT(pos3 == 2u || pos3 == 5u || pos3 == 8u);
        SEQAN_ASSERT_NEQ(pos3, pos1);
        SEQAN_ASSERT_NEQ(pos3, pos2);

        SEQAN_ASSERT_NOT(find(fd, pt));
    }

    // Very short haystack and TNeedle != TText
    {
        typedef String<Dna>                                 THaystack;
        typedef String<Rna>                                 TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;
        typedef Finder<TIndex, TFinderSpec>                 TFinder;
        typedef Finder<TIndex>                              TFinderExact;
        typedef typename IndexPattern<TNeedle, TIndex, TPatternSpecTmp>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;
        typedef typename DefaultIndexPattern<TNeedle, TIndex>::Type TPatternExactSpec;
        typedef Pattern<TNeedle, TPatternExactSpec>         TPatternExact;
        
        THaystack txt("GTGTG");
        TNeedle nl("GAGAGAGAGAGAG");
        TIndex hstk(txt);
        configureIndex(hstk);
        TFinder fd(hstk);
        TPattern pt(nl, -2);

        //TODO(krugel) Does not work because Finder<FinderPartitioning> depends on TNeedle = TText
        //StringSet<TNeedle> needles;
        //appendValue(needles, nl);
        //TPattern tempPattern;
        //preparePatterns(fd, tempPattern, needles, 2);

        // First try exact search
        TFinderExact fd2(hstk);
        TPatternExact pt2(nl);

        SEQAN_ASSERT_NOT(find(fd2, pt2));
        
        SEQAN_ASSERT_NOT(find(fd, pt));
    }

    // Test with long needles and a Dna Alphabet of all find-related functions:
    // find(), findBegin(), position(), beginPosition(), endPosition(), infix()
    {
        typedef String<Dna>                                 THaystack;
        typedef String<Dna>                                 TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;
        typedef Finder<TIndex, TFinderSpec>                 TFinder;
        typedef typename IndexPattern<TNeedle, TIndex, TPatternSpecTmp>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;

        THaystack txt("taaaataaaatacaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaat");
        TNeedle nl("taaaataaaatacaataaaataaaatataataaaataaaataaaat");
        TIndex hstk(txt);
        configureIndex(hstk);
        TFinder fd(hstk);
        TPattern pt(nl, -2);

        StringSet<TNeedle> needles;
        appendValue(needles, nl);
        TPattern tempPattern;
        preparePatterns(fd, tempPattern, needles, 2);

        // The results are compared to a truth table
        TMatches truth;
        truth.insert(Match(44u, 0u, 45u, -2, -2, "taaaataaaatacaataaaataaaataaaataaaataaaataaaa"));
        truth.insert(Match(45u, 1u, 46u, -1, -2, "aaaataaaatacaataaaataaaataaaataaaataaaataaaat"));
        truth.insert(Match(45u, 0u, 46u, -1, -1, "taaaataaaatacaataaaataaaataaaataaaataaaataaaat"));
        truth.insert(Match(46u, 0u, 47u, -2, -2, "taaaataaaatacaataaaataaaataaaataaaataaaataaaata"));
        truth.insert(Match(60u, 15u, 61u, -2, -2, "taaaataaaataaaataaaataaaataaaataaaataaaataaaat"));
        truth.insert(Match(65u, 20u, 66u, -2, -2, "taaaataaaataaaataaaataaaataaaataaaataaaataaaat"));
        truth.insert(Match(70u, 25u, 71u, -2, -2, "taaaataaaataaaataaaataaaataaaataaaataaaataaaat"));
        
        TMatchesIter it;
        while (find(fd, pt)) {
            while (findBegin(fd, pt)) {
                it = truth.find(Match(position(fd), beginPosition(fd), endPosition(fd), getScore(pt), getBeginScore(pt), infix(fd)));
                SEQAN_ASSERT_MSG(it != truth.end(), "Found a false positive match or the same match twice!");
                // DPBacktracking already finds complete infixes (and not only endPositions) and therefore returns the score of the infix here.
                SEQAN_ASSERT(getScore(pt) == it->score || getScore(pt) == it->beginScore);
                SEQAN_ASSERT_EQ(getBeginScore(pt), it->beginScore);
                SEQAN_ASSERT_EQ(infix(fd),         it->infix);
                truth.erase(it);
            }
        }
        SEQAN_ASSERT_EQ_MSG(truth.size(), 0u, "Some matches haven't been found.");
    }
    
    // Test with high search tolerance.
    {
        typedef String<char>                                THaystack;
        typedef String<char>                                TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;
        typedef Finder<TIndex, TFinderSpec>                 TFinder;
        typedef typename IndexPattern<TNeedle, TIndex, TPatternSpecTmp>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;

        THaystack txt("1234567890abcdefghijklXYZ1234567890abcdefghijklXYZ1234567890abcdefghijklXYZ1234567890");
        TNeedle nl("aPcQeRgSijk");
        TIndex hstk(txt);
        configureIndex(hstk);
        TFinder fd(hstk);
        TPattern pt(nl, -4);

        StringSet<TNeedle> needles;
        appendValue(needles, nl);
        TPattern tempPattern;
        preparePatterns(fd, tempPattern, needles, 4);

        // The results are compared to a truth table
        TMatches truth;
        truth.insert(Match(20u, 10u, 21u, -4, -4, "abcdefghijk"));
        truth.insert(Match(45u, 35u, 46u, -4, -4, "abcdefghijk"));
        truth.insert(Match(70u, 60u, 71u, -4, -4, "abcdefghijk"));
        
        TMatchesIter it;
        while (find(fd, pt)) {
            while (findBegin(fd, pt)) {
                it = truth.find(Match(position(fd), beginPosition(fd), endPosition(fd), getScore(pt), getBeginScore(pt), infix(fd)));
                SEQAN_ASSERT_MSG(it != truth.end(), "Found a false positive match or the same match twice!");
                // DPBacktracking already finds complete infixes (and not only endPositions) and therefore returns the score of the infix here.
                SEQAN_ASSERT(getScore(pt) == it->score || getScore(pt) == it->beginScore);
                SEQAN_ASSERT_EQ(getBeginScore(pt), it->beginScore);
                SEQAN_ASSERT_EQ(infix(fd),         it->infix);
                truth.erase(it);
            }
        }
        // q-sample cannot work correctly because needles are smaller than Q + stepSize - 1
        SEQAN_ASSERT_EQ_MSG(truth.size(), 0u, "Some matches haven't been found.");
    }

    // Test of Metafunctions, setHaystack(), setNeedle(), setScoreLimit(), needle(), haystack()
    {
        typedef String<char>                                THaystack;
        typedef String<char>                                TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;
        typedef Finder<TIndex, TFinderSpec>                 TFinder;
        typedef typename IndexPattern<TNeedle, TIndex, TPatternSpecTmp>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;

        THaystack txtEmpty("");
        THaystack txt("Dies ist der Haystack des Tests. Ja, das ist er wirklich!");
        TNeedle nl("des");
        TIndex hstkEmpty(txtEmpty);
        TIndex hstk(txt);
        configureIndex(hstk);
        TFinder fd(hstk); // fd(hstkEmpty);
        TPattern pt;
        
        setHost(fd, hstk);
        setNeedle(pt, nl);
        setScoreLimit(pt, -1);
        
        typename Haystack<TFinder>::Type &  hstk2 = haystack(fd);
        typename Host<TFinder>::Type &      hstk3 = host(fd);
        typename Container<TFinder>::Type & hstk4 = host(fd);
        typename Needle<TPattern>::Type     nl2 = needle(pt);
        typename Host<TPattern>::Type       nl3 = host(pt);
        
        SEQAN_ASSERT_EQ(indexText(hstk2), txt);
        SEQAN_ASSERT_EQ(indexText(hstk3), txt);
        SEQAN_ASSERT_EQ(indexText(hstk4), txt);
        SEQAN_ASSERT_EQ(nl2, nl);
        SEQAN_ASSERT_EQ(nl3, nl);
        SEQAN_ASSERT_EQ(needle(reinterpret_cast<TPattern const &>(pt)), nl);

        StringSet<TNeedle> needles;
        appendValue(needles, nl);
        TPattern tempPattern;
        preparePatterns(fd, tempPattern, needles, 1);

        // The results are compared to a truth table
        TMatches truth;
        truth.insert(Match(3u, 2u, 4u, -1, -1, "es"));
        truth.insert(Match(3u, 1u, 4u, -1, -1, "ies"));
        truth.insert(Match(10u, 9u, 11u, -1, -1, "de"));
        truth.insert(Match(11u, 9u, 12u, -1, -1, "der"));
        truth.insert(Match(23u, 22u, 24u, -1, -1, "de"));
        truth.insert(Match(24u, 23u, 25u, 0, -1, "es"));
        truth.insert(Match(24u, 22u, 25u, 0, 0, "des"));
        truth.insert(Match(24u, 21u, 25u, 0, -1, " des"));
        truth.insert(Match(25u, 22u, 26u, -1, -1, "des "));
        truth.insert(Match(28u, 27u, 29u, -1, -1, "es"));
        truth.insert(Match(28u, 26u, 29u, -1, -1, "Tes"));
        truth.insert(Match(39u, 37u, 40u, -1, -1, "das"));

        TMatchesIter it;
        while (find(fd, pt)) {
            while (findBegin(fd, pt)) {
                it = truth.find(Match(position(fd), beginPosition(fd), endPosition(fd), getScore(pt), getBeginScore(pt), infix(fd)));
                SEQAN_ASSERT_MSG(it != truth.end(), "Found a false positive match or the same match twice!");
                // DPBacktracking already finds complete infixes (and not only endPositions) and therefore returns the score of the infix here.
                SEQAN_ASSERT(getScore(pt) == it->score || getScore(pt) == it->beginScore);
                SEQAN_ASSERT_EQ(getBeginScore(pt), it->beginScore);
                SEQAN_ASSERT_EQ(infix(fd),         it->infix);
                truth.erase(it);
            }
        }
        SEQAN_ASSERT_EQ_MSG(truth.size(), 0u, "Some matches haven't been found.");
    }

    // Test of clear() and recycling the finder and pattern.
    {
        typedef String<char>                                THaystack;
        typedef String<char>                                TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;
        typedef Finder<TIndex, TFinderSpec>                 TFinder;
        typedef typename IndexPattern<TNeedle, TIndex, TPatternSpecTmp>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;

        THaystack txt = "123XXXabaXXX45aba123";
        TNeedle nl = "XXXaba";
        TIndex hstk(txt);
        configureIndex(hstk);
        TFinder fd(hstk);
        TPattern pt(nl, -2);

        StringSet<TNeedle> needles;
        appendValue(needles, nl);
        TPattern tempPattern;
        preparePatterns(fd, tempPattern, needles, 2);

        // The results are compared to a truth table
        TMatches truth, truth2;
        truth.insert(Match(6u, 3u, 7u, -2, -2, "XXXa"));
        truth.insert(Match(7u, 4u, 8u, -1, -2, "XXab"));
        truth.insert(Match(7u, 3u, 8u, -1, -1, "XXXab"));
        truth.insert(Match(7u, 2u, 8u, -1, -2, "3XXXab"));
        truth.insert(Match(8u, 5u, 9u, 0, -2, "Xaba"));
        truth.insert(Match(8u, 4u, 9u, 0, -1, "XXaba"));
        truth.insert(Match(8u, 3u, 9u, 0, 0, "XXXaba"));
        truth.insert(Match(8u, 2u, 9u, 0, -1, "3XXXaba"));
        truth.insert(Match(8u, 1u, 9u, 0, -2, "23XXXaba"));
        truth.insert(Match(9u, 4u, 10u, -1, -2, "XXabaX"));
        truth.insert(Match(9u, 3u, 10u, -1, -1, "XXXabaX"));
        truth.insert(Match(9u, 2u, 10u, -1, -2, "3XXXabaX"));
        truth.insert(Match(10u, 3u, 11u, -2, -2, "XXXabaXX"));
        truth.insert(Match(14u, 9u, 15u, -2, -2, "XXX45a"));
        truth.insert(Match(16u, 11u, 17u, -2, -2, "X45aba"));
        truth.insert(Match(16u, 10u, 17u, -2, -2, "XX45aba"));
        truth.insert(Match(16u, 9u, 17u, -2, -2, "XXX45aba"));
        truth2 = truth; // copy to use it again after clear()
        
        TMatchesIter it;
        while (find(fd, pt)) {
            while (findBegin(fd, pt)) {
                it = truth.find(Match(position(fd), beginPosition(fd), endPosition(fd), getScore(pt), getBeginScore(pt), infix(fd)));
                SEQAN_ASSERT_MSG(it != truth.end(), "Found a false positive match or the same match twice!");
                // DPBacktracking already finds complete infixes (and not only endPositions) and therefore returns the score of the infix here.
                SEQAN_ASSERT(getScore(pt) == it->score || getScore(pt) == it->beginScore);
                SEQAN_ASSERT_EQ(getBeginScore(pt), it->beginScore);
                SEQAN_ASSERT_EQ(infix(fd),         it->infix);
                truth.erase(it);
            }
        }
        SEQAN_ASSERT_EQ_MSG(truth.size(), 0u, "Some matches haven't been found.");
        
        clear(fd);
        SEQAN_ASSERT_EQ(indexText(hstk), "123XXXabaXXX45aba123");

        while (find(fd, pt)) {
            while (findBegin(fd, pt)) {
                it = truth2.find(Match(position(fd), beginPosition(fd), endPosition(fd), getScore(pt), getBeginScore(pt), infix(fd)));
                SEQAN_ASSERT_MSG(it != truth2.end(), "Found a false positive match or the same match twice!");
                // DPBacktracking already finds complete infixes (and not only endPositions) and therefore returns the score of the infix here.
                SEQAN_ASSERT(getScore(pt) == it->score || getScore(pt) == it->beginScore);
                SEQAN_ASSERT_EQ(getBeginScore(pt), it->beginScore);
                SEQAN_ASSERT_EQ(infix(fd),         it->infix);
                truth2.erase(it);
            }
        }

        SEQAN_ASSERT_EQ_MSG(truth.size(), 0u, "Some matches haven't been found.");
    }

    // Search in a fasta file, also with clear() and setScoreLimit().
    {
        typedef String<char>                                THaystack;
        typedef String<char>                                TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;
        typedef Finder<TIndex, TFinderSpec>                 TFinder;
        typedef typename IndexPattern<TNeedle, TIndex, TPatternSpecTmp>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;

        String<char> testFile = getAbsolutePath("/tests/find_index_approx/fasta-small.fasta");
        SeqFileIn seqIO(toCString(testFile));
        String<char> id;
        THaystack txt;
        readRecord(id, txt, seqIO);

        //String<char, FileReader<Fasta> > fileTxt(testFile);
        //THaystack txt(fileTxt);

        TNeedle nl = "ACGTACGT";
        TIndex hstk(txt);
        configureIndex(hstk);
        TFinder fd(hstk);
        TPattern pt(nl, -1);

        StringSet<TNeedle> needles;
        appendValue(needles, nl);
        TPattern tempPattern;
        preparePatterns(fd, tempPattern, needles, 1);

        // The results are compared to a truth table
        TMatches truth;
        truth.insert(Match(970u, 964u, 971u, -1, -1, "ACTACGT"));
        truth.insert(Match(2158u, 2151u, 2159u, -1, -1, "ACTTACGT"));
        
        TMatchesIter it;
        while (find(fd, pt)) {
            while (findBegin(fd, pt)) {
                it = truth.find(Match(position(fd), beginPosition(fd), endPosition(fd), getScore(pt), getBeginScore(pt), infix(fd)));
                SEQAN_ASSERT_MSG(it != truth.end(), "Found a false positive match or the same match twice!");
                // DPBacktracking already finds complete infixes (and not only endPositions) and therefore returns the score of the infix here.
                SEQAN_ASSERT(getScore(pt) == it->score || getScore(pt) == it->beginScore);
                SEQAN_ASSERT_EQ(getBeginScore(pt), it->beginScore);
                SEQAN_ASSERT_EQ(infix(fd),         it->infix);
                truth.erase(it);
            }
        }

        // q-sample cannot work correctly because needles are smaller than Q + stepSize - 1
        SEQAN_ASSERT_EQ_MSG(truth.size(), 0u, "Some matches haven't been found.");

        clear(fd);
        
        TNeedle nl2 = "ACGTACGTTAC";
        setNeedle(pt, nl2);
        setScoreLimit(pt, -2);
        
        clear(needles);
        appendValue(needles, nl);
        clear(tempPattern);
        preparePatterns(fd, tempPattern, needles, 2);

        // The results are compared to a truth table
        TMatches truth2;
        truth2.insert(Match(696u, 687u, 697u, -2, -2, "ACTTACTTAC"));
        truth2.insert(Match(710u, 701u, 711u, -2, -2, "ACTTACTTAC"));
        truth2.insert(Match(1386u, 1377u, 1387u, -2, -2, "ACTACCTTAC"));

        //TMatchesIter it; // Already declared above
        while (find(fd, pt)) {
            while (findBegin(fd, pt)) {
                it = truth2.find(Match(position(fd), beginPosition(fd), endPosition(fd), getScore(pt), getBeginScore(pt), infix(fd)));
                SEQAN_ASSERT_MSG(it != truth2.end(), "Found a false positive match or the same match twice!");
                // DPBacktracking already finds complete infixes (and not only endPositions) and therefore returns the score of the infix here.
                SEQAN_ASSERT(getScore(pt) == it->score || getScore(pt) == it->beginScore);
                SEQAN_ASSERT_EQ(getBeginScore(pt), it->beginScore);
                SEQAN_ASSERT_EQ(infix(fd),         it->infix);
                truth2.erase(it);
            }
        }
        SEQAN_ASSERT_EQ_MSG(truth.size(), 0u, "Some matches haven't been found.");
    }
}

// Test the special case of giving a predefined number of pieces for the intermediate partitioning pattern.
template <typename TIndexSpec, typename TFinderSpec, typename TPatternSpecTmp>
void testIndexApproxNumberOfPieces() {
    typedef std::set<Match>                                 TMatches;
    typedef TMatches::iterator                              TMatchesIter;

    // Test with high search tolerance.
    // (Just copied from above, except for the third parameter in the pattern constructor.)
    {
        typedef String<char>                                THaystack;
        typedef String<char>                                TNeedle;
        typedef Index<THaystack, TIndexSpec>                TIndex;
        typedef Finder<TIndex, TFinderSpec>                 TFinder;
        typedef typename IndexPattern<TNeedle, TIndex, TPatternSpecTmp>::Type TPatternSpec;
        typedef Pattern<TNeedle, TPatternSpec>              TPattern;

        THaystack txt("1234567890abcdefghijklXYZ1234567890abcdefghijklXYZ1234567890abcdefghijklXYZ1234567890");
        TNeedle nl("aPcQeRgSijk");
        TIndex hstk(txt);
        configureIndex(hstk);
        TFinder fd(hstk);
        TPattern pt(nl, -4, EditDistanceScore(), 2);  // Use two pieces and search each with tolerance 2.

        // The results are compared to a truth table
        TMatches truth;
        truth.insert(Match(20u, 10u, 21u, -4, -4, "abcdefghijk"));
        truth.insert(Match(45u, 35u, 46u, -4, -4, "abcdefghijk"));
        truth.insert(Match(70u, 60u, 71u, -4, -4, "abcdefghijk"));
        
        TMatchesIter it;
        while (find(fd, pt)) {
            while (findBegin(fd, pt)) {
                it = truth.find(Match(position(fd), beginPosition(fd), endPosition(fd), getScore(pt), getBeginScore(pt), infix(fd)));
                SEQAN_ASSERT_MSG(it != truth.end(), "Found a false positive match or the same match twice!");
                // DPBacktracking already finds complete infixes (and not only endPositions) and therefore returns the score of the infix here.
                SEQAN_ASSERT(getScore(pt) == it->score || getScore(pt) == it->beginScore);
                SEQAN_ASSERT_EQ(getBeginScore(pt), it->beginScore);
                SEQAN_ASSERT_EQ(infix(fd),         it->infix);
                truth.erase(it);
            }
        }
        SEQAN_ASSERT_EQ_MSG(truth.size(), 0u, "Some matches haven't been found.");
    }
}

// Dynamic programming backtracking

SEQAN_DEFINE_TEST(test_find_index_approx_dpbacktracking_esa) {
    testIndexApprox<
        IndexEsa<>,
        DPBacktracking<>,
        DPBacktracking<>
    >();
}

SEQAN_DEFINE_TEST(test_find_index_approx_dpbacktracking_wotd) {
    testIndexApprox<
        IndexWotd<>,
        DPBacktracking<>,
        DPBacktracking<>
    >();
}

SEQAN_DEFINE_TEST(test_find_index_approx_dpbacktracking_sttd64) {
    testIndexApprox<
        IndexSttd64<>,
        DPBacktracking<>,
        DPBacktracking<>
    >();
}

// Partitioning into exact search

SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intoexactsearch_esa) {
    testIndexApprox<
        IndexEsa<>,
        FinderPartitioning<>,
        Partitioning<IntoExactSearch>
    >();
}

SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intoexactsearch_wotd) {
    testIndexApprox<
        IndexWotd<>,
        FinderPartitioning<>,
        Partitioning<IntoExactSearch>
    >();
}

SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intoexactsearch_fm) {
    testIndexApprox<
        FMIndex<>,
        FinderPartitioning<>,
        Partitioning<IntoExactSearch>
    >();
}

SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intoexactsearch_qgram) {
    testIndexApprox<
        IndexQGram<UngappedShape<2> >,
        FinderPartitioning<>,
        Partitioning<IntoExactSearch, EditDistanceScore>
    >();
}

SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intoexactsearch_qsample) {
    testIndexApprox<
        IndexQGram<UngappedShape<3> >,
        FinderPartitioning<>,
        Partitioning<IntoExactSearch, EditDistanceScore>
    >();
}

SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intoexactsearch_sttd64) {
    testIndexApprox<
        IndexSttd64<>,
        FinderPartitioning<>,
        Partitioning<IntoExactSearch>
    >();
}

SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intoexactsearch_lz) {
    testIndexApprox<
        IndexLZ<>,
        FinderPartitioning<>,
        Partitioning<IntoExactSearch, EditDistanceScore>
    >();
}

SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intoexactsearch_sadakane) {
    testIndexApprox<
        IndexSadakane<>,
        FinderPartitioning<>,
        Partitioning<IntoExactSearch, EditDistanceScore>
    >();
}

SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intoexactsearch_qgram2l) {
    testIndexApprox<
        IndexQGram2L<UngappedShape<2> >,
        FinderPartitioning<>,
        Partitioning<IntoExactSearch, EditDistanceScore>
    >();
}

SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intoexactsearch_digest) {
    // TODO(krugel) Digest Config
    const unsigned PARTITION_SIZE = 5;
    const unsigned OUTBUF_SIZE = 5;
    const unsigned INBUF_SIZE = 10;
    const unsigned TAIL_LENGTH = 3;
    const unsigned PREFIX_LENGTH = 2;

    typedef DigestConfig<PARTITION_SIZE, OUTBUF_SIZE, INBUF_SIZE, TAIL_LENGTH, PREFIX_LENGTH> TConfig;

    testIndexApprox<
        IndexDigest<TConfig>,
        FinderPartitioning<>,
        Partitioning<IntoExactSearch, EditDistanceScore>
    >();
}

SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intoexactsearch_digest_prepare) {
    // TODO(krugel) Digest Config
    const unsigned PARTITION_SIZE = 5;
    const unsigned OUTBUF_SIZE = 5;
    const unsigned INBUF_SIZE = 10;
    const unsigned TAIL_LENGTH = 3;
    const unsigned PREFIX_LENGTH = 2;

    typedef DigestConfig<PARTITION_SIZE, OUTBUF_SIZE, INBUF_SIZE, TAIL_LENGTH, PREFIX_LENGTH> TConfig;

    testIndexApprox<
        IndexDigest<TConfig>,
        FinderPartitioning<Default, Default, True>,  // DoPreparePatterns
        Partitioning<IntoExactSearch, EditDistanceScore>
    >();
}

// ...with explicit verify algorithm, here dynamic programming
SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intoexactsearch_esa_dpsearch) {
    testIndexApprox<
        IndexEsa<>,
        FinderPartitioning<>,
        Partitioning<IntoExactSearch, EditDistanceScore, Default, DPSearch<EditDistanceScore> >
    >();
}

// ...with explicit verify algorithm, here Myers (default)
SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intoexactsearch_esa_myers) {
    testIndexApprox<
        IndexEsa<>,
        FinderPartitioning<>,
        Partitioning<IntoExactSearch, EditDistanceScore, Default, Myers<FindInfix> >
    >();
}

// Intermediate partitioning

SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intermediate_esa) {
    testIndexApprox<
        IndexEsa<>,
        FinderPartitioning<DPBacktracking<> >,
        Partitioning<Intermediate, EditDistanceScore, DPBacktracking<> >
    >();
}

SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intermediate_wotd) {
    testIndexApprox<
        IndexWotd<>,
        FinderPartitioning<DPBacktracking<> >,
        Partitioning<Intermediate, EditDistanceScore, DPBacktracking<> >
    >();
}

SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intermediate_sttd64) {
    testIndexApprox<
        IndexSttd64<>,
        FinderPartitioning<DPBacktracking<> >,
        Partitioning<Intermediate, EditDistanceScore, DPBacktracking<> >
    >();
}

// ...with specified number of pieces
SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_intermediate_esa_numberofpieces) {
    testIndexApproxNumberOfPieces<
        IndexEsa<>,
        FinderPartitioning<DPBacktracking<> >,
        Partitioning<Intermediate, EditDistanceScore, DPBacktracking<> >
    >();
}

// ...with partitioning into exact search to find the pieces
SEQAN_DEFINE_TEST(test_find_index_approx_partitioning_hierarchical_esa) {
    testIndexApprox<
        IndexEsa<>,
        FinderPartitioning<FinderPartitioning<> >,
        Partitioning<PartitioningHierarchical, EditDistanceScore, Partitioning<IntoExactSearch> >
    >();
}

template <typename TFinderSpec, typename TPatternSpec, typename THaystack, typename TSegmentOrString>
void test_find_index_on_segments_Helper(THaystack &haystack, const TSegmentOrString &needle) {
    Finder<THaystack, TFinderSpec> finder(haystack);
    Pattern<TSegmentOrString, TPatternSpec> pattern(needle);

    bool didFind = false;
    while (find(finder, pattern)) {
        findBegin(finder, pattern);
        didFind = true;
    }
    SEQAN_ASSERT(didFind);

    // TODO(holtgrew): Some kind of assertion on the results.
}

template <typename TIndexSpec, typename TFinderSpec, typename TPatternSpec>
void test_find_index_on_segments_Helper2() {
    // TODO(holtgrew): Should be const.
    CharString kText = "CGATCGAT";
    Index<CharString, TIndexSpec> kHaystack(kText);
    CharString kNeedle = "GATC";

    test_find_index_on_segments_Helper<TFinderSpec, TPatternSpec>(kHaystack, kNeedle);

    Segment<CharString, PrefixSegment> myPrefix(prefix(kNeedle, 1));
    test_find_index_on_segments_Helper<TFinderSpec, TPatternSpec>(kHaystack, myPrefix);
    test_find_index_on_segments_Helper<TFinderSpec, TPatternSpec>(kHaystack, prefix(kNeedle, 1));

    Segment<CharString, InfixSegment> myInfix(infix(kNeedle, 1, 2));
    test_find_index_on_segments_Helper<TFinderSpec, TPatternSpec>(kHaystack, myInfix);
    test_find_index_on_segments_Helper<TFinderSpec, TPatternSpec>(kHaystack, infix(kNeedle, 1, 2));

    Segment<CharString, SuffixSegment> mySuffix(suffix(kNeedle, 1));
    test_find_index_on_segments_Helper<TFinderSpec, TPatternSpec>(kHaystack, mySuffix);
    test_find_index_on_segments_Helper<TFinderSpec, TPatternSpec>(kHaystack, suffix(kNeedle, 1));
}

// Test string search code on segments.  At the moment, this only
// tests whether the code compiles and runs through without any
// obvious errors.  No checks are done on the results.
SEQAN_DEFINE_TEST(test_find_index_on_segments) {
    test_find_index_on_segments_Helper2<IndexEsa<>, DPBacktracking<>, DPBacktracking<> >();
    test_find_index_on_segments_Helper2<IndexEsa<>, FinderPartitioning<>, Partitioning<> >();
    test_find_index_on_segments_Helper2<IndexEsa<>, FinderPartitioning<DPBacktracking<> >, Partitioning<Intermediate, EditDistanceScore, DPBacktracking<> > >();
    test_find_index_on_segments_Helper2<IndexEsa<>, FinderPartitioning<FinderPartitioning<> >, Partitioning<PartitioningHierarchical, EditDistanceScore, Partitioning<IntoExactSearch> > >();
}

template <typename TPatternSpec>
void test_pattern_copycon() {
    typedef Pattern<CharString, TPatternSpec> TPattern;
    TPattern p1("Some needle");
    TPattern p2(p1);
    TPattern const p3(p2);
    TPattern const p4(p3);
}

template <typename TPatternSpec>
void test_pattern_assign() {
    typedef Pattern<CharString, TPatternSpec> TPattern;
    TPattern p1("Some needle");
    TPattern p2;
    TPattern const p3(p1);
    TPattern p4;
    p2 = p1;
    p4 = p3;
}

#ifdef SEQAN_CXX11_STANDARD

template <typename TPatternSpec>
void test_pattern_movecon() {
    typedef Pattern<CharString, TPatternSpec> TPattern;
    TPattern p1("Some needle");
    TPattern p2(std::move(p1));
    TPattern const p3(std::move(p2));
    TPattern const p4(std::move(p3));
}

template <typename TPatternSpec>
void test_pattern_moveassign() {
    typedef Pattern<CharString, TPatternSpec> TPattern;
    TPattern p1("Some needle");
    TPattern p2;
    TPattern const p3(p1);
    TPattern p4;
    p2 = std::move(p1);
    p4 = std::move(p3);
}

template <typename TPatternSpec>
void test_pattern_set_host()
{
    typedef Pattern<DnaString, TPatternSpec> TPattern;

    {  // Set to lvalue reference.
        TPattern p;
        DnaString ndl = "AGTGGATAGAGAT";
        setHost(p, ndl);

        SEQAN_ASSERT(&host(p) == &ndl);
    }

    {  // Set to lvalue with different alphabet causing the source to be changed.
        DnaString ndl = "AGTGGATAGAGAT";
        TPattern p;
        setHost(p, ndl);
        SEQAN_ASSERT_EQ(&host(p), &ndl);

        CharString ndl2 = "AT CARLS";
        setHost(p, ndl2);
        SEQAN_ASSERT(&host(p) == &ndl);
        SEQAN_ASSERT_EQ(ndl, "ATACAAAA");
        SEQAN_ASSERT_EQ(ndl2, "AT CARLS");
    }

    {  // Set to const lvalue reference with different alphabet causing the source to be changed.
        DnaString const ndl = "AGTGGATAGAGAT";
        TPattern p;
        setHost(p, ndl);
        SEQAN_ASSERT(&host(p) != &ndl);
        SEQAN_ASSERT_EQ(host(p), ndl);

        CharString ndl2 = "AT CARLS";
        setHost(p, ndl2);
        SEQAN_ASSERT(&host(p) != &ndl);
        SEQAN_ASSERT_EQ(ndl, "AGTGGATAGAGAT");
        SEQAN_ASSERT_EQ(ndl2, "AT CARLS");
        SEQAN_ASSERT_EQ(host(p), "ATACAAAA");
    }

    {  // Set with rvalue reference.
        DnaString ndl = "AGTGGATAGAGAT";
        TPattern p;
        setHost(p, std::move(ndl));
        SEQAN_ASSERT(&host(p) != &ndl);
        SEQAN_ASSERT_EQ(host(p), "AGTGGATAGAGAT");
    }

    {  // Set to rvalue reference after setting holder dependent.
        DnaString ndl = "AGTGGATAGAGAT";
        TPattern p;
        setHost(p, ndl);
        SEQAN_ASSERT(&host(p) == &ndl);

        CharString ndl2 = "AT CARLS";
        setHost(p, std::move(ndl2));
        SEQAN_ASSERT(&host(p) == &ndl);
        SEQAN_ASSERT_EQ(ndl, "ATACAAAA");
        SEQAN_ASSERT_EQ(host(p), "ATACAAAA");
    }

    {  // Set to rvalue reference before setting holder dependent.
        DnaString ndl = "AGTGGATAGAGAT";
        TPattern p;
        setHost(p, std::move(ndl));
        SEQAN_ASSERT(&host(p) != &ndl);
        SEQAN_ASSERT_EQ(host(p), "AGTGGATAGAGAT");

        CharString ndl2 = "AT CARLS";
        setHost(p, ndl2);
        SEQAN_ASSERT_EQ(host(p), "ATACAAAA");
    }

    {  // Set to rvalue with different alphabet.
        TPattern p;
        setHost(p, "AT CARLS");
        SEQAN_ASSERT_EQ(host(p), "ATACAAAA");
    }
}

#endif  // SEQAN_CXX11_STANDARD

SEQAN_DEFINE_TEST(test_pattern_copycon) {
    // Test whether the needle is preserved in copying a pattern.
    // See http://trac.mi.fu-berlin.de/seqan/ticket/318
    test_pattern_copycon<DPBacktracking<> >();
    test_pattern_copycon<Partitioning<IntoExactSearch> >();
    test_pattern_copycon<Partitioning<Intermediate, EditDistanceScore, DPBacktracking<> > >();
    test_pattern_copycon<Partitioning<PartitioningHierarchical, EditDistanceScore, Partitioning<IntoExactSearch> > >();
}

SEQAN_DEFINE_TEST(test_pattern_assign) {
    // Test whether the needle is preserved in assigning a pattern.
    // See http://trac.mi.fu-berlin.de/seqan/ticket/318
    test_pattern_assign<DPBacktracking<> >();
    test_pattern_assign<Partitioning<IntoExactSearch> >();
    test_pattern_assign<Partitioning<Intermediate, EditDistanceScore, DPBacktracking<> > >();
    test_pattern_assign<Partitioning<PartitioningHierarchical, EditDistanceScore, Partitioning<IntoExactSearch> > >();
}

#ifdef SEQAN_CXX11_STANDARD
SEQAN_DEFINE_TEST(test_pattern_movecon) {
    test_pattern_movecon<DPBacktracking<> >();
    test_pattern_movecon<Partitioning<IntoExactSearch> >();
    test_pattern_movecon<Partitioning<Intermediate, EditDistanceScore, DPBacktracking<> > >();
    test_pattern_movecon<Partitioning<PartitioningHierarchical, EditDistanceScore, Partitioning<IntoExactSearch> > >();
}

SEQAN_DEFINE_TEST(test_pattern_moveassign) {
    test_pattern_moveassign<DPBacktracking<> >();
    test_pattern_moveassign<Partitioning<IntoExactSearch> >();
    test_pattern_moveassign<Partitioning<Intermediate, EditDistanceScore, DPBacktracking<> > >();
    test_pattern_moveassign<Partitioning<PartitioningHierarchical, EditDistanceScore, Partitioning<IntoExactSearch> >>();
}

// TODO(rrahn): Should be a typed test for all pattern classes.
SEQAN_DEFINE_TEST(test_pattern_set_host) {
    test_pattern_set_host<DPBacktracking<> >();
    test_pattern_set_host<Partitioning<> >();
    test_pattern_set_host<Partitioning<Intermediate, EditDistanceScore, DPBacktracking<> > >();
    test_pattern_set_host<Partitioning<PartitioningHierarchical, EditDistanceScore, Partitioning<IntoExactSearch> >>();
}
#endif  // SEQAN_CXX11_STANDARD

}  // namespace seqan

#endif //#ifndef TESTS_FIND_INDEX_APPROX_TEST_FIND_INDEX_APPROX_H_
