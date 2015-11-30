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
// Author: Andre Dau <dau@in.tum.de>
// Author: Johannes Krugel <krugel@in.tum.de>
// ==========================================================================
// Tests for the SeqAn module index2.
// ==========================================================================

// TODO(krugel) Test pattern constructor etc. like in /tests/find/test_find.cpp

#ifndef TESTS_INDEX_TEST_INDEX_DIGEST_H_
#define TESTS_INDEX_TEST_INDEX_DIGEST_H_

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/index2.h>

#include "../index/test_index_helpers.h"

namespace seqan {

const unsigned PARTITION_SIZE = 5;
const unsigned OUTBUF_SIZE = 5;
const unsigned INBUF_SIZE = 10;
const unsigned TAIL_LENGTH = 3;
const unsigned PREFIX_LENGTH = 32;

typedef DigestConfig<PARTITION_SIZE, OUTBUF_SIZE, INBUF_SIZE, TAIL_LENGTH, PREFIX_LENGTH> TDigestConfig;

const unsigned PARTITION_SIZE_2 = 64 * 1024 * 1024; // 64 * 1024 * 1024, 
const unsigned OUTBUF_SIZE_2    = 100 * 1024; // 128 * 1024 * 1024, 
const unsigned INBUF_SIZE_2     = 10 * 1024; // 10 * 1024
const unsigned TAIL_LENGTH_2    = 1000; // 1000
const unsigned PREFIX_LENGTH_2  = 32; // 32

// TODO Longer prefix

typedef DigestConfig<PARTITION_SIZE_2, OUTBUF_SIZE_2, INBUF_SIZE_2, TAIL_LENGTH_2, PREFIX_LENGTH_2> TDigestConfigBig;


template <typename TIndex, typename TLexical, typename TSize>
inline void
getLex(TLexical & lex, TSize pos1, TSize pos2) {
    typedef typename TIndex::TDigestSuffix                  TDigestSuffix;
    
    TDigestSuffix left, right;
    left.i1 = pos1;
    right.i1 = pos2;
    _getPrefix<TIndex::PREFIX_LENGTH>(left.i2, lex.text1, left.i1);
    _getPrefix<TIndex::PREFIX_LENGTH>(right.i2, lex.text2, right.i1);
    
    compare(lex, left, right);
}

SEQAN_DEFINE_TEST(test_index_digest_bits)
{
    typedef String<char>                                    TText;
    typedef Index<TText, IndexDigest<TDigestConfig> >       TIndex;
    typedef TIndex::TPrefix                                 TPrefix;

    TText txt("mississippi");

    TPrefix pref;
    String<char> prefBits;

    _getPrefix<PREFIX_LENGTH>(pref, txt, 0u);
    _getBitRepresentation(prefBits, pref);
    // ASCII values of        'm'        'i'        's'        's'
    SEQAN_ASSERT_EQ(prefBits, "01101101" "01101001" "01110011" "01110011");

    clear(pref);
    clear(prefBits);
    _getPrefix<PREFIX_LENGTH>(pref, txt, 5u);
    _getBitRepresentation(prefBits, pref);
    // ASCII values of        's'        's'        'i'        'p'
    SEQAN_ASSERT_EQ(prefBits, "01110011" "01110011" "01101001" "01110000");

    TText sffx = suffix(txt, 5u);
    String<char> bits;
    for (DigestBitIterator<TText> itBit(sffx); !atEnd(itBit); ++itBit) {
        appendValue(bits, *itBit ? '1' : '0');
    }
    //std::cerr << prefBits << std::endl << bits << std::endl;
    SEQAN_ASSERT(isPrefix(prefix(prefBits, PREFIX_LENGTH), bits));
    SEQAN_ASSERT_EQ(bits, "01110011" "01110011" "01101001" "01110000" "01110000" "01101001");
}

SEQAN_DEFINE_TEST(test_index_digest_bits_dna)
{
    typedef String<Dna>                                     TText;
    typedef Index<TText, IndexDigest<TDigestConfig> >       TIndex;
    typedef TIndex::TPrefix                                 TPrefix;

    TText txt("ACGTTGCAACGTTGCAACGTTGCA");

    TPrefix pref;
    String<char> prefBits;

    _getPrefix<PREFIX_LENGTH>(pref, txt, 0u);
    _getBitRepresentation(prefBits, pref);
    // Binary encoding
    SEQAN_ASSERT_EQ(prefBits, "0001101111100100" "0001101111100100");

    TText sffx = suffix(txt, 5u);
    String<char> bits;
    for (DigestBitIterator<TText> itBit(sffx); !atEnd(itBit); ++itBit) {
        appendValue(bits, *itBit ? '1' : '0');
    }
    //std::cerr << prefBits << std::endl << bits << std::endl;
    SEQAN_ASSERT_EQ(bits, "1001000" "001101111100100" "0001101111100100");

}

SEQAN_DEFINE_TEST(test_index_digest_lexical)
{
    typedef String<char>                                    TText;
    typedef Index<TText, IndexDigest<TDigestConfig> >       TIndex;
    typedef typename Size<TText>::Type                      TSize;
    typedef Lexical<LexicalSuffixDigest<True, PREFIX_LENGTH, TSize, TText, TText> > TLexical;

    TText txt("mississippi");

    TLexical lex(txt, txt);
    getLex<TIndex>(lex, 0u, 5u); // mississippi < ssippi
    SEQAN_ASSERT(isLess(lex));
    SEQAN_ASSERT(!isEqual(lex) && !isPrefix(lex) && !hasPrefix(lex));

    getLex<TIndex>(lex, 5u, 0u); // ssippi < mississippi
    SEQAN_ASSERT(isGreater(lex));
    SEQAN_ASSERT(!isEqual(lex) && !isPrefix(lex) && !hasPrefix(lex));
    
    getLex<TIndex>(lex, 1u, 10u); // ississippi > i
    SEQAN_ASSERT(isGreater(lex));
    SEQAN_ASSERT(hasPrefix(lex));
    SEQAN_ASSERT(!isEqual(lex) && !isPrefix(lex));
}

// A test for strings.
SEQAN_DEFINE_TEST(test_index_digest_construct_find)
{
    typedef String<char>                                    TText;
    typedef Index<TText, IndexDigest<TDigestConfig> >       TIndex;
    typedef TIndex::TPartitions                             TPartitions;
    typedef typename Iterator<TPartitions>::Type            TPartitionsIterator;
    typedef typename Reference<TPartitionsIterator>::Type   TPartitionReference;
    
    TText txt("mississippi");

    TIndex index(txt);
    index.clearRedundantFibres = false;

    indexCreate(index, DigestPartitions());
    SEQAN_ASSERT_GEQ(length(index.partitions), length(txt) / PARTITION_SIZE);

    // The partitions should cover the whole text.
    // (i1 = start position, i2 = end position, i3 = length)
    SEQAN_ASSERT_EQ(index.partitions[0u].i1, 0u);
    SEQAN_ASSERT_EQ(index.partitions[length(index.partitions) - 1u].i2, length(txt));
    
    for (TPartitionsIterator it = begin(index.partitions); !atEnd(it, index.partitions); ++it) {
        TPartitionReference partition = value(it);
        SEQAN_ASSERT_GEQ_MSG(partition.i3, 1u, "The partition is too short.");
        SEQAN_ASSERT_LEQ_MSG(partition.i3, PARTITION_SIZE, "The partition is too long.");
        // Total length minus tail length.
        SEQAN_ASSERT_LEQ_MSG(partition.i2 - partition.i1 - partition.i3, TAIL_LENGTH, "The tail is too long.");
        //std::cerr << partition.i1 << ".." << partition.i2 << " (" << partition.i3 << ")" << std::endl;
    }
    
    indexCreate(index, DigestSuffixArrays());
    // We would have to include core/tests/index/test_index_helpers.h
    // Futhermore, the data structure acually is not a suffix array, because the entries already include the offset.
    //isSuffixArray(index.suffixArrays[0], txt);
    
    SEQAN_ASSERT_EQ(length(index.suffixArrays), length(index.partitions));
    SEQAN_ASSERT_EQ(length(index.suffixArrays[0u]), index.partitions[0u].i3);
    
    indexCreate(index, DigestSuffixTrees());
    
    SEQAN_ASSERT_GEQ(length(index.suffixTrees), 1u);
    SEQAN_ASSERT_EQ(length(index.suffixTrees), length(index.dividers));

    String<char> ndl("si");
    Finder<TIndex> finder(index);
    Pattern<String<char>, PatternDigest<TIndex> > pattern(ndl);
    
    SEQAN_ASSERT(_goFirstMatchingTree(finder, pattern));
    SEQAN_ASSERT_EQ(finder.currentTree, 3u);
    SEQAN_ASSERT(_goNextMatchingTree(finder, pattern));
    SEQAN_ASSERT_EQ(finder.currentTree, 4u);
    SEQAN_ASSERT_NOT(_goNextMatchingTree(finder, pattern));
    
    clear(finder);
    clear(pattern);
    
    SEQAN_ASSERT(find(finder, pattern));
    SEQAN_ASSERT_EQ(beginPosition(finder), 6u);
    SEQAN_ASSERT(find(finder, pattern));
    SEQAN_ASSERT_EQ(beginPosition(finder), 3u);
    SEQAN_ASSERT_NOT(find(finder, pattern));
}

// A test for strings.
SEQAN_DEFINE_TEST(test_index_digest_construct_find_2)
{
    typedef String<char>                                    TText;
    typedef Index<TText, IndexDigest<TDigestConfig> >       TIndex;
    typedef TIndex::TPartitions                             TPartitions;
    typedef typename Iterator<TPartitions>::Type            TPartitionsIterator;
    typedef typename Reference<TPartitionsIterator>::Type   TPartitionReference;
    
    TText txt("abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef");

    TIndex index(txt);
    index.clearRedundantFibres = false;

    indexCreate(index, DigestPartitions());
    SEQAN_ASSERT_GEQ(length(index.partitions), length(txt) / PARTITION_SIZE);

    // The partitions should cover the whole text.
    // (i1 = start position, i2 = end position, i3 = length)
    SEQAN_ASSERT_EQ(index.partitions[0u].i1, 0u);
    SEQAN_ASSERT_EQ(index.partitions[length(index.partitions) - 1u].i2, length(txt));
    
    for (TPartitionsIterator it = begin(index.partitions); !atEnd(it, index.partitions); ++it) {
        TPartitionReference partition = value(it);
        SEQAN_ASSERT_GEQ_MSG(partition.i3, 1u, "The partition is too short.");
        SEQAN_ASSERT_LEQ_MSG(partition.i3, PARTITION_SIZE, "The partition is too long.");
        // Total length minus tail length.
        SEQAN_ASSERT_LEQ_MSG(partition.i2 - partition.i1 - partition.i3, TAIL_LENGTH, "The tail is too long.");
        //std::cerr << partition.i1 << ".." << partition.i2 << " (" << partition.i3 << ")" << std::endl;
    }
    
    indexCreate(index, DigestSuffixArrays());
    // We would have to include core/tests/index/test_index_helpers.h
    // Futhermore, the data structure acually is not a suffix array, because the entries already include the offset.
    //isSuffixArray(index.suffixArrays[0], txt);
    
    SEQAN_ASSERT_EQ(length(index.suffixArrays), length(index.partitions));
    SEQAN_ASSERT_EQ(length(index.suffixArrays[0u]), index.partitions[0u].i3);
    
    indexCreate(index, DigestSuffixTrees());
    
    SEQAN_ASSERT_GEQ(length(index.suffixTrees), 1u);
    SEQAN_ASSERT_EQ(length(index.suffixTrees), length(index.dividers));

    String<char> ndl("abcdefghijklmnopqrstuvwxyzabcdefg");
    Finder<TIndex> finder(index);
    Pattern<String<char>, PatternDigest<TIndex> > pattern(ndl);
    
    SEQAN_ASSERT(find(finder, pattern));
    SEQAN_ASSERT_EQ(beginPosition(finder), 26u);
    SEQAN_ASSERT(find(finder, pattern));
    SEQAN_ASSERT_EQ(beginPosition(finder), 0u);
    SEQAN_ASSERT_NOT(find(finder, pattern));
}

/*
// A test for strings.
SEQAN_DEFINE_TEST(test_index_digest_construct_find_file)
{
    // TODO(krugel) Packed<> is slower?      (and doesn't work with Skew7?)
    typedef Dna                                             TAlphabet;
    typedef String<TAlphabet>                               TText;
    typedef Index<TText, IndexDigest<TDigestConfigBig> >    TIndex;
    typedef TIndex::TPartitions                             TPartitions;
    typedef typename Iterator<TPartitions>::Type            TPartitionsIterator;
    typedef typename Reference<TPartitionsIterator>::Type   TPartitionReference;
    
    String<char> testFile = SEQAN_PATH_TO_ROOT();
    append(testFile, "/tests/find_index_approx/uniform-10M.txt"); //////////////////////////////////////////////////
    
    String<TAlphabet, FileReader<Fasta> > fileTxt(testFile);
    TText txt(fileTxt);

    TIndex index(txt);
    index.clearRedundantFibres = false;

    indexCreate(index, DigestPartitions());
    SEQAN_ASSERT_GEQ(length(index.partitions), length(txt) / PARTITION_SIZE_2);

    // The partitions should cover the whole text.
    // (i1 = start position, i2 = end position, i3 = length)
    SEQAN_ASSERT_EQ(index.partitions[0u].i1, 0u);
    SEQAN_ASSERT_EQ(index.partitions[length(index.partitions) - 1u].i2, length(txt));
    
    for (TPartitionsIterator it = begin(index.partitions); !atEnd(it, index.partitions); ++it) {
        TPartitionReference partition = value(it);
        SEQAN_ASSERT_GEQ_MSG(partition.i3, 1u, "The partition is too short.");
        SEQAN_ASSERT_LEQ_MSG(partition.i3, PARTITION_SIZE_2, "The partition is too long.");
        // Total length minus tail length.
        SEQAN_ASSERT_LEQ_MSG(partition.i2 - partition.i1 - partition.i3, TAIL_LENGTH_2, "The tail is too long.");
    }

    indexCreate(index, DigestSuffixArrays(), BwtWalk<>());
    // We would have to include ../../../core/tests/index/test_index_helpers.h
    // Futhermore, the data structure acually is not a suffix array, because the entries already include the offset.
    //isSuffixArray(index.suffixArrays[0], txt);
    
    SEQAN_ASSERT_EQ(length(index.suffixArrays), length(index.partitions));
    SEQAN_ASSERT_EQ(length(index.suffixArrays[0u]), index.partitions[0u].i3);
    
    indexCreate(index, DigestSuffixTrees());

    SEQAN_ASSERT_GEQ(length(index.suffixTrees), 1u);
    SEQAN_ASSERT_EQ(length(index.suffixTrees), length(index.dividers));

    // Fasta.txt
    //String<TAlphabet> ndl("TATCATCACT");
    //Finder<TIndex> finder(index);
    //Pattern<String<TAlphabet>, PatternDigest<TIndex> > pattern(ndl);

    //SEQAN_ASSERT(find(finder, pattern));
    //SEQAN_ASSERT_EQ(beginPosition(finder), 10834u);
    //SEQAN_ASSERT(find(finder, pattern));
    //SEQAN_ASSERT_EQ(beginPosition(finder), 152533u);
    //SEQAN_ASSERT(find(finder, pattern));
    //SEQAN_ASSERT_EQ(beginPosition(finder), 2330u);
    //SEQAN_ASSERT(find(finder, pattern));
    //SEQAN_ASSERT_EQ(beginPosition(finder), 150781u);
    //SEQAN_ASSERT_NOT(find(finder, pattern));

    // uniform-10M
    String<TAlphabet> ndl("AGTGCAGGTT");
    Finder<TIndex> finder(index);
    Pattern<String<TAlphabet>, PatternDigest<TIndex> > pattern(ndl);

    SEQAN_ASSERT(find(finder, pattern));
    SEQAN_ASSERT(find(finder, pattern));
    SEQAN_ASSERT(find(finder, pattern));
    //SEQAN_ASSERT(!find(finder, pattern));
}
*/
/*
SEQAN_DEFINE_TEST(test_index_digest_open_save)
{
    {
		typedef String<unsigned int>                        TPositions;
		typedef IndexDigest<TDigestConfig> 					TIndexSpec;

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
            SEQAN_ASSERT(configureSave(index0, toCString(filename)));
            
            TFinder finder0(index0);
            TNeedle ndl0("ist");	
            TPattern pattern0(ndl0);
            find(finder0, pattern0);
            
            // Save the index to disk.
            SEQAN_ASSERT(save(index0, toCString(filename)));
        }
        

        // Reopen the index and check if it still works.
        TIndex index;
        SEQAN_ASSERT(open(index, toCString(filename)));
        
        std::cout << "st: " << length(index.suffixTrees) << std::endl;
        
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
        // TODO(krugel) This fails with q-sample index
        SEQAN_ASSERT_EQ(length(pos), 2u);
        SEQAN_ASSERT_EQ(pos[0], 5u);
        SEQAN_ASSERT_EQ(pos[1], 31u);
    }
}
*/

}  // namespace seqan

#endif  // TESTS_INDEX_TEST_INDEX_DIGEST_H_
