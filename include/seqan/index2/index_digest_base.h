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
// Author: Andre Dau <dau@in.tum.de>
// Author: Johannes Krugel <krugel@in.tum.de>
// ==========================================================================
// Suffix forest data structure and construction algorithm optimized for
// secondary memory (DiGeST - Disk based Genomic Suffix Tree).
// 
// [BSTU 2008] M. Barsky, U. Stege, A. Thomo, C. Upton
// "A new method for indexing genomes using on-disk suffix trees"
// 17th ACM Conference on Information and Knowledge Management (CIKM'08),
// ACM, 2008, 649-658, http://dx.doi.org/10.1145/1458082.1458170
// ==========================================================================

#ifndef SEQAN_INDEX_INDEX_DIGEST_BASE_H_
#define SEQAN_INDEX_INDEX_DIGEST_BASE_H_

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

struct DigestDividers_;
struct DigestSuffixTrees_;
struct DigestSuffixArrays_;
struct DigestPartitions_;

typedef FibreText                       DigestText;
typedef FibreRawText                    DigestRawText;
typedef Tag<DigestPartitions_> const    DigestPartitions;
typedef Tag<DigestSuffixArrays_> const  DigestSuffixArrays;
typedef Tag<DigestDividers_> const      DigestDividers;
typedef Tag<DigestSuffixTrees_> const   DigestSuffixTrees;

// ----------------------------------------------------------------------------
// Config
// ----------------------------------------------------------------------------

// Config-Struct using 64 Bit variables for config params (should only be used if necessary).
template <  __uint64 PARTITION_SIZE_= 32 * 1024 * 1024, 
            __uint64 OUTBUF_SIZE_   = 1 * 1024, 
            __uint64 INBUF_SIZE_    = 64 * 1024, 
            __uint64 TAIL_LENGTH_   = 1000,
            unsigned PREFIX_LENGTH_ = 32,
            typename TPrefixElem_   = unsigned>
struct DigestConfigLarge {
    typedef __uint64                                        Type;
    typedef TPrefixElem_                                    TPrefixElem;

    static const Type PARTITION_SIZE        = PARTITION_SIZE_;
    static const Type OUTBUF_SIZE           = OUTBUF_SIZE_;
    static const Type INBUF_SIZE            = INBUF_SIZE_;
    static const Type TAIL_LENGTH           = TAIL_LENGTH_;
    static const TPrefixElem PREFIX_LENGTH  = PREFIX_LENGTH_;
    
    // Configuration for suffix arrays
    template <typename TFile = ExternalConfigLarge<>::TFile, unsigned PAGESIZE = ExternalConfigLarge<>::PAGESIZE, unsigned FRAMES = ExternalConfigLarge<>::FRAMES>
    struct SuffixArrayConfig {
        typedef ExternalConfigLarge<TFile, PAGESIZE, FRAMES> Type;
    };
};

// Default config struct using 32 Bit variables for config params.
template <  unsigned PARTITION_SIZE_= 32 * 1024 * 1024,
            unsigned OUTBUF_SIZE_   = 1 * 1024,
            unsigned INBUF_SIZE_    = 64 * 1024,
            unsigned TAIL_LENGTH_   = 1000,
            unsigned PREFIX_LENGTH_ = 32,
            typename TPrefixElem_   = unsigned>
struct DigestConfig {
    typedef unsigned                                        Type;
    typedef TPrefixElem_                                    TPrefixElem;

    static const Type PARTITION_SIZE        = PARTITION_SIZE_;
    static const Type OUTBUF_SIZE           = OUTBUF_SIZE_;
    static const Type INBUF_SIZE            = INBUF_SIZE_;
    static const Type TAIL_LENGTH           = TAIL_LENGTH_;
    static const TPrefixElem PREFIX_LENGTH  = PREFIX_LENGTH_;

    // Configuration for suffix arrays
    template <typename TFile = ExternalConfig<>::TFile, unsigned PAGESIZE = ExternalConfig<>::PAGESIZE, unsigned FRAMES = ExternalConfig<>::FRAMES>
    struct SuffixArrayConfig {
        typedef ExternalConfig<TFile, PAGESIZE, FRAMES> Type;
    };
};

// Define index tag with default config
template <typename TConfig = DigestConfig<> >
struct IndexDigest {};

// ----------------------------------------------------------------------------
// Internal node of suffix tree
// ----------------------------------------------------------------------------

template<typename TBufSize_, typename TSize>
struct InternalNodeDigest_ {
public:
    typedef TBufSize_                                       TBufSize;

    // Invalid index value indicating that there is no child
    static const TBufSize NONE = ~(TBufSize)0;
    // static const TBufSize NONE = MaxValue<unsigned>::VALUE;

    // TODO(krugel) Using the meta function MaxValue yields a compile error in VS2015: 
    //     error C2131: expression did not evaluate to a constant
    //     failure was caused by non - constant arguments or reference to a non - constant symbol
    // To fix this issue (bug in VS?), change alphabet_math.h
    //     struct MaximumValueUnsigned_ { static const T VALUE = ~(T)0; };

    // Index of left child in the tree
    TBufSize leftChild;
    // Index of right child in the tree
    TBufSize rightChild;
    // Left index of leaf interval (suffix positions) represented by this node
    TBufSize leftmostLeaf;
    // String depth of this node in bits within the suffix tree (length of the processed prefix from the root to this node)
    TSize repLength;

    InternalNodeDigest_() {}

    InternalNodeDigest_(TSize repLength, TBufSize leftmostLeaf, TBufSize leftChild, TBufSize rightChild):
        leftChild(leftChild), rightChild(rightChild), leftmostLeaf(leftmostLeaf), repLength(repLength) { }
};

// ----------------------------------------------------------------------------
// Index
// ----------------------------------------------------------------------------

template <typename TObject, typename TConfig_>
class Index<TObject, IndexDigest<TConfig_> > {
public:
    typedef TConfig_                                        TConfig;
    typedef typename Fibre<Index, DigestText>::Type         TText;
    typedef typename Value<TText>::Type                     TTextValue;
    typedef typename TConfig::Type                          TBufSize;
    
    // Container type to store the bit values of the prefix,
    // must be unsigned and at least as big as TTextValue.
    typedef typename TConfig::TPrefixElem                   TPrefixElem;
    
    static const TBufSize PARTITION_SIZE    = TConfig::PARTITION_SIZE;
    static const TBufSize OUTBUF_SIZE       = TConfig::OUTBUF_SIZE;
    static const TBufSize INBUF_SIZE        = TConfig::INBUF_SIZE;
    static const TBufSize TAIL_LENGTH       = TConfig::TAIL_LENGTH;
    static const unsigned PREFIX_LENGTH     = TConfig::PREFIX_LENGTH;
    static const unsigned PREFIX_ARRAY_SIZE = (PREFIX_LENGTH - 1u) / BitsPerValue<TPrefixElem>::VALUE + 1u;
    
    typedef typename Cargo<Index>::Type                     TCargo;
    typedef Holder<TText>                                   THolder;
    typedef typename Size<TText>::Type                      TSize;
    typedef typename Position<TText>::Type                  TPosition;
    typedef InternalNodeDigest_<TBufSize, TSize>            TNode;

    typedef String<TNode, Array<OUTBUF_SIZE> >              TNodes;
    typedef String<TPosition, Array<OUTBUF_SIZE> >          TLeaves;

    // A suffix tree consists of
    // 1. the actual nodes of the tree,
    // 2. the leaves (suffix positions) = suffix array, and
    // 3. the offset of this suffix array within the (virtual) global suffix array.
    typedef Triple<TNodes, TLeaves, TSize>                  TSuffixTree;

    // Array of suffix trees.
    // Only one suffix tree will be loaded into main memory at a time.
    //typedef String<TSuffixTree, MMap<ExternalConfig<File<>, 1u, 2u> > > TSuffixTrees;  // Does not work for longer texts
    typedef String<TSuffixTree, External<ExternalConfig<File<>, 1u, 3u> > > TSuffixTrees;
    //typedef String<TSuffixTree> TSuffixTrees;
    
    // The prefix is stored as a list of containers.
    // The characters of the prefix are concatenated in binary form and 
    // the concatenation is then split up and stored in the container types.
    typedef String<TPrefixElem, Array<PREFIX_ARRAY_SIZE> >  TPrefix;

    // A DiGeST suffix consists of
    // 1. the starting position and
    // 2. a short prefix of the suffix to speed up computations later on.
    typedef Pair<TPosition, TPrefix>                        TDigestSuffix;

    // Configuration for the suffix arrays.
    typedef typename TConfig::template SuffixArrayConfig<File<>, INBUF_SIZE, 2u>::Type TSuffixArrayConfig;

    // A suffix array.
    //typedef String<TDigestSuffix, MMap<TSuffixArrayConfig> > TSuffixArray;  // Does not work for longer texts
    typedef String<TDigestSuffix, External<TSuffixArrayConfig> > TSuffixArray;
    //typedef String<TDigestSuffix>                         TSuffixArray;
    // List of all suffix arrays.
    typedef String<TSuffixArray>                            TSuffixArrays;

    // The divider array.
    typedef String<TDigestSuffix>                           TDividers;

    // A partition consists of
    // 1. the start position,
    // 2. the end position (including the tail), and
    // 3. the length (excluding the tail).
    typedef Triple<TPosition, TPosition, TSize>             TPartition;
    typedef String<TPartition>                              TPartitions;

    THolder         text;
    TCargo          cargo;        // user-defined cargo
    TPartitions     partitions;
    TSuffixArrays   suffixArrays;
    TDividers       dividers;
    TSuffixTrees    suffixTrees;
    bool clearRedundantFibres; // Some fibres can be cleared as soon as they are no longer needed
    // (E.g. the suffix arrays can be removed as soon as the suffix trees are created.)

    Index() {}
    
    Index(Index &other):
        text(other.text),
        cargo(other.cargo),
        partitions(other.partitions),
        suffixArrays(other.suffixArrays),
        dividers(other.dividers),
        suffixTrees(other.suffixTrees),
        clearRedundantFibres(other.clearRedundantFibres) {}
    
    Index(Index const &other):
        text(other.text),
        cargo(other.cargo),
        partitions(other.partitions),
        suffixArrays(other.suffixArrays),
        dividers(other.dividers),
        suffixTrees(other.suffixTrees),
        clearRedundantFibres(other.clearRedundantFibres) {}
    
    template <typename TText_>
    Index(TText_ &_text):
        text(_text),
        clearRedundantFibres(true) { }
    
    //template <typename TText_>
    //Index(TText_ const &_text):
    //    text(_text),
    //    clearRedundantFibres(true) { }

    template <typename TText_>
    Index(TText_ &_text, bool _clearRedundantFibres):
        text(_text),
        clearRedundantFibres(_clearRedundantFibres) { }
    
};

// ----------------------------------------------------------------------------
// Lexical
// ----------------------------------------------------------------------------

template <typename WithLcp, unsigned PREFIX_LENGTH, typename TSize, typename TText1, typename TText2>
struct LexicalSuffixDigest {};

// Compute and store the result of the lexical comparison of two suffixes.
// Also computes the length of the longest common prefix (lcp) in bit and the first bit after the lcp.
template <typename WithLcp, unsigned PREFIX_LENGTH, typename TSize, typename TText1, typename TText2>
struct Lexical<LexicalSuffixDigest<WithLcp, PREFIX_LENGTH, TSize, TText1, TText2> > {
public:
    typedef typename Size<Lexical>::Type                    TLexicalSize;
    typedef typename Position<TText1>::Type                 TPosition1;
    typedef typename Position<TText2>::Type                 TPosition2;

    enum TResult {
        EQUAL = 1u,
        LESS = 2u,
        GREATER = 4u,
        LEFT_IS_PREFIX = 8u,
        RIGHT_IS_PREFIX = 16u
    };

    TText1 const & text1;
    TText2 const & text2;
    TResult data_compare;
    TLexicalSize data_lcp;
    bool bitAfterLcp2;

    Lexical(TText1 const & text1_, TText2 const & text2_):
        text1(text1_), text2(text2_) { /*std::cerr << "[" << text2 << "]" << std::endl;*/ }

    Lexical(Lexical const & other):
        text1(other.text1),
        text2(other.text2),
        data_compare(other.data_compare),
        data_lcp(other.data_lcp),
        bitAfterLcp2(other.bitAfterLcp2) { }

    template <typename TDigestSuffix1, typename TDigestSuffix2>
    Lexical(TText1 const & text1, TDigestSuffix1 const & left, TText2 const & text2, TDigestSuffix2 const & right):
        text1(text1), text2(text2)
    {
        compare(*this, left, right);
    }

    ~Lexical() {}
};

// Comparator to determine the greater of two suffixes, given by iterators to suffix arrays.
// Will be used in a priority queue.
template <typename TText, unsigned PREFIX_LENGTH>
struct GreaterSuffixDigest {
    typedef typename Size<TText>::Type                      TSize;
    
    TText const & text;

    GreaterSuffixDigest(TText & _text):
        text(_text) { }
    
    template <typename TSuffixArrayIterator>
    bool operator() (const TSuffixArrayIterator * it1, const TSuffixArrayIterator * it2) const {
        Lexical<LexicalSuffixDigest<False, PREFIX_LENGTH, TSize, TText, TText> > lex(text, **it1, text, **it2);
        return isGreater(lex) || hasPrefix(lex);
    }
};


// ============================================================================
// Metafunctions
// ============================================================================

template <typename TText, typename TSpec>
struct SAValue<Index<TText, IndexDigest<TSpec> > > {
    typedef typename Size<TText>::Type                      Type;
};

template <typename WithLcp, unsigned PREFIX_LENGTH, typename TSize, typename TText1, typename TText2>
struct Size<Lexical<LexicalSuffixDigest<WithLcp, PREFIX_LENGTH, TSize, TText1, TText2> > > {
    typedef TSize                                           Type;
};

template <typename WithLcp, unsigned PREFIX_LENGTH, typename TSize, typename TText1, typename TText2>
struct Size<Lexical<LexicalSuffixDigest<WithLcp, PREFIX_LENGTH, TSize, TText1, TText2> > const> {
    typedef TSize                                           Type;
};

template<typename TText, typename TSpec>
struct Fibre<Index<TText, IndexDigest<TSpec> >, DigestPartitions> {
    typedef Index<TText, IndexDigest<TSpec> >               TIndex;
    typedef typename TIndex::TPartitions                    Type;
};

template<typename TText, typename TSpec>
struct Fibre<Index<TText, IndexDigest<TSpec> >, DigestSuffixArrays> {
    typedef Index<TText, IndexDigest<TSpec> >               TIndex;
    typedef typename TIndex::TSuffixArrays                  Type;
};

template<typename TText, typename TSpec>
struct Fibre<Index<TText, IndexDigest<TSpec> >, DigestDividers> {
    typedef Index<TText, IndexDigest<TSpec> >               TIndex;
    typedef typename TIndex::TDividers                      Type;
};

template<typename TText, typename TSpec>
struct Fibre<Index<TText, IndexDigest<TSpec> >, DigestSuffixTrees> {
    typedef Index<TText, IndexDigest<TSpec> >               TIndex;
    typedef typename TIndex::TSuffixTrees                   Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Create partitions
// ----------------------------------------------------------------------------

// Determine the end position of the partition after appending a tail to ensure
// no suffix is a prefix of another suffix.
template <unsigned TAIL_LENGTH, typename TText, typename TPosition>
inline TPosition
_getTailEndPosition(TText const & text, TPosition const & posEnd, TPosition const & posMax) {
    ignoreUnusedVariableWarning(text);
    TPosition endTail = posAdd(posEnd, TAIL_LENGTH);
    return _min(endTail, posMax);
}

// Create partitions for one text (default version).
// Partitions a text into smaller pieces and appends a tail to each piece
// in order to ensure no suffix is a prefix of another.
template <typename TText, typename TSpec>
inline bool
_createPartitions(Index<TText, IndexDigest<TSpec> > & index) {
    typedef Index<TText, IndexDigest<TSpec> >               TIndex;
    typedef typename TIndex::TSize                          TSize;
    typedef typename TIndex::THolder                        THolder;
    typedef typename Reference<THolder>::Type               TTextReference;
    typedef typename TIndex::TPartitions                    TPartitions;
    typedef typename Iterator<TPartitions>::Type            TPartitionsIterator;
    typedef typename Reference<TPartitionsIterator>::Type   TPartitionReference;
    typedef typename TIndex::TPosition                      TPosition;

    TTextReference text = indexText(index);
    TSize textLength = length(text);
    TPosition pos = beginPosition(text);
    TPosition endPos = endPosition(text);

    TSize numPartitions = textLength / TIndex::PARTITION_SIZE;

    if (textLength % TIndex::PARTITION_SIZE != 0u) {
        ++numPartitions;
    }

    resize(index.partitions, numPartitions, Exact());

    for (TPartitionsIterator it = begin(index.partitions); !atEnd(it, index.partitions); ++it) {
        TPosition end = posAdd(pos, TIndex::PARTITION_SIZE);
        if (end > endPos) end = endPos;
        TPosition endTail = _getTailEndPosition<TIndex::TAIL_LENGTH>(text, end, endPos);
        TPartitionReference partition = value(it);
        partition.i1 = pos;
        partition.i2 = endTail;
        partition.i3 = posSub(end, pos);
        pos = end;
    }
    return true;
}

// Create partitions for a set of strings.
// Partitions the strings from a string set into smaller pieces and appends to each piece a tail
// in order to ensure no suffix is a prefix of another.
template <typename TText, typename TSetSpec, typename TSpec>
inline bool
_createPartitions(Index<StringSet<TText, TSetSpec>, IndexDigest<TSpec> > & index) {
    typedef StringSet<TText, TSetSpec>                      TStringSet;
    typedef Index<TStringSet, IndexDigest<TSpec> >          TIndex;
    typedef typename TIndex::TSize                          TSize;
    typedef typename Position<TIndex>::Type                 TPosition;
    typedef typename TIndex::THolder                        THolder;
    typedef typename Reference<THolder>::Type               TTextReference;
    typedef typename TIndex::TPartitions                    TPartitions;
    typedef typename Iterator<TPartitions>::Type            TPartitionsIterator;
    typedef typename Reference<TPartitionsIterator>::Type   TPartitionReference;
    typedef typename Iterator<TStringSet>::Type             TStringSetIterator;
    typedef typename Reference<TStringSetIterator>::Type    TStringReference;

    TTextReference stringSet = indexText(index);
    TSize totalNumPartitions = 0u;

    for (TStringSetIterator itSet = begin(stringSet); !atEnd(itSet, stringSet); ++itSet) {
        TStringReference text = value(itSet);
        TSize textLength = length(text);
        TSize numPartitions = textLength / TIndex::PARTITION_SIZE;

        if (textLength % TIndex::PARTITION_SIZE != 0u) {
            ++numPartitions;
        }
        totalNumPartitions += numPartitions;
    }

    resize(index.partitions, totalNumPartitions, Exact());

    TPartitionsIterator itPart = begin(index.partitions);

    for (TStringSetIterator itSet = begin(stringSet); !atEnd(itSet, stringSet); ++itSet) {
        TStringReference text = value(itSet);
        TSize textLength = length(text);
        TPosition pos = beginPosition(text);
        TPosition endPos = endPosition(text);

        while (posLess(pos, endPos)) {
            TPosition end = posAdd(pos, TIndex::PARTITION_SIZE);
            if (posLess(endPos, end)) end = endPos;
            TPosition endTail = _getTailEndPosition<TIndex::TAIL_LENGTH>(text, end, endPos);
            TPartitionReference partition = value(itPart);
            partition.i1 = pos;
            partition.i2 = endTail;
            partition.i3 = posSub(end, pos);
            pos = end;
            ++itPart;
        }
    }
    return true;
}

// ----------------------------------------------------------------------------
// Create suffix arrays
// ----------------------------------------------------------------------------

// Given a text and position within this text this function computes the prefix
// of the corresponding suffix starting at this position.
template <unsigned PREFIX_LENGTH, typename TText, typename TPosition, typename TPrefix>
inline void
_getPrefix(TPrefix & prefixArray, TText const & text, TPosition pos) {
    typedef typename Value<TPrefix>::Type                   TPrefixElem;
    typedef typename Iterator<TPrefix>::Type                TPrefixIter;
    typedef typename Value<TText>::Type                     TValue;
    typedef typename Size<TText>::Type                      TSize;
    typedef typename MakeSigned_<TSize>::Type               TSizeSigned;
    typedef typename Iterator<const TText>::Type            TTextIter;

    typedef typename Value<TText>::Type                     TValue;
    typedef typename MakeUnsigned_<TValue>::Type            TMask;      // Yields e.g. unsigned char or SimpleType<unsigned char, Dna_>.
    typedef typename Value<TMask>::Type                     TMaskValue; // Yields in both above cases unsigned char.

    static const TSizeSigned BITS_PER_VALUE       = BitsPerValue<TValue>::VALUE;
    static const TSizeSigned BITS_PER_PREFIX_ELEM = BitsPerValue<TPrefixElem>::VALUE;
    static const TSize       PREFIX_ARRAY_SIZE    = LENGTH<TPrefix>::VALUE;
    
    //TSize textLength = length(text);
    TSizeSigned paOffset = BITS_PER_PREFIX_ELEM - BITS_PER_VALUE;
    resize(prefixArray, PREFIX_ARRAY_SIZE, 0u, Exact());
    //std::cerr << "PREFIX_ARRAY_SIZE = " << (unsigned int) PREFIX_ARRAY_SIZE << std::endl;
    //std::cerr << "capacity          = " << (unsigned int) capacity(prefixArray) << std::endl;
    //std::cerr << "length            = " << (unsigned int) length(prefixArray) << std::endl;
    //std::cerr << "PREFIX_LENGTH     = " << (unsigned int) PREFIX_LENGTH << std::endl;
    //std::cerr << "BITS_PER_VALUE    = " << (unsigned int) BITS_PER_VALUE << std::endl;
    //std::cerr << "BITS_PER_PREFIX_E = " << (unsigned int) BITS_PER_PREFIX_ELEM << std::endl;
    //std::cerr << std::endl;

    // Points to the current position in the text (initially the beginning of the suffix starting at pos)
    
    TPrefixIter prefixIter = begin(prefixArray);
    TPrefixIter prefixEndIter = end(prefixArray);
    TTextIter textEndIter = end(text);
    for(TTextIter textIter = iter(text, pos); textIter != textEndIter; ++textIter) {
        //std::cerr << "i                 = " << (unsigned int) i << std::endl;
        //std::cerr << "paIndex           = " << (unsigned int) paIndex << std::endl;
        //std::cerr << "paOffset          = " << (int) paOffset << std::endl;

        //if (posAtEnd(pos, text)) return;
        //TMaskValue val2 = (TMaskValue) text[pos];
        //TMaskValue val = convert<TMaskValue>(*textIter);
        TMaskValue val = convert<TMaskValue>(*textIter);
        if (paOffset <= 0) {
            //SEQAN_ASSERT_LT(paIndex, PREFIX_ARRAY_SIZE);
            value(prefixIter) |= (val >> -paOffset);
            //++paIndex;
            ++prefixIter;
            if (prefixIter == prefixEndIter) break;
            paOffset += BITS_PER_PREFIX_ELEM;
            if (paOffset < BITS_PER_PREFIX_ELEM) {
                //SEQAN_ASSERT_LT(paIndex, PREFIX_ARRAY_SIZE);
                value(prefixIter) |= (val << paOffset);
            }
        } else {
            //SEQAN_ASSERT_LT(paIndex, PREFIX_ARRAY_SIZE);
            value(prefixIter) |= (val << paOffset);
        }
        paOffset -= BITS_PER_VALUE;
        //pos = posNext(pos);
    }
}

// Computes the suffix arrays for the partitions.
template <typename TText_, typename TSpec, typename TAlgorithm>
inline bool
_createSuffixArrays(Index<TText_, IndexDigest<TSpec> > & index, TAlgorithm const &) {
    typedef Index<TText_, IndexDigest<TSpec> >              TIndex;
    typedef typename TIndex::TText                          TText;
    typedef typename Infix<TText>::Type                     TInfix;
    typedef typename TIndex::THolder                        THolder;
    typedef typename Reference<THolder>::Type               TTextReference;
    typedef typename TIndex::TPrefix                        TPrefix;
    typedef typename TIndex::TSize                          TSize;
    typedef typename TIndex::TPosition                      TPosition;
    typedef typename TIndex::TSuffixArray                   TSuffixArray;
    typedef typename Iterator<TSuffixArray>::Type           TSuffixArrayIterator;
    typedef typename Reference<TSuffixArrayIterator>::Type  TSuffixArrayElementReference;
    //typedef typename Value<TSuffixArrayIterator>::Type      TSuffixArrayElement;
    typedef typename TIndex::TPartitions                    TPartitions;
    typedef typename Iterator<TPartitions>::Type            TPartitionsIterator;
    typedef typename Reference<TPartitionsIterator>::Type   TPartitionReference;
    typedef typename TIndex::TSuffixArrays                  TSuffixArrays;
    typedef typename Iterator<TSuffixArrays>::Type          TSuffixArraysIterator;
    typedef typename Reference<TSuffixArraysIterator>::Type TSuffixArrayReference;
    typedef String<TSize>                                   TSuffixArrayRaw;
    typedef typename Iterator<TSuffixArrayRaw>::Type        TSuffixArrayRawIterator;

    static const TSize PREFIX_ARRAY_SIZE = TIndex::PREFIX_ARRAY_SIZE;
    ignoreUnusedVariableWarning(PREFIX_ARRAY_SIZE);

    indexRequire(index, DigestPartitions());

    //std::cerr << "_createSuffixArrays" << std::endl;
    
    TSize numPartitions = length(index.partitions);
    resize(index.suffixArrays, numPartitions, Exact());
    TTextReference text = indexText(index);
    
    TSuffixArraysIterator sasIter = begin(index.suffixArrays);
    for (TPartitionsIterator pIter = begin(index.partitions); !atEnd(pIter, index.partitions); ++pIter, ++sasIter) {
        TPartitionReference partitionSpecs = value(pIter);
        TPosition beginPos = partitionSpecs.i1;
        TPosition endPos = partitionSpecs.i2;
        TPosition lengthWithoutTail = partitionSpecs.i3;
        TInfix partition = infix(text, beginPos, endPos);
        TSize lengthWithTail = length(partition);
        
        // We create the suffix array for the partition including the tail,
        // so that the ordering is valid also globally for the whole string.
        TSuffixArrayRaw saRaw;
        //std::cerr << position(pIter)+1 << "/" << numPartitions << " ";
        //std::cerr << "r";
        resize(saRaw, lengthWithTail, Exact());
        
        //std::cerr << "c";
        createSuffixArray(saRaw, partition, TAlgorithm());
        
        TSuffixArrayReference sa = value(sasIter);

        //std::cerr << "f";

        resize(sa, lengthWithoutTail, Exact());
        //mmapAdvise(sa, MAP_SEQUENTIAL);
        //std::cerr << "i";
        TSuffixArrayIterator saIter = begin(sa);
        for (TSuffixArrayRawIterator saRawIter = begin(saRaw); !atEnd(saRawIter, saRaw); ++saRawIter) {
            // Finally, we only store the values that refer to the partition without tail.
            if (*saRawIter < lengthWithoutTail) {
                TSuffixArrayElementReference saElem = *saIter;
                // We store the global starting position of the suffix,
                // together with a short binary prefix of the suffix.
                saElem.i1 = posAdd(beginPos, *saRawIter);
                saElem.i2 = TPrefix(0u);
                _getPrefix<TIndex::PREFIX_LENGTH>(saElem.i2, text, saElem.i1);
                ++saIter;
            }
        }
        //std::cerr << "d ";
        
/*-
        reserve(sa, lengthWithoutTail, Exact());

        for (TSuffixArrayRawIterator saRawIter = begin(saRaw); !atEnd(saRawIter, saRaw); ++saRawIter) {
            // Finally, we only store the values that refer to the partition without tail.
            if (*saRawIter < lengthWithoutTail) {
                TSuffixArrayElement saElem;
                // We store the global starting position of the suffix,
                // together with a short binary prefix of the suffix.
                saElem.i1 = posAdd(beginPos, *saRawIter);
                saElem.i2 = TPrefix(0u);
                _getPrefix<TIndex::PREFIX_LENGTH>(saElem.i2, text, saElem.i1);
                appendValue(sa, saElem);
                SEQAN_ASSERT_LEQ(length(saElem.i2), PREFIX_ARRAY_SIZE);
            }
        }
//*/
    }    

    SEQAN_ASSERT_LEQ(length(front(front(index.suffixArrays)).i2), PREFIX_ARRAY_SIZE);
    SEQAN_ASSERT_LEQ(length(back(back(index.suffixArrays)).i2), PREFIX_ARRAY_SIZE);
    
    if (index.clearRedundantFibres) {
        clear(index.partitions);
    }
    //std::cerr << "return";
    return true;
}

// ----------------------------------------------------------------------------
// Debug output
// ----------------------------------------------------------------------------

template<typename T>
std::string
convertToStr(T value) {
    std::stringstream ss;
    ss << value;
    return ss.str();
}

template <typename TString, typename TStringDest>
inline void
_getBitRepresentation(TStringDest & dest, TString const & string) {
    typedef typename Iterator<const TString>::Type          TIterator;
    typedef typename Value<TString>::Type                   TValue;
    typedef typename MakeUnsigned_<TValue>::Type            TMask;      // Yields e.g. unsigned char or SimpleType<unsigned char, Dna_>.
    typedef typename Value<TMask>::Type                     TMaskValue; // Yields in both above cases unsigned char.

    static const unsigned BITS_PER_VALUE = BitsPerValue<TValue>::VALUE;
    static const TMaskValue MASK_HIGHEST_BIT = 1u << (BITS_PER_VALUE - 1u);
    
    for (TIterator it = begin(string); !atEnd(it, string); ++it) {
        TMaskValue val = convert<TMaskValue>(*it);
        TMaskValue mask = MASK_HIGHEST_BIT;
        while (mask != 0u) {
            appendValue(dest, ((mask & val) != 0u) ? '1' : '0');
            mask >>= 1u;
        }
    }
}

template <typename TBufSize, typename TSize>
inline void
_printPreorder(TBufSize root, String<InternalNodeDigest_<TBufSize, TSize> > const & nodes) {
    typedef InternalNodeDigest_<TBufSize, TSize>            TNode;
    
    std::cerr << std::endl << "NODE: " <<root<<std::endl;
    std::cerr << "RepLength: " << nodes[root].repLength << std::endl; 
    std::cerr << "LeftmostLeaf: " << nodes[root].leftmostLeaf << std::endl; 
    std::cerr << "LeftChild: " << nodes[root].leftChild << std::endl; 
    std::cerr << "RightChild: " << nodes[root].rightChild << std::endl; 
    
    if (nodes[root].leftChild != TNode::NONE) {
        _printPreorder(nodes[root].leftChild, nodes);
    }
    if (nodes[root].rightChild != TNode::NONE) {
        _printPreorder(nodes[root].rightChild, nodes);
    }
}

template <typename TSuffixTree>
inline void
_printSuffixTree(TSuffixTree const & suffixTree) {
    //_printPreorder(0u, suffixTree.i1);
    std::cerr << std::endl << "Leaves:" << std::endl;
    for (unsigned i = 0u; i < length(suffixTree.i2); ++i) {
        std::cerr << suffixTree.i2[i] << " " <<std::endl;
    }
}

template <typename TSuffixTrees>
inline void
_printSuffixTrees(TSuffixTrees const & suffixTrees) {
    for (unsigned i = 0u; i < length(suffixTrees); ++i) {
        std::cerr << "-----------------------------------------------" << std::endl;
        std::cerr << "------------------- NEW TREE ------------------" << std::endl;
        std::cerr << "-----------------------------------------------" << std::endl;
        _printSuffixTree(suffixTrees[i]);
    }
}

// ----------------------------------------------------------------------------
// Lexical comparison
// ----------------------------------------------------------------------------

template <typename TPrefix>
inline bool
_getFirstBit(TPrefix const & pref) {
    typedef typename Value<TPrefix>::Type                   TPrefixElem;

    static const unsigned BITS_PER_PREFIX_ELEM = BitsPerValue<TPrefixElem>::VALUE;
    static const TPrefixElem MASK_HIGHEST_BIT_PREFIX = 1u << (BITS_PER_PREFIX_ELEM - 1u); // 10000..0

    return pref[0u] & MASK_HIGHEST_BIT_PREFIX;
}

template <typename WithLcp, unsigned PREFIX_LENGTH, typename TSize, typename TText1, typename TText2>
inline bool
_getBitAfterLcp1(Lexical<LexicalSuffixDigest<WithLcp, PREFIX_LENGTH, TSize, TText1, TText2> > const & lex) {
    return !lex.bitAfterLcp2;
}

template <typename WithLcp, unsigned PREFIX_LENGTH, typename TSize, typename TText1, typename TText2>
inline bool
_getBitAfterLcp2(Lexical<LexicalSuffixDigest<WithLcp, PREFIX_LENGTH, TSize, TText1, TText2> > const & lex) {
    return lex.bitAfterLcp2;
}

// lcp counted in characters
template <typename TText1, typename TText2, typename TSize>
inline TSize
_lcpLength(TText1 const & text1, TSize pos1, TText2 const & text2, TSize pos2) {
    typedef typename Iterator<const TText1>::Type           TTextIter1;
    typedef typename Iterator<const TText2>::Type           TTextIter2;
    
    //static unsigned int debugNumCalls = 0u;
    //static unsigned int debugLcpLengthSum = 0u;
    //++debugNumCalls;

    TTextIter1 iter1 = iter(text1, pos1);
    TTextIter2 iter2 = iter(text2, pos2);
    
    TSize maxLength = std::min(length(text1) - pos1, length(text2) - pos2);
    TSize lcpRest2;
    for (lcpRest2 = 0u; lcpRest2 < maxLength; ++lcpRest2, ++iter1, ++iter2) {
        if (*iter1 != *iter2) {
            break;
        }
    }

    //debugLcpLengthSum += lcpRest2;
    //if (debugNumCalls % 100000u == 0u) {
    //    std::cerr << debugLcpLengthSum << " / " << debugNumCalls << " ~ " << double(debugLcpLengthSum) / double(debugNumCalls) << ", ";
    //}

    return lcpRest2;
    //return lcpLength(suffix(text, pos1), suffix(text, pos2));
}

// lcp counted in bits
template <typename WithLcp, unsigned PREFIX_LENGTH, typename TSize, typename TText1, typename TPosition1, typename TText2, typename TPosition2, typename TPrefix>
inline void
compare_(
        Lexical<LexicalSuffixDigest<WithLcp, PREFIX_LENGTH, TSize, TText1, TText2> > & lexical,
        Pair<TPosition1, TPrefix> const & suffix1,
        Pair<TPosition2, TPrefix> const & suffix2,
        bool debug = false) {

    ignoreUnusedVariableWarning(debug);
    //static unsigned int debugNumCalls = 0u;
    //static unsigned int debugNumCallsExpensive = 0u;
    //++debugNumCalls;
    //if (debugNumCalls % 1000000u == 0u) {
    //    std::cerr << debugNumCallsExpensive << " / " << debugNumCalls << " = " << double(debugNumCallsExpensive) / double(debugNumCalls) << ", ";
    //}

    typedef Lexical<LexicalSuffixDigest<WithLcp, PREFIX_LENGTH, TSize, TText1, TText2> > TLexical;
    typedef typename Value<TPrefix>::Type                   TPrefixElem;
    //typedef typename Size<TPrefix>::Type                    TPrefixSize;
    typedef typename Iterator<const TPrefix>::Type          TPrefixIter;
    typedef typename Value<TText1>::Type                    TValue;
    typedef typename MakeUnsigned_<TValue>::Type            TMask;      // Yields e.g. unsigned char or SimpleType<unsigned char, Dna_>.
    typedef typename Value<TMask>::Type                     TMaskValue; // Yields in both above cases unsigned char.
    typedef typename Comparator<TPrefix>::Type              TPrefixComparator;

    static const TSize       PREFIX_ARRAY_SIZE = LENGTH<TPrefix>::VALUE;
    static const unsigned BITS_PER_PREFIX_ELEM = BitsPerValue<TPrefixElem>::VALUE;
    static const unsigned BITS_PER_VALUE       = BitsPerValue<TValue>::VALUE;
    static const unsigned PREFIX_LENGTH_STRING = PREFIX_LENGTH / BITS_PER_VALUE;
    static const TPrefixElem MASK_HIGHEST_BIT_PREFIX = 1u << (BITS_PER_PREFIX_ELEM - 1u); // 10000..0
    static const TMaskValue  MASK_HIGHEST_BIT        = 1u << (BITS_PER_VALUE - 1u);       // 10000..0
    static const bool withLcp = WithLcp::VALUE;

    ignoreUnusedVariableWarning(PREFIX_ARRAY_SIZE);

    // Lengths in bits
    TSize length1 = posSub(length(lexical.text1), suffix1.i1) * BITS_PER_VALUE;
    TSize length2 = posSub(length(lexical.text2), suffix2.i1) * BITS_PER_VALUE;
    TSize prefixLength = _min(_min(PREFIX_LENGTH, length1), length2);
    
    // Shortcuts
    TPrefix const & prefix1 = suffix1.i2;
    TPrefix const & prefix2 = suffix2.i2;
    // if (debug) std::cerr << __LINE__ << ";";

    // Speedup if we don't need the exact lcp, we can simply compare the prefixes first
    if (! withLcp && length1 > PREFIX_LENGTH && length2 > PREFIX_LENGTH) {
        TPrefixComparator prefixComp(prefix1, prefix2);
        lexical.data_lcp = 0;
        if (isLess(prefixComp)) {
            lexical.data_compare = TLexical::LESS;
            return;
        } else if (isGreater(prefixComp)) {
            lexical.data_compare = TLexical::GREATER;
            return;
        }
    }

    //++debugNumCallsExpensive;

    // Compare prefixes.
    TPrefixIter prefixIter1 = begin(prefix1), prefixIter2 = begin(prefix2);
    
    // if (debug) std::cerr << __LINE__ << ";";
    TPrefixElem mask = 0u;
    for (lexical.data_lcp = 0u; lexical.data_lcp < prefixLength;) {
        // If the whole prefix element is identical skip bitwise comparison
        //SEQAN_ASSERT_LT((unsigned int) i, PREFIX_ARRAY_SIZE);
        if (lexical.data_lcp + BITS_PER_PREFIX_ELEM <= prefixLength && value(prefixIter1) == value(prefixIter2)) {
            lexical.data_lcp += BITS_PER_PREFIX_ELEM;
            ++prefixIter1;
            ++prefixIter2;
            if (lexical.data_lcp == prefixLength) {
                mask = MASK_HIGHEST_BIT_PREFIX;
            }
        } else {
            // Try to find the differing bit.
            TPrefixElem diff = value(prefixIter1) ^ value(prefixIter2);
            // std::cerr << "suffix1: " << value(prefixIter1) << " | suffix2: " << value(prefixIter2) << std::endl;
            // std::cerr << "Diff: " << diff << std::endl;
            for (mask = MASK_HIGHEST_BIT_PREFIX; lexical.data_lcp < prefixLength; mask >>= 1u, ++lexical.data_lcp) {
                // Bit differ at this position
                if ((diff & mask) != 0u) {
                    // std::cerr << "CALC LCP: " << lcp << std::endl;
                    lexical.bitAfterLcp2 = (mask & value(prefixIter2)) != 0u;
                    lexical.data_compare = lexical.bitAfterLcp2 ? TLexical::LESS : TLexical::GREATER;
                    return;
                }
            }
            // Can only be reached if both prefixes are identical.
        }
    }
    // if (debug) std::cerr << __LINE__ << ";";

    // Either suffix1 or suffix2 is shorter than PREFIX_LENGTH_STRING.
    if (prefixLength < PREFIX_LENGTH) {
        if (length2 == length1) {
            lexical.data_compare = TLexical::EQUAL;
            return;
        } else if (length2 > length1){
            // If not the whole prefix of suffix2 was compared (because suffix1 was shorter than
            // the maximum prefix length), we can extract the bit after the lcp (which is length1)
            // directly from the prefix of suffix2.
            // Mask has to be non-zero which is garantueed.
            // (otherwise the very first if of the previous for loop would have been true.)
            lexical.bitAfterLcp2 = (mask & value(prefixIter2)) != 0u;
            lexical.data_compare = TLexical::LEFT_IS_PREFIX;
            return;
        } else {
            lexical.bitAfterLcp2 = (mask & value(prefixIter1)) == 0u;
            lexical.data_compare = TLexical::RIGHT_IS_PREFIX;
            return;
        }
    }
    
    // Prefixes were identical and suffix1 and suffix2 are at least of size PREFIX_LENGTH_STRING.
    // Now we have to compare using the underlying text.

    // Start positions of the suffices.
    TPosition1 pos1 = suffix1.i1;
    TPosition2 pos2 = suffix2.i1;

    // TODO(krugel) Use iterators
    //typedef typename Iterator<const TText1, Rooted>::Type   TTextIter1;
    //typedef typename Iterator<const TText2>::Type           TTextIter2;
    
    //TTextIter1 iter1(iter(lexical.text1, pos1));
    //TTextIter2 iter2(iter(lexical.text2, pos2));

    // if (debug) std::cerr << __LINE__ << ";";
    // Skip prefixes which were already compared.
    bool pos1Check = posAddAndCheck(pos1, PREFIX_LENGTH_STRING, lexical.text1);
    bool pos2Check = posAddAndCheck(pos2, PREFIX_LENGTH_STRING, lexical.text2);
    lexical.data_lcp = PREFIX_LENGTH_STRING * BITS_PER_VALUE;

    if (pos1Check && pos2Check) {
        // The lcp in number of matching characters
        TSize lcpRest = _lcpLength(lexical.text1, pos1, lexical.text2, pos2);
        //TSize lcpRest = lcpLength(suffix(lexical.text1, pos1), suffix(lexical.text2, pos2));

        // The lcp in bits
        lexical.data_lcp += lcpRest * BITS_PER_VALUE;

        // Skip lcp
        pos1Check = posAddAndCheck(pos1, lcpRest, lexical.text1);
        pos2Check = posAddAndCheck(pos2, lcpRest, lexical.text2);
    }

    if (!pos1Check && !pos2Check) {
    // if (debug) std::cerr << __LINE__ << ";";
        lexical.data_compare = TLexical::EQUAL;
        return;
    } else if (!pos1Check) {
    // if (debug) std::cerr << __LINE__ << ";";
        // Now suffix2 is longer than suffix1 which has length PREFIX_LENGTH_STRING.
        // --> Extract bit after lcp from first character after prefix.
        lexical.bitAfterLcp2 = ((TMaskValue) value(lexical.text2, pos2) & MASK_HIGHEST_BIT) != 0u;
        lexical.data_compare = TLexical::LEFT_IS_PREFIX;
        return;
    } else if (!pos2Check) {
    // if (debug) std::cerr << __LINE__ << ";";
        lexical.bitAfterLcp2 = ((TMaskValue) value(lexical.text1, pos1) & MASK_HIGHEST_BIT) == 0u;
        lexical.data_compare = TLexical::RIGHT_IS_PREFIX;
        return;
    }

    // Determine bit in the first character in which suffix1 and suffix2 differ.
    TValue diffSuffix2 = value(lexical.text2, pos2);
    TMaskValue diff = (TMaskValue) value(lexical.text1, pos1) ^ (TMaskValue) diffSuffix2;

    TMaskValue mask2 = MASK_HIGHEST_BIT;
    while ((diff & mask2) == 0u) {
        ++lexical.data_lcp;
        mask2 >>= 1u;
    }
    lexical.bitAfterLcp2 = ((TMaskValue) diffSuffix2 & mask2) != 0u;
    lexical.data_compare = lexical.bitAfterLcp2 ? TLexical::LESS : TLexical::GREATER;
}

// ----------------------------------------------------------------------------
// Merge suffix arrays
// ----------------------------------------------------------------------------

// Insert one suffix into the suffix tree.
// Note: The suffix type has no partition id (this suffix type is only used directly in the merge function and the priority queue).
// Returns true, if further inserts may exceed the capacity limit.
// OUTBUF_SIZE = maximal number of leaves and internal nodes of this suffix tree
template <unsigned OUTBUF_SIZE, typename TNodes, typename TLeaves, typename TSize, typename TDigestSuffix>
inline bool
_insertSuffix(
        Triple<TNodes, TLeaves, TSize> & tree, // suffix tree 
        TDigestSuffix const & suffix,   // suffix (position, prefix)
        TSize suffixLength,             // length of the suffix
        TSize lcp,                      // lcp of the suffix with the last inserted element
        bool bitAfterLcp,               // first bit after lcp with last inserted element
        TSize & numLeaves,              // number of leaves
        TSize & numNodes,               // number of internal nodes
        bool debug = false) {           // output debug information?
    ignoreUnusedVariableWarning(debug);

    typedef typename Value<TNodes>::Type                TNode;
    typedef typename Reference<TNodes>::Type            TNodeReference;
    typedef typename TNode::TBufSize                    TBufSize;
    //typedef typename Iterator<TNodes, Rooted>::Type     TNodesIter;

    // Shortcuts
    TNodes & nodes = tree.i1;
    TLeaves & leaves = tree.i2;

    //if (debug) {
    //    std::cerr << "New leaf: " << suffix.i1 << std::endl;
    //    std::cerr << "suffixLength: " << suffixLength << std::endl;
    //    std::cerr << "LCP: " << lcp << std::endl;
    //    std::cerr << "First bit: " << (bitAfterLcp ? 1u : 0u) << std::endl;
    //}

    // TODO(krugel) Store boundary path: String<nodeId>, resize, appendValue, binary search (using nodes[nodeId].repLength)
    
    // Traverse tree on the right border up to depth lcp
    TBufSize oldParentId = 0u, oldChildId = 0u;
    
    // Follow the boundary path up to depth lcp.
    while ((value(nodes, oldChildId).repLength) < lcp) {
        oldParentId = oldChildId;
        TNodeReference node = value(nodes, oldChildId);
        if (node.rightChild != TNode::NONE) {
            oldChildId = node.rightChild;
        } else if (node.leftChild != TNode::NONE) {
            oldChildId = node.leftChild;
        } else {
            SEQAN_ASSERT_MSG(false, "Unexpectedly discovered a leaf node.");
        }
    }

    // Node at the depth lcp or child of the edge at depth lcp
    TNodeReference oldChild = nodes[oldChildId];
    // Parent node of the previous node
    TNodeReference oldParent = nodes[oldParentId];

/*
    // Traverse tree on the right border up to depth lcp
    TNodesIter parentIter = begin(nodes);
    TNodesIter childIter = begin(nodes);
    
    // Follow the boundary path up to depth lcp.
    while ((value(childIter).repLength) < lcp) {
        setPosition(parentIter, position(childIter));
        TNodeReference node = value(childIter);
        if (node.rightChild != TNode::NONE) {
            setPosition(childIter, node.rightChild);
        } else if (node.leftChild != TNode::NONE) {
            setPosition(childIter, node.leftChild);
        } else {
            SEQAN_ASSERT_MSG(false, "Unexpectedly discovered a leaf node.");
        }
    }

    // Node at the depth lcp or child of the edge at depth lcp
    TNodeReference oldChild = value(childIter);
    // Parent node of the previous node
    TNodeReference oldParent = value(parentIter);
*/

    // Insert new leaf (suffix position)
    TBufSize newLeafId = numLeaves;
    appendValue(leaves, suffix.i1);
    ++numLeaves;

    // std::cerr << std::endl << "Leaves:" << std::endl;
    // for (TSize i = 0u; i < length(leaves); ++i) {
    //     std::cerr << leaves[i] << " " <<std::endl;
    // }

    // If the last inserted suffix and the new suffix are equal there is nothing more to do
    // because we store the leftmost child in each node.
    
    //if (debug) {
        //std::cerr << "oldChild.repLength = " << oldChild.repLength << std::endl;
    //}

    //if (oldChild.repLength != suffixLength) {
        if (oldChild.repLength == lcp) {
            // Special case: path ends in a node -> new node which is a child of this node.
            // Create new node.
            TBufSize newNodeId = numNodes;
            TNode node((TBufSize) suffixLength, newLeafId, TNode::NONE, TNode::NONE);
            appendValue(nodes, node);
            ++numNodes;

            // Link to new node
            if (bitAfterLcp) {
                oldChild.rightChild = newNodeId;
            } else {
                oldChild.leftChild = newNodeId;
            }
        } else {
            // Path ends on an edge -> split edge.
            // Create new parent node which has a new child node and the old child node as children.
            TBufSize newParentId = numNodes;
            TBufSize newChildId = newParentId + 1u;
            TNode newChild(suffixLength, newLeafId, TNode::NONE, TNode::NONE);
            TNode newParent(lcp, oldChild.leftmostLeaf, oldChildId, newChildId);
            appendValue(nodes, newParent);
            ++numNodes;
            appendValue(nodes, newChild);
            ++numNodes;

            // The old parent node has to point to the new parent instead of the old child
            if (oldParent.leftChild == oldChildId) {
                oldParent.leftChild = newParentId;
            } else {
                oldParent.rightChild = newParentId;
            }
            //std::cerr << "ids: " << "op=" << oldParentId << "(" << oldParent.repLength << ")" << ", oc=" << oldChildId << "(" << oldChild.repLength << ")"  << ", np=" << newParentId << "(" << newParent.repLength << ")"  << ", nc=" << newChildId << "(" << newChild.repLength << ")"  << std::endl;
        }
    //}

    // In each step at most two nodes and exactly one leaf are added.
    return numNodes >= OUTBUF_SIZE - 1u || numLeaves >= OUTBUF_SIZE;
}

// OUTBUF_SIZE = maximal number of leaves and internal nodes of this suffix tree.
template <unsigned OUTBUF_SIZE, typename TNodes, typename TLeaves, typename TSize>
inline void
_initSuffixTree(Triple<TNodes, TLeaves, TSize> & tree, TSize offset) {
    typedef typename Value<TNodes>::Type                    TNode;
    
    reserve(tree.i1, OUTBUF_SIZE, Exact()); // nodes
    reserve(tree.i2, OUTBUF_SIZE, Exact()); // leaves
    
    // Insert root node
    appendValue(tree.i1, TNode(0u, 0u, TNode::NONE, TNode::NONE));
    tree.i3 = offset;
}

// TODO(krugel) Try other heaps
// - PriorityType<unsigned int, CompareNodeQGram2L<unsigned int, TNodes, TOccCounter>, PriorityHeap>
// - std::priority_queue<unsigned int, std::vector<unsigned int>, CompareNodeQGram2L<unsigned int, TNodes, TOccCounter> >
// - boost::mutable_queue<unsigned int, std::vector<unsigned int>, CompareNodeQGram2L<unsigned int, TNodes, TOccCounter> > 
//                                  http://www.boost.org/doc/libs/1_39_0/boost/pending/mutable_queue.hpp
// - Boost.heap:                    http://www.boost.org/doc/libs/1_49_0/doc/html/heap.html
// - boost::d_ary_heap_indirect:    http://codereview.stackexchange.com/questions/4148/slow-mutable-priority-queue
// - heapplus.h:                    http://blog.smellthedata.com/2009/08/mutable-heaps-in-c-for-updates-to.html

// Merge suffix arrays into suffix trees.
template <typename TText, typename TSpec>
inline bool
_createSuffixTrees(Index<TText, IndexDigest<TSpec> > & index) {
    typedef Index<TText, IndexDigest<TSpec> >               TIndex;
    typedef typename TIndex::TDigestSuffix                  TDigestSuffix;
    typedef typename TIndex::TSuffixArrays                  TSuffixArrays;
    typedef typename Iterator<TSuffixArrays>::Type          TSuffixArraysIterator;
    typedef typename TIndex::TSuffixArray                   TSuffixArray;
    typedef typename Iterator<TSuffixArray, Rooted>::Type   TSuffixArrayIterator;
    typedef GreaterSuffixDigest<TText, TIndex::PREFIX_LENGTH> TGreater;
    //typedef typename Reference<TSuffixArrayIterator>::Type  TSuffixReference;
    typedef typename TIndex::TSize                          TSize;
    typedef PriorityType<TSuffixArrayIterator *, TGreater>  TPriorityQueue;
    //typedef typename TIndex::TSuffixTrees                   TSuffixTrees;
    //typedef typename Reference<TSuffixTrees>::Type          TSuffixTreesReference;
    typedef typename TIndex::TSuffixTree                    TSuffixTree;
    //typedef typename TIndex::TPrefixElem                    TPrefixElem;
    typedef typename TIndex::TPosition                      TPosition;
    typedef typename TIndex::TBufSize                       TBufSize;
    typedef typename TIndex::TTextValue                     TValue;
    //typedef typename Suffix<TText>::Type                    TTextSuffix;
    typedef Lexical<LexicalSuffixDigest<True, TIndex::PREFIX_LENGTH, TSize, TText, TText> > TLexical;

    static const TBufSize OUTBUF_SIZE = TIndex::OUTBUF_SIZE;
    static const TSize INVALID_SUFFIX = MaxValue<TSize>::VALUE;
    static const TSize PREFIX_ARRAY_SIZE = TIndex::PREFIX_ARRAY_SIZE;
    ignoreUnusedVariableWarning(PREFIX_ARRAY_SIZE);

    indexRequire(index, DigestSuffixArrays());

    //std::cerr << "_createSuffixTrees" << std::endl;
    //mmapAdvise(index.suffixTrees, MAP_SEQUENTIAL);

    // Comparator for the heap
    TGreater greaterDigest(indexText(index));
    // Heap
    TPriorityQueue pq(greaterDigest);
    
    // Stores the current positions within the suffix arrays
    String<TSuffixArrayIterator> saIters;

    TSize numPartitions = length(index.suffixArrays);
    reserve(saIters, numPartitions, Exact());
    
    TSize numLeaves = 0u, numNodes = 0u;
    
    // True, if current suffix tree is full 
    bool full = true;
    // Last inserted suffix
    TDigestSuffix last;
    last.i1 = INVALID_SUFFIX;
    // Absolute end position of the text
    TPosition endPos = endPosition(indexText(index));

    // Init heap and current positions
    for (TSuffixArraysIterator sasIter = begin(index.suffixArrays); !atEnd(sasIter, index.suffixArrays); ++sasIter) {
        TSuffixArrayIterator sa = begin(*sasIter);
        appendValue(saIters, sa);
        //std::cerr << "init: " << suffix(indexText(index), (*sa).i1) << std::endl;

        SEQAN_ASSERT_LEQ(length(value(sa).i2), PREFIX_ARRAY_SIZE);
        push(pq, &back(saIters));
    }

    TSize numTrees = static_cast<TSize>(static_cast<double>(length(value(index.text))) / static_cast<double>(TIndex::OUTBUF_SIZE)) * 2u;

    TSize suffixArrayOffset = 0u;
    // Insert suffices one after another
    while (!empty(pq)) {
        // Current suffix
        TSuffixArrayIterator * currPointer = top(pq);

        // TODO(krugel) Use reference/pointer instead, also for last
        TDigestSuffix curr = value(*currPointer);
        SEQAN_ASSERT_LEQ(length(curr.i2), PREFIX_ARRAY_SIZE);

        ++(*currPointer);
        if (!atEnd(*currPointer)) {
            // Instead of pop() and push() we simply update the top element and restore the heap order
            top(pq) = currPointer;
            adjustTop(pq);
        } else {
            pop(pq);
        }

        //typedef typename Suffix<TText>::Type                    TTextSuffix;
        //TTextSuffix suffString = suffix(indexText(index), curr.i1);
        //TTextSuffix lastString = suffix(indexText(index), last.i1);
        //std::cerr << "Insert: " << suffString << std::endl;
        //if (curr.i1 == 1048571 || prefix(suffString, 7) == "abcdefg") {
        //    std::cerr << length(index.suffixTrees) << ";" << prefix(suffString, 8) << "," << prefix(lastString, 8) << std::endl;
        //    String<char> bitRep;
        //    _getBitRepresentation(bitRep, suffString);
        //    std::cerr << "Suffix: " << suffString << std::endl;
        //    std::cerr << "Bits: " << bitRep << std::endl;
        //}

        // If full create new suffix tree
        if (full) {
            if (last.i1 != INVALID_SUFFIX) {
                suffixArrayOffset += length(back(index.suffixTrees).i2);
                appendValue(index.dividers, last);
                
                //String<char> prefBits;
                //_getBitRepresentation(prefBits, last.i2);
                //std::cerr << "Divider = " << prefBits << std::endl;
                //std::cerr << "Divider = " << prefix(suffix(indexText(index), last.i1), 20u) << "... (" << length(back(index.suffixTrees).i1) << "+" << length(back(index.suffixTrees).i2) << ")" << front(back(index.suffixTrees).i2) << std::endl;
            }
            
            if (numTrees < 256
                || (numTrees <     1024 && (length(index.suffixTrees)+1) %    16 == 0)
                || (numTrees <   4*1024 && (length(index.suffixTrees)+1) %    64 == 0)
                || (numTrees <  16*1024 && (length(index.suffixTrees)+1) %   256 == 0)
                || (numTrees <  64*1024 && (length(index.suffixTrees)+1) %  1024 == 0)
                || (numTrees < 256*1024 && (length(index.suffixTrees)+1) %  4096 == 0)
                || (                       (length(index.suffixTrees)+1) % 16384 == 0)
                )
            {
                //std::cerr << length(index.suffixTrees)+1u << "/" << numTrees << ", ";
            }

            TSuffixTree newTree;
            _initSuffixTree<OUTBUF_SIZE>(newTree, suffixArrayOffset);
            last.i1 = INVALID_SUFFIX;
            appendValue(index.suffixTrees, newTree);

            //TSuffixTreesReference newTreeRef = back(index.suffixTrees);
            //_initSuffixTree<OUTBUF_SIZE>(newTreeRef, suffixArrayOffset);
            numLeaves = 0u;
            numNodes = 1u;
        }

        TSize lcp = 0u;
        bool bitAfterLcp = false;

        // Calculate lcp and first bit after lcp
        if (last.i1 == INVALID_SUFFIX) {
            bitAfterLcp = _getFirstBit(curr.i2);
        } else {
            TLexical lex(indexText(index), last, indexText(index), curr);
            bitAfterLcp = _getBitAfterLcp2(lex);
            lcp = lcpLength(lex); // counted in bit

            // Must not happen
            //if (isGreater(lex)) {
            //    compare_(lex, last, curr, true);
            //    std::cerr << "isGreater!" << std::endl;
            //    std::cerr << "last = " << last.i1 << ": " << prefix(suffix(indexText(index), last.i1), 60u) << std::endl;
            //    std::cerr << "curr = " << curr.i1 << ": " << prefix(suffix(indexText(index), curr.i1), 60u) << std::endl;
            //    std::cerr << "|trees| = " << length(index.suffixTrees) << std::endl;
            //    std::cerr << "|nodes| = " << length(back(index.suffixTrees).i1) << std::endl;
            //    std::cerr << "|leaves| = " << length(back(index.suffixTrees).i2) << std::endl;
            //}
            //SEQAN_ASSERT_MSG(isLess(lex), "The suffixes are not sorted correctly, consider increasing TAIL_LENGTH.");
            if (! isLess(lex)) {
                std::cerr << "The suffixes are not sorted correctly, consider increasing TAIL_LENGTH. The current suffix is ignored." << std::endl;
                continue;
            }
        }

        // Suffix length in bit
        TSize currLength = posSub(endPos, curr.i1) * BitsPerValue<TValue>::VALUE;

        // Insert suffix into tree
        full = _insertSuffix<OUTBUF_SIZE>(back(index.suffixTrees), curr, currLength, lcp, bitAfterLcp, numLeaves, numNodes);
        //, curr.i1 == 1048571 || prefix(suffString, 8) == "CCGGCAAA");
        //, prefix(suffString, 4) == "CTTA" && debugThis);

        last = curr;
    }

    if (last.i1 != INVALID_SUFFIX) {
        appendValue(index.dividers, last);
    }

    if (index.clearRedundantFibres) {
        // TODO(krugel) Clear
        //clear(index.suffixArrays); // Using this sometimes causes errors like:
        // "runtime error: *** glibc detected ***: double free or corruption (out)"
    }
    
    //_printSuffixTrees(index.suffixTrees);
    return true;
}

// ----------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
inline void
clear(Index<TText, IndexDigest<TSpec> > & index) {
    clear(index.partitions);
    clear(index.suffixArrays);
    clear(index.dividers);
    clear(index.suffixTrees);
}

// ----------------------------------------------------------------------------
// getFibre
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, IndexDigest<TSpec> >, DigestPartitions>::Type &
getFibre(Index<TText, IndexDigest<TSpec> > &index, DigestPartitions) {
    return index.partitions;
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, IndexDigest<TSpec> > const, DigestPartitions>::Type &
getFibre(Index<TText, IndexDigest<TSpec> > const &index, DigestPartitions) {
    return index.partitions;
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, IndexDigest<TSpec> >, DigestSuffixArrays>::Type &
getFibre(Index<TText, IndexDigest<TSpec> > &index, DigestSuffixArrays) {
    return index.suffixArrays;
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, IndexDigest<TSpec> > const, DigestSuffixArrays>::Type &
getFibre(Index<TText, IndexDigest<TSpec> > const &index, DigestSuffixArrays) {
    return index.suffixArrays;
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, IndexDigest<TSpec> >, DigestDividers>::Type &
getFibre(Index<TText, IndexDigest<TSpec> > &index, DigestDividers) {
    return index.dividers;
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, IndexDigest<TSpec> > const, DigestDividers>::Type &
getFibre(Index<TText, IndexDigest<TSpec> > const &index, DigestDividers) {
    return index.dividers;
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, IndexDigest<TSpec> >, DigestSuffixTrees>::Type &
getFibre(Index<TText, IndexDigest<TSpec> > &index, DigestSuffixTrees) {
    return index.suffixTrees;
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, IndexDigest<TSpec> > const, DigestSuffixTrees>::Type &
getFibre(Index<TText, IndexDigest<TSpec> > const &index, DigestSuffixTrees) {
    return index.suffixTrees;
}

// ----------------------------------------------------------------------------
// indexCreate
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
inline bool
indexCreate(Index<TText, IndexDigest<TSpec> > & index, DigestPartitions const, Default const) {
    return _createPartitions(index);
}

template <typename TText, typename TSpec>
inline bool
indexCreate(Index<TText, IndexDigest<TSpec> > & index, DigestSuffixArrays const, Default const) {
    // Times for 10 MB Dna
    // - Shawarma<DeepShallow> doesn't work
    // - LarssonSadakane 85 30
    // - ManberMyers > 300
    // - BWTWalk 98 22
    // - Skew3 191 71
    // - Skew7 160 62

    // TODO(krugel) Choose default based on alphabet type: Dna --> BwtWalk
    return _createSuffixArrays(index, LarssonSadakane());
}

template <typename TText, typename TSpec, typename TAlgorithm>
inline bool
indexCreate(Index<TText, IndexDigest<TSpec> > & index, DigestSuffixArrays const, TAlgorithm const) {
    return _createSuffixArrays(index, TAlgorithm());
}

template <typename TText, typename TSpec>
inline bool
indexCreate(Index<TText, IndexDigest<TSpec> > & index, DigestDividers const, Default const) {
    return _createSuffixTrees(index);
}

template <typename TText, typename TSpec>
inline bool
indexCreate(Index<TText, IndexDigest<TSpec> > & index, DigestSuffixTrees const, Default const) {
    return _createSuffixTrees(index);
}

// ----------------------------------------------------------------------------
// Save
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
inline bool
configureSave(Index<TText, IndexDigest<TSpec> > & index, const char * fileName) {
    // TODO(krugel) Add function parameter openMode
    int openMode = OPEN_RDWR | OPEN_CREATE; // We don't use OPEN_APPEND, so the file gets cleared.
    
    String<char> name;

    name = fileName; append(name, ".st");
    if (! open(getFibre(index, DigestSuffixTrees()), toCString(name), openMode)) return false;

    return true;
}

template <typename TText, typename TSpec>
inline bool
save(Index<TText, IndexDigest<TSpec> > & index, const char * fileName, int openMode) {
    String<char> name;

    name = fileName; append(name, ".txt");
    if (! save(getFibre(index, DigestText()), toCString(name), openMode)) return false;

    name = fileName; append(name, ".div");
    if (! save(getFibre(index, DigestDividers()), toCString(name), openMode)) return false;

    // We don't have to store the suffix trees, because they are in a memory-mapped string and already on disk.

    return true;
}

template <typename TText, typename TSpec>
inline bool
save(Index<TText, IndexDigest<TSpec> > & index, const char * fileName) {
    return save(index, fileName, OPEN_WRONLY | OPEN_CREATE);
}

// ----------------------------------------------------------------------------
// Open
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
inline bool
open(Index<TText, IndexDigest<TSpec> > & index, const char * fileName, int openMode) {
    String<char> name;

    name = fileName; append(name, ".txt");
    if (! open(getFibre(index, DigestText()), toCString(name), openMode)) return false;

    name = fileName; append(name, ".st");
    if (! open(getFibre(index, DigestSuffixTrees()), toCString(name), openMode)) return false;

    name = fileName; append(name, ".div");
    if (! open(getFibre(index, DigestDividers()), toCString(name), openMode)) return false;

    return true;
}

template <typename TText, typename TSpec>
inline bool
open(Index<TText, IndexDigest<TSpec> > & index, const char * fileName) {
    return open(index, fileName, OPEN_RDONLY);
}

}  // namespace seqan

#endif  // SEQAN_INDEX_INDEX_DIGEST_BASE_H_
