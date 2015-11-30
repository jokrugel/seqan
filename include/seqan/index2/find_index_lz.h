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
// Author: Tobias Stadler <tobiasstadler@mytum.de>
// ==========================================================================
// Finder class for the LZ index.
// 
// [Nav2004] G. Navarro
//   "Indexing text using the Ziv-Lempel trie",
//   Journal of Discrete Algorithms, 2004, 2, 87-114,
//   http://dx.doi.org/10.1016/S1570-8667(03)00066-2
// ==========================================================================

#ifndef SEQAN_INDEX_FIND_INDEX_LZ_H_
#define SEQAN_INDEX_FIND_INDEX_LZ_H_

namespace seqan {

template<typename TSpec = void>
struct FinderLz {};

namespace impl
{
    template<typename TValue, typename TSize, typename TBitBin, typename TPatternText>
    inline TSize findReverseText(const LZRevTrie<TValue, TSize, TBitBin> &me, const TPatternText &pattern) {
        typename Position<TPatternText>::Type i = seqan::length(pattern);

        TSize preorderId = 0;
        do {
            preorderId = getNodeChild(me, preorderId, convert<TValue>(pattern[i - 1]));

            --i;
        } while(preorderId != 0 && i > 0);

        return preorderId;
    }
}

// ============================================================================
// Metafunctions
// ============================================================================

template<typename TText, typename TBitBin, size_t NODE_NODE_END_POSITION_SAMPLE_RATE>
struct DefaultFinder<Index<TText, IndexLZ<TBitBin, NODE_NODE_END_POSITION_SAMPLE_RATE> > > {
    typedef FinderLz<> Type;
};

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template < typename TText, typename TBitBin, size_t NODE_NODE_END_POSITION_SAMPLE_RATE, typename TSpecFinder >
class Finder< Index<TText, IndexLZ<TBitBin, NODE_NODE_END_POSITION_SAMPLE_RATE> >, FinderLz<TSpecFinder>  >
{
protected:
    typedef Index<TText, IndexLZ<TBitBin, NODE_NODE_END_POSITION_SAMPLE_RATE> >    TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type                                  TSA;
    typedef typename Iterator<TSA const, Standard>::Type                           TIterator;
    typedef typename Size<TIndex>::Type                                            TSAValue;

public:
    Holder<TIndex> index;
    Pair<TIterator> range;
    TIterator data_iterator;
    TSAValue data_length;

    TSA occurrences;

    Finder()
    {
        SEQAN_CHECKPOINT

        clear(*this);
    }

    Finder(TIndex &_index):
    index(_index)
    {
        SEQAN_CHECKPOINT

        clear(*this);
    }
};

// ============================================================================
// Functions
// ============================================================================

template < typename TText, typename TBitBin, size_t LZTRIE_NODE_END_POSITION_SAMPLE_RATE, typename TSpecFinder >
inline void
clear(Finder< Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_POSITION_SAMPLE_RATE> >, FinderLz<TSpecFinder>  > & me)
{
    SEQAN_CHECKPOINT

    typedef Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_POSITION_SAMPLE_RATE> > TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type TSA;
    typedef typename Iterator<TSA const, Standard>::Type TIterator;

    me.range.i1 = me.range.i2 = TIterator();
    me.data_length = 0;

    clear(me.occurrences);
}

template<typename TText, typename TBitBin, size_t LZTRIE_NODE_END_POSITION_SAMPLE_RATE, typename TSpecFinder, typename TPatternText, typename TSpecFinder2>
inline void _findFirstIndex(Finder<Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_POSITION_SAMPLE_RATE> > , FinderLz<TSpecFinder> > &me, const TPatternText &pattern, const TSpecFinder2) {
    SEQAN_CHECKPOINT

    typedef typename Value<TText>::Type TValue;
    typedef typename SAValue<TText>::Type TSAValue;

    Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_POSITION_SAMPLE_RATE> > &index = haystack(me);
    const LZData<TValue, TSAValue, TBitBin, LZTRIE_NODE_END_POSITION_SAMPLE_RATE> &data = index.data;

    const TSAValue patternLength = length(pattern);

    indexRequire(index, FibreSA());

    if(patternLength == 0 || data.textLength < patternLength) {
        me.range.i1 = me.range.i2 = typename Iterator<typename Fibre<Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_POSITION_SAMPLE_RATE> >, FibreSA>::Type>::Type();

        return;
    }
    // find occurences of type 1
    TSAValue preorderId = findReverseText(data.reverseTrie, pattern);
    if(preorderId != 0) {
        String<TSAValue> stack1;
        reserve(stack1, impl::getNodeSubtreeSize(data.reverseTrie, preorderId), Exact());

        appendValue(stack1, preorderId, Insist());

        // for each node in the subtree of preorder id
        do {
            const TSAValue i = stack1[length(stack1) - 1], _i = getForwardTriePreorderId(data.reverseTrie, i);
            resize(stack1, length(stack1) - 1);

            // check if there is a corresponding node in LZTrie
            if(_i < length(data.forwardTrie)) {
                reserve(me.occurrences, length(me.occurrences) + impl::getNodeSubtreeSize(data.forwardTrie, _i));

                String<TSAValue> stack2;
                reserve(stack2, impl::getNodeSubtreeSize(data.forwardTrie, _i), Exact());

                appendValue(stack2, _i, Insist());

                const TSAValue blockOffset = impl::getNodeDepth(data.forwardTrie, _i) - patternLength;

                // for each node in the subtree of top(stack)
                do {
                    const TSAValue j = stack2[length(stack2) - 1];
                    resize(stack2, length(stack2) - 1);

                    const TSAValue occurencePosition = impl::getNodeEndTextPosition(data.forwardTrie, impl::getNodeId(data.forwardTrie, j) - 1) + blockOffset + 1;

                    // report the match
                    if(occurencePosition + patternLength <= data.textLength) {
                        appendValue(me.occurrences, occurencePosition);
                    }

                    TSAValue position;
                    TSAValue childId = impl::getNodeFirstChild(data.forwardTrie, j, position);
                    while(childId > 0) {
                        appendValue(stack2, childId, Insist());

                        childId = impl::getNodeNextChild(data.forwardTrie, position);
                    }
                } while(!empty(stack2));
            }

            TSAValue position;
            TSAValue childId = impl::getNodeFirstChild(data.reverseTrie, i, position);
            while(childId > 0) {
                appendValue(stack1, childId, Insist());

                childId = impl::getNodeNextChild(data.reverseTrie, position);
            }
        } while(!empty(stack1));
    }

    if(patternLength > 1) {
        // find occurences of type 2
        TSAValue suffixId = impl::getNodeChild(data.reverseTrie, (TSAValue) 0, convert<TValue>(pattern[patternLength - 1]));
        for(TSAValue i = patternLength - 1; i > 0 && suffixId != 0; --i) {
            const TSAValue suffixPreorderId = getForwardTriePreorderId(data.reverseTrie, suffixId);
            if(suffixPreorderId < length(data.forwardTrie)) {
                const TSAValue suffixSubtreeSize = impl::getNodeSubtreeSize(data.forwardTrie, suffixPreorderId);

                const TSAValue prefixId = findReverseText(data.reverseTrie, prefix(pattern, i));
                if(prefixId != 0) {
                    String<TSAValue> stack;
                    appendValue(stack, prefixId);

                    // for each node in the subtree of prefixId
                    do {
                        const TSAValue j = stack[length(stack) - 1];
                        resize(stack, length(stack) - 1, Insist());

                        const TSAValue prefixPreorderId = getForwardTriePreorderId(data.reverseTrie, j);
                        if(prefixPreorderId < length(data.forwardTrie)) {
                            const TSAValue prefixNodeId = impl::getNodeId(data.forwardTrie, prefixPreorderId);
                            const TSAValue nextPhrasePreorderId = impl::getNodePreorderId(data.forwardTrie, prefixNodeId + 1);
                            if(suffixPreorderId <= nextPhrasePreorderId && nextPhrasePreorderId < suffixPreorderId + suffixSubtreeSize) {
                                const TSAValue occurencePosition = impl::getNodeEndTextPosition(data.forwardTrie, prefixNodeId) - i + 1;

                                // report the match
                                if(occurencePosition + patternLength <= data.textLength) {
                                    appendValue(me.occurrences, occurencePosition);
                                }
                            }
                        }

                        TSAValue position;
                        TSAValue childId = impl::getNodeFirstChild(data.reverseTrie, j, position);
                        while(childId > 0) {
                            appendValue(stack, childId);

                            childId = impl::getNodeNextChild(data.reverseTrie, position);
                        }
                    } while(!empty(stack));
                }
            }

            suffixId = impl::getNodeChild(data.reverseTrie, suffixId, convert<TValue>(pattern[i - 1]));
        }

        if(patternLength > 2) {
            impl::__String<bool, TSAValue, TBitBin> used(1);
            impl::resize(used, patternLength * (patternLength + 1) / 2, Exact());

            String<TSAValue> infixPreorderIds;
            resize(infixPreorderIds, patternLength, Exact());

            // store the preorder id for the phrases infix(pattern, i, i+1)
            for(typename Size<TPatternText>::Type i = 0; i < patternLength; ++i) {
                infixPreorderIds[i] = impl::getNodeChild(data.forwardTrie, (TSAValue) 0, convert<TValue>(pattern[i]));
            }

            for(TSAValue j = 0; j < patternLength; ++j) {
                for(TSAValue i = 0; i <= j; ++i) {
                    const TSAValue S0 = impl::getNodeId(data.forwardTrie, infixPreorderIds[i]);

                    if(S0 != 0) {
                        const TSAValue usedPosition = i * patternLength - (i * (i - 1)) / 2 + (j - i);
                        if(!impl::getValue(used, usedPosition)) {
                            impl::assignValue(used, usedPosition, true);

                            TSAValue S = S0, n = j + 1, _n = n;

                            // get the longest run of connsecutive phrases that are fully contained in the pattern, starting at S0
                            preorderId = 0;
                            while(_n < patternLength) {
                                preorderId = impl::getNodeChild(data.forwardTrie, preorderId, convert<TValue>(pattern[_n]));

                                if(impl::getNodeId(data.forwardTrie, preorderId) == S + 1) {
                                    impl::assignValue(used, n * patternLength - (n * (n - 1)) / 2 + (_n - n), true);

                                    ++S;
                                    n = _n + 1;

                                    preorderId = 0;
                                } else if(preorderId == 0) {
                                    break;
                                }

                                ++_n;
                            }

                            if(preorderId != 0 && n < patternLength) {
                                // check if the the remaining characters at the end are a prefix of the successive phrase
                                TSAValue nextId = impl::getNodePreorderId(data.forwardTrie, S + 1);
                                do {
                                    if(nextId == preorderId) {
                                        ++S;
                                        n = patternLength;

                                        break;
                                    }

                                    nextId = impl::getNodeParent(data.forwardTrie, nextId);
                                } while(nextId != 0);
                            }

                            const size_t span = S - S0 + (i > 0 ? 2 : 1);
                            if(n == patternLength && span >= 3) {
                                preorderId = impl::getNodePreorderId(data.forwardTrie, S0 - 1);

                                // check if the remaining characters at the begin ar a suffix of the preceeding phrase
                                TSAValue k = i;
                                for(; preorderId != 0 && k > 0 && impl::getNodeValue(data.forwardTrie, preorderId) == pattern[k - 1]; --k) {
                                    preorderId = impl::getNodeParent(data.forwardTrie, preorderId);
                                }

                                if(k == 0) {
                                    const TSAValue occurencePosition = impl::getNodeEndTextPosition(data.forwardTrie, S0 - 1) + 1 - i;

                                    // report the match
                                    if(occurencePosition + patternLength <= data.textLength) {
                                        appendValue(me.occurrences, occurencePosition);
                                    }
                                }
                            }
                        }

                        if(j + 1 < patternLength) {
                            // store the the preorder id for the phrase infix(pattern, i, j + 1) at infixPreorderIds[i]
                            // note infixPreorderIds[i] already contains the preorder id for phrase infix(pattern, i, j)
                            infixPreorderIds[i] = impl::getNodeChild(data.forwardTrie, infixPreorderIds[i], convert<TValue>(pattern[j + 1]));
                        }
                    }
                }
            }
        }
    }

    me.range.i1 = begin(me.occurrences);
    me.range.i2 = end(me.occurrences);
}

}  // namespace seqan

#endif // SEQAN_INDEX_FIND_INDEX_LZ_H_
