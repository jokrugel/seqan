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
// Data structure and construction of the LZ index.
// 
// [Nav2004] G. Navarro
//   "Indexing text using the Ziv-Lempel trie",
//   Journal of Discrete Algorithms, 2004, 2, 87-114,
//   http://dx.doi.org/10.1016/S1570-8667(03)00066-2
// ==========================================================================

#ifndef SEQAN_INDEX_INDEX_LZ_H_
#define SEQAN_INDEX_INDEX_LZ_H_

#include <climits>
#include <cmath>
#include <cstring>

#include <algorithm>
#include <fstream>

namespace seqan {

namespace impl
{
    // ============================================================================
    // Private: Tags, Classes, Enums
    // ============================================================================

    const unsigned char RankPerByteValue[] =
    {
        8, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 4, 7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3,
        7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2,
        7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2,
        6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1,
        7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2,
        6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1,
        6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1,
        5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1, 4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0
    };

    template<typename TValue, typename TSize, typename TBitBin>
    class __String {
    public:
        unsigned char BITS_PER_VALUE;

        TSize data_capacity;
        TSize data_length;

        TBitBin *data;

        __String(const unsigned char bitsPerValue = 0):
        BITS_PER_VALUE(bitsPerValue),
        data_capacity(0),
        data_length(0),
        data(0)
        {
        }

        __String(const __String &other):
        BITS_PER_VALUE(other.BITS_PER_VALUE),
        data_capacity(0),
        data_length(0),
        data(0)
        {
            resize(*this, other.data_length, Generous());

            memcpy(data, other.data, (size_t) ceil(other.data_length * (size_t) other.BITS_PER_VALUE / (double) BitsPerValue<TBitBin>::VALUE) * BytesPerValue<TBitBin>::VALUE);
        }

        ~__String() {
            delete[] data;
        }
    };

    template<typename TKey, typename TValue, typename TBitBin>
    class HashTable {
    public:
        static const TKey PRIME1 = 4294967279u;
        static const TKey PRIME2 = 4294967197u;

        __String<TValue, TKey, TBitBin> data;

        HashTable(const unsigned char bitsPerValue = 0):
        data(bitsPerValue)
        {
        }
    };

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE = BitsPerValue<TBitBin>::VALUE>
    class ParenthesesTrie {
    public:
        __String<bool, TSize, TBitBin> parentheses;
        String<TSize> rankSamples;

        String<TValue> values;

        TSize MIN_NEAR_VALUE;
        TSize MIN_FAR_VALUE;

        HashTable<TSize, TSize, TBitBin> nearParentParenthesisPositions;
        HashTable<TSize, TSize, TBitBin> farParentParenthesisPositions;

        HashTable<TSize, TSize, TBitBin> nearClosingParenthesisPositions;
        HashTable<TSize, TSize, TBitBin> farClosingParenthesisPositions;

        ParenthesesTrie(const TSize minNearValue = 1, const TSize minFarValue = 1):
        parentheses(1),
        rankSamples(),
        values(),
        MIN_NEAR_VALUE(minNearValue),
        MIN_FAR_VALUE(minFarValue),
        nearParentParenthesisPositions(bitsPerValue(minFarValue - 1)),
        farParentParenthesisPositions(BitsPerValue<TSize>::VALUE),
        nearClosingParenthesisPositions(bitsPerValue(minFarValue - 1)),
        farClosingParenthesisPositions(BitsPerValue<TSize>::VALUE)
        {
        }
    };

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    class LZNode {
    public:
        TValue character;

        TChildrenSize children_capacity;
        TChildrenSize children_length;

        TIDValue *children;

        LZNode(const TValue &_character = TValue()):
        character(_character),
        children_capacity(0),
        children_length(0),
        children(0)
        {
        }
    };

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSITION_SAMPLE_RATE>
    class LZTrie {
    public:
        ParenthesesTrie<TValue, TSize, TBitBin> trie;

        __String<TSize, TSize, TBitBin> ids;
        __String<TSize, TSize, TBitBin> rids;

        String<TSize> endPositionSamples;

        bool hasPseudoEndingCharacter;

        LZTrie():
        trie(),
        ids(),
        rids(),
        endPositionSamples(),
        hasPseudoEndingCharacter(false)
        {
        }
    };

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    class LZRevNode {
    public:
        TIDValue forwardTriePreorderId;

        TValue character;

        TChildrenSize children_length;
        TChildrenSize children_capacity;

        TIDValue *children;

        LZRevNode(const TIDValue _forwardTriePreorderId = 0, const TValue &_character = TValue()):
        forwardTriePreorderId(_forwardTriePreorderId),
        character(_character),
        children_length(0),
        children_capacity(0),
        children(0)
        {
        }
    };

    template<typename TValue, typename TSize, typename TBitBin>
    class LZRevTrie {
    public:
        ParenthesesTrie<TValue, TSize, TBitBin> trie;

        __String<TSize, TSize, TBitBin> ids;

        LZRevTrie():
        trie(),
        ids()
        {
        }
    };

    // ============================================================================
    // Private: Functions
    // ============================================================================

    // ----------------------------------------------------------------------------
    // String operations
    // ----------------------------------------------------------------------------

    template<typename TValue, typename TSpec, typename TStream>
    inline void
    open(TStream &in, String<TValue, TSpec> &me)
    {
        typename Size<String<TValue, TSpec> >::Type data_length;
        in.read((char*) &data_length, sizeof(typename Size<String<TValue, TSpec> >::Type));

        resize(me, data_length, Exact());
        for(typename Iterator<String<TValue, TSpec> >::Type it = begin(me); it != end(me); ++it) {
            TValue _value;
            in.read((char*) &_value, sizeof(TValue));

            *it = _value;
        }
    }

    template<typename TValue, typename TSpec, typename TStream>
    inline void
    save(TStream &out, const String<TValue, TSpec> &me)
    {
        const typename Size<String<TValue, TSpec> >::Type data_length = length(me);
        out.write((const char*) &data_length, sizeof(typename Size<String<TValue, TSpec> >::Type));

        for(typename Iterator<const String<TValue, TSpec> >::Type it = begin(me); it != end(me); ++it) {
            const TValue _value = *it;

            out.write((const char*) &_value, sizeof(TValue));
        }
    }

    template<typename TValue>
    inline size_t
    size(const String<TValue> &me)
    {
        size_t _size = 0;

        _size += sizeof(TValue*); // me.data_begin
        _size += sizeof(TValue*); // me.data_end

        _size += sizeof(typename Size<String<TValue> >::Type); // me.data_capacity

        _size += me.data_capacity * sizeof(TValue);

        return _size;
    }

    // ----------------------------------------------------------------------------
    // rank, select and excess operations
    // ----------------------------------------------------------------------------

    template<typename TSize, typename TBitBin>
    inline TSize
    rank(const __String<bool, TSize, TBitBin> &me, const String<TSize> &rankSamples, const TSize rankSampleRate, const TSize _position)
    {
        SEQAN_ASSERT_LT_MSG(_position, length(me) * me.BITS_PER_VALUE, "");

        TSize i = (_position / rankSampleRate) * rankSampleRate;
        const unsigned char *byteString = getByteString(me) + i / BitsPerValue<unsigned char>::VALUE;

        TSize _rank = rankSamples[_position / rankSampleRate];
        for(; (i + BitsPerValue<unsigned char>::VALUE) < _position; i += BitsPerValue<unsigned char>::VALUE) {
            _rank += RankPerByteValue[*byteString];

            ++byteString;
        }

        for(; i < _position; ++i) {
            if(getValue(me, i) == false) {
                ++_rank;
            }
        }

        return _rank;
    }

    template<typename TSize, typename TBitBin>
    inline TSize
    select(const __String<bool, TSize, TBitBin> &me, const String<TSize> &rankSamples, const TSize rankSampleRate, TSize numberOfOccurences)
    {
        SEQAN_ASSERT_LT_MSG(numberOfOccurences, length(me) * me.BITS_PER_VALUE, "");

        TSize l = 0, r = (TSize) length(rankSamples);
        while(l < r) {
            const TSize m = l + (r - l) / 2;

            const TSize _rank = rankSamples[m];
            if(_rank < numberOfOccurences) {
                l = m + 1;
            } else if(numberOfOccurences < _rank) {
                r = m;
            } else {
                TSize _position = m * rankSampleRate - 1;
                while(getValue(me, _position) == true) {
                    --_position;
                }

                return _position;
            }
        }

        // r > 0, because me.rankSamples[0] == 0
        numberOfOccurences -= rankSamples[r - 1];

        TSize _position = (r - 1) * rankSampleRate;
        while(numberOfOccurences > 0) {
            if(getValue(me, _position) == false)
                --numberOfOccurences;

            ++_position;
        }

        return _position - 1;
    }

    template<typename TSize, typename TBitBin>
    inline TSize
    excess(const __String<bool, TSize, TBitBin> &me, const String<TSize> &rankSamples, const TSize rankSampleRate, const TSize _position)
    {
        return 2 * rank(me, rankSamples, rankSampleRate, _position) - _position;
    }

    // ----------------------------------------------------------------------------
    // __String operations
    // ----------------------------------------------------------------------------

    template<typename TValue, typename TSize, typename TBitBin>
    inline void
    clear(__String<TValue, TSize, TBitBin> &me)
    {
        me.data_length = 0;
    }

    template<typename TTargetValue, typename TSize, typename TBitBin, typename TSourceValue, typename TExpand>
    inline void
    appendValue(__String<TTargetValue, TSize, TBitBin> &me, const TSourceValue &_value, const Tag<TExpand> tag)
    {
        const TSize _position = length(me);

        resize(me, _position + 1, tag);
        assignValue(me, _position, _value);
    }

    template<typename TTargetValue, typename TSize, typename TBitBin, typename TSourceValue>
    inline void
    assignValue(__String<TTargetValue, TSize, TBitBin> &me, const TSize _position, const TSourceValue &_value)
    {
        SEQAN_ASSERT_LEQ_MSG(_position, length(me), "");
        // TODO(krugel): changed from SEQAN_ASSERT_LT_MSG to SEQAN_ASSERT_LEQ_MSG. Correct?

        TSize i = _position * me.BITS_PER_VALUE;
        for(unsigned char j = 0; j < me.BITS_PER_VALUE; ++j) {
            if(isBitSet(_value, j)) {
                setBit<TBitBin, TSize>(me.data, i + j);
            } else {
                clearBit<TBitBin, TSize>(me.data, i + j);
            }
        }
    }

    template<typename TValue, typename TSize, typename TBitBin>
    inline TValue
    getValue(const __String<TValue, TSize, TBitBin> &me, const TSize _position)
    {
        SEQAN_ASSERT_LEQ_MSG(_position, length(me), "");
        // TODO(krugel): changed from SEQAN_ASSERT_LT_MSG to SEQAN_ASSERT_LEQ_MSG. Correct?

        TValue _value = 0;

        TSize i = _position * me.BITS_PER_VALUE;
        for(unsigned char j = 0; j < me.BITS_PER_VALUE; ++j) {
            if(isBitSet<TBitBin, TSize>(me.data, i + j)) {
                setBit(_value, j);
            }
        }

        return _value;
    }

    template<typename TValue, typename TSize, typename TBitBin>
    const unsigned char*
    getByteString(const __String<TValue, TSize, TBitBin> &me)
    {
        return (unsigned char*) me.data;
    }

    template<typename TValue,typename TSize,    typename TBitBin>
    inline TSize
    capacity(const __String<TValue, TSize, TBitBin> &me)
    {
        return me.data_capacity;
    }

    template<typename TValue, typename TSize, typename TBitBin>
    inline TSize
    length(const __String<TValue, TSize, TBitBin> &me)
    {
        return me.data_length;
    }

    template<typename TValue, typename TSize, typename TBitBin,    typename TExpand>
    inline TSize
    reserve(__String<TValue, TSize, TBitBin> &me, const TSize new_capacity, const Tag<TExpand> /*tag*/)
    {
        if(capacity(me) < new_capacity) {
            const size_t new_bitbin_capacity = (TSize) ceil((new_capacity * (size_t) me.BITS_PER_VALUE) / (double) BitsPerValue<TBitBin>::VALUE);
            const size_t old_bitbin_length = (TSize) ceil((length(me) * (size_t) me.BITS_PER_VALUE) / (double) BitsPerValue<TBitBin>::VALUE);

            TBitBin *data = new TBitBin[new_bitbin_capacity];
            memset(data + old_bitbin_length, 0, (new_bitbin_capacity - old_bitbin_length) * BytesPerValue<TBitBin>::VALUE);
            memcpy(data, me.data, old_bitbin_length * BytesPerValue<TBitBin>::VALUE);

            delete[] me.data;

            me.data = data;
            me.data_capacity = new_bitbin_capacity * (size_t) BitsPerValue<TBitBin>::VALUE / me.BITS_PER_VALUE;
        }

        return capacity(me);
    }

    template<typename TValue, typename TSize, typename TBitBin, typename TExpand>
    inline TSize
    resize(__String<TValue, TSize, TBitBin> &me, const TSize new_length, const Tag<TExpand> tag)
    {
        if(capacity(me) < new_length) {
            reserve(me, new_length, tag);
        }

        me.data_length = new_length;

        return length(me);
    }

    template<typename TValue, typename TSize, typename TBitBin, typename TStream>
    inline void
    open(TStream &in, __String<TValue, TSize, TBitBin> &me)
    {
        in.read((char*) &me.BITS_PER_VALUE, sizeof(unsigned char));

        TSize data_length;
        in.read((char*) &data_length, sizeof(TSize));

        resize(me, data_length, Exact());
        for(TSize i = 0; i < (TSize) ceil(me.data_length * (size_t) me.BITS_PER_VALUE / (double) BitsPerValue<TBitBin>::VALUE); ++i) {
            in.read((char*) (me.data + i), sizeof(TBitBin));
        }
    }

    template<typename TValue, typename TSize, typename TBitBin, typename TStream>
    inline void
    save(TStream &out, const __String<TValue, TSize, TBitBin> &me)
    {
        out.write((const char*) &me.BITS_PER_VALUE, sizeof(unsigned char));

        out.write((const char*) &me.data_length, sizeof(TSize));

        for(TSize i = 0; i < (TSize) ceil(me.data_length * (size_t) me.BITS_PER_VALUE / (double) BitsPerValue<TBitBin>::VALUE); ++i) {
            out.write((const char*) (me.data + i), sizeof(TBitBin));
        }
    }

    template<typename TValue, typename TSize, typename TBitBin>
    inline size_t
    size(const __String<TValue, TSize, TBitBin> &me)
    {
        size_t _size = 0;

        _size += sizeof(unsigned char); // me.BITS_PER_VALUE

        _size += sizeof(TSize); // me.data_capacity
        _size += sizeof(TSize); // me.data_length

        _size += sizeof(TBitBin*) + (size_t) ceil(me.data_capacity * (size_t) me.BITS_PER_VALUE / (double) BitsPerValue<TBitBin>::VALUE) * sizeof(TBitBin); // me.data

        return _size;
    }

    // ----------------------------------------------------------------------------
    // HashTable operations
    // ----------------------------------------------------------------------------

    template<typename TKey, typename TValue, typename TBitBin>
    inline void
    clear(HashTable<TKey, TValue, TBitBin> &me)
    {
        clear(me.data);
    }

    template<typename TKey, typename TValue, typename TBitBin>
    inline TKey
    reserve(HashTable<TKey, TValue, TBitBin> &me, const TKey new_capacity, const Generous)
    {
        const TKey new_data_size = ((TKey) 1) << bitsPerValue(1.5 * new_capacity);

        if(length(me.data) < new_data_size) {
            resize(me.data, new_data_size, Exact());
        }

        return length(me.data);
    }

    template<typename TKey, typename TValue, typename TBitBin>
    inline void
    insertValue(HashTable<TKey, TValue, TBitBin> &me, const TKey &_key, const TValue &_value, const Insist)
    {
        const TKey bitMask = length(me.data) - 1;

        TKey _position = (_key * me.PRIME1) & bitMask;
        if (getValue(me.data, _position) != 0) {
            do {
                _position = (_position + me.PRIME2) & bitMask;
            } while(getValue(me.data, _position) != 0);
        }

        assignValue(me.data, _position, _value);
    }

    template<typename TKey, typename TValue, typename TBitBin>
    inline TValue
    getFirstValue(const HashTable<TKey, TValue, TBitBin> &me, const TKey &_key, TKey &_position)
    {
        _position = (_key * me.PRIME1) & (length(me.data) - 1);

        return getValue(me.data, _position);
    }

    template<typename TKey, typename TValue, typename TBitBin>
    inline TValue
    getNextValue(const HashTable<TKey, TValue, TBitBin> &me, TKey &_position)
    {
        _position = (_position + me.PRIME2) & (length(me.data) - 1);

        return getValue(me.data, _position);
    }

    template<typename TKey, typename TValue, typename TBitBin, typename TStream>
    inline void
    open(TStream &in, HashTable<TKey, TValue, TBitBin> &me)
    {
        open(in, me.data);
    }

    template<typename TKey, typename TValue, typename TBitBin, typename TStream>
    inline void
    save(TStream &out, const HashTable<TKey, TValue, TBitBin> &me)
    {
        save(out, me.data);
    }

    template<typename TKey, typename TValue, typename TBitBin>
    inline size_t
    size(const HashTable<TKey, TValue, TBitBin> &me)
    {
        return size(me.data);
    }

    // ----------------------------------------------------------------------------
    // ParenthesesTrie operations
    // ----------------------------------------------------------------------------

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE>
    inline void
    clear(ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me)
    {
        clear(me.parentheses);

        clear(me.rankSamples);

        clear(me.values);

        clear(me.nearClosingParenthesisPositions);
        clear(me.farClosingParenthesisPositions);
        clear(me.nearParentParenthesisPositions);
        clear(me.farParentParenthesisPositions);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE, typename TExpand>
    inline TSize
    reserve(ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me, const TSize new_capacity, const Tag<TExpand> tag)
    {
        const TSize new_parentheses_capacity = 2 * new_capacity;

        if(capacity(me.parentheses) < new_parentheses_capacity) {
            reserve(me.parentheses, new_parentheses_capacity, tag);
            reserve(me.values, new_capacity, tag);
        }

        return capacity(me.parentheses) / 2;
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE>
    inline TSize
    length(const ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me)
    {
        return length(me.parentheses) / 2;
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE>
    inline TValue
    getNodeValue(const ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me, const TSize preorderId)
    {
        return me.values[preorderId];
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE>
    inline TSize
    getNodeParent(const ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me, const TSize preorderId)
    {
        SEQAN_ASSERT_LT_MSG(preorderId, length(me), "");

        TSize _position = select(me.parentheses, me.rankSamples, (TSize) RANK_SAMPLE_RATE, preorderId + 1);

        TSize localExcess = 1, beginPosition = (me.MIN_NEAR_VALUE <= _position) ? _position - me.MIN_NEAR_VALUE : 0;
        for(TSize i = _position; i != beginPosition && beginPosition + localExcess <= i; --i) {
            if(getValue(me.parentheses, i - 1) == false) {
                --localExcess;

                if(localExcess == 0) {
                    return rank(me.parentheses, me.rankSamples, (TSize) RANK_SAMPLE_RATE, i - 1);
                }
            } else {
                ++localExcess;
            }
        }

        const TSize myExcess = excess(me.parentheses, me.rankSamples, (TSize) RANK_SAMPLE_RATE, _position) - 1;

        // check if (parent) opening paranthesis is at "near"
        TSize p, minDistance = ~(TSize) 0, distance = getFirstValue(me.nearParentParenthesisPositions, _position, p);
        while(distance != 0) {
            if((distance < minDistance || minDistance == ~(TSize) 0) && distance <= _position && excess(me.parentheses, me.rankSamples, (TSize) RANK_SAMPLE_RATE, _position - distance) == myExcess) {
                minDistance = distance;
            }

            distance = getNextValue(me.nearParentParenthesisPositions, p);
        }

        if(minDistance == ~(TSize) 0) {
            // (parent) opening paranthesis must be "far"
            distance = getFirstValue(me.farParentParenthesisPositions, _position, p);
            while(distance != 0) {
                if((distance < minDistance || minDistance == ~(TSize) 0) && distance <= _position && excess(me.parentheses, me.rankSamples, (TSize) RANK_SAMPLE_RATE, _position - distance) == myExcess) {
                    minDistance = distance;
                }

                distance = getNextValue(me.farParentParenthesisPositions, p);
            }
        }

        // return the parents preorder id
        return rank(me.parentheses, me.rankSamples, (TSize) RANK_SAMPLE_RATE, _position - minDistance);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE>
    inline TSize
    getNodeFirstChild(const ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me, const TSize preorderId, TSize &_position)
    {
        SEQAN_ASSERT_LT_MSG(preorderId, length(me), "");

        _position = select(me.parentheses, me.rankSamples, (TSize) RANK_SAMPLE_RATE, preorderId + 1) + 1;
        // check if there is a child
        if(getValue(me.parentheses, _position) == false) {
            return preorderId + 1;
        }

        return 0;
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE>
    inline TSize
    __getClosingParenthesisPosition(const ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me, const TSize _position)
    {
        TSize localExcess = 1, endPosition = _position + me.MIN_NEAR_VALUE <= length(me.parentheses) ? _position + me.MIN_NEAR_VALUE : length(me.parentheses);
        for(TSize i = _position + 1; i < endPosition && i + localExcess <= endPosition; ++i) {
            if(getValue(me.parentheses, i) == true) {
                --localExcess;

                if(localExcess == 0) {
                    return i;
                }
            } else {
                ++localExcess;
            }
        }

        const TSize myExcess = excess(me.parentheses, me.rankSamples, (TSize) RANK_SAMPLE_RATE, _position);

        // check if closing parenthesis is "near"
        TSize p, minDistance = ~(TSize) 0, distance = getFirstValue(me.nearClosingParenthesisPositions, _position, p);
        while(distance != 0) {
            if((distance < minDistance || minDistance == ~(TSize) 0) && _position + distance + 1 < length(me.parentheses) && excess(me.parentheses, me.rankSamples, (TSize) RANK_SAMPLE_RATE, _position + distance + 1) == myExcess) {
                minDistance = distance;
            }

            distance = getNextValue(me.nearClosingParenthesisPositions, p);
        }

        if(minDistance == ~(TSize) 0) {
            // closing parenthesis must be "far"
            distance = getFirstValue(me.farClosingParenthesisPositions, _position, p);
            while(distance != 0) {
                if((distance < minDistance || minDistance == ~(TSize) 0) && _position + distance + 1 < length(me.parentheses) && excess(me.parentheses, me.rankSamples, (TSize) RANK_SAMPLE_RATE, _position + distance + 1) == myExcess) {
                    minDistance = distance;
                }

                distance = getNextValue(me.farClosingParenthesisPositions, p);
            }
        }

        return _position + minDistance;
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE>
    inline TSize
    getNodeNextChild(const ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me, TSize &_position)
    {
        _position = __getClosingParenthesisPosition(me, _position) + 1;
        // check if there is a next child
        if(getValue(me.parentheses, _position) == false) {
            return rank(me.parentheses, me.rankSamples, (TSize) RANK_SAMPLE_RATE, _position);
        }

        return 0;
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE>
    inline TSize
    getNodeSubtreeSize(const ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me, const TSize preorderId) {
        const TSize _position = select(me.parentheses, me.rankSamples, (TSize) RANK_SAMPLE_RATE, preorderId + 1);

        return (__getClosingParenthesisPosition(me, _position) - _position + 1) / 2;
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE>
    inline void
    __createNodeChildrenQueryStructure(ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me, TSize &_position)
    {
        // get the closing parenthesis position
        const TSize positionOfOpeningParenthesis = _position++;
        while(getValue(me.parentheses, _position) == false) {
            __createNodeChildrenQueryStructure(me, _position);

            ++_position;
        }

        // store the closing parenthesis position
        const TSize distanceToClosingParenthesis = _position - positionOfOpeningParenthesis;
        if(me.MIN_NEAR_VALUE <= distanceToClosingParenthesis) {
            if(me.MIN_FAR_VALUE <= distanceToClosingParenthesis) {
                insertValue(me.farClosingParenthesisPositions, positionOfOpeningParenthesis, distanceToClosingParenthesis, Insist());
            } else {
                insertValue(me.nearClosingParenthesisPositions, positionOfOpeningParenthesis, distanceToClosingParenthesis, Insist());
            }
        }
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE>
    inline void
    createNodeChildrenQueryStructure(ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me, const TSize numberOfNearClosingPositions, const TSize numberOfFarClosingPositions)
    {
        me.nearClosingParenthesisPositions = HashTable<TSize, TSize, TBitBin>(bitsPerValue(me.MIN_FAR_VALUE - 1));
        reserve(me.nearClosingParenthesisPositions, numberOfNearClosingPositions, Generous());

        me.farClosingParenthesisPositions = HashTable<TSize, TSize, TBitBin>(bitsPerValue(length(me.parentheses) - 1));
        reserve(me.farClosingParenthesisPositions, numberOfFarClosingPositions, Generous());

        TSize _position = 1;
        while(getValue(me.parentheses, _position) == false) {
            __createNodeChildrenQueryStructure(me, _position);

            ++_position;
        }
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE>
    inline void
    __createNodeParentQueryStructure(ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me, TSize &_position)
    {
        const TSize positionOfParentParenthesis = _position++;
        while(getValue(me.parentheses, _position) == false) {
            // store the parent opening parenthesis position
            const TSize distanceToParentParenthesis = _position - positionOfParentParenthesis;
            if(me.MIN_NEAR_VALUE <= distanceToParentParenthesis) {
                if(me.MIN_FAR_VALUE <= distanceToParentParenthesis) {
                    insertValue(me.farParentParenthesisPositions, _position, distanceToParentParenthesis, Insist());
                } else {
                    insertValue(me.nearParentParenthesisPositions, _position, distanceToParentParenthesis, Insist());
                }
            }

            __createNodeParentQueryStructure(me, _position);

            ++_position;
        }
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE>
    inline void
    createNodeParentQueryStructure(ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me, const TSize numberOfNearParentPositions, const TSize numberOfFarParentPositions)
    {
        me.nearParentParenthesisPositions = HashTable<TSize, TSize, TBitBin>(bitsPerValue(me.MIN_FAR_VALUE - 1));
        reserve(me.nearParentParenthesisPositions, numberOfNearParentPositions, Generous());

        me.farParentParenthesisPositions = HashTable<TSize, TSize, TBitBin>(bitsPerValue(length(me.parentheses) - 1));
        reserve(me.farParentParenthesisPositions, numberOfFarParentPositions, Generous());

        TSize _position = 0;
        __createNodeParentQueryStructure(me, _position);
    }


    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE>
    inline void
    createRankQueryStructure(ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me)
    {
        reserve(me.rankSamples, (TSize) ceil(length(me.parentheses) / (double) RANK_SAMPLE_RATE), Exact());

        const unsigned char *byteString = getByteString(me.parentheses);

        TSize _rank = 0;

        appendValue(me.rankSamples, _rank, Insist());
        for(TSize i = BitsPerValue<unsigned char>::VALUE; i < length(me.parentheses); i += BitsPerValue<unsigned char>::VALUE) {
            _rank += RankPerByteValue[*byteString];

            if(i % RANK_SAMPLE_RATE == 0) {
                appendValue(me.rankSamples, _rank, Insist());
            }

            ++byteString;
        }
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE, typename TStream>
    inline void
    open(TStream &in, ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me)
    {
        open(in, me.parentheses);

        open(in, me.rankSamples);

        open(in, me.values);

        in.read((char*) &me.MIN_NEAR_VALUE, sizeof(TSize));
        in.read((char*) &me.MIN_FAR_VALUE, sizeof(TSize));

        open(in, me.nearClosingParenthesisPositions);
        open(in, me.farClosingParenthesisPositions);
        open(in, me.nearParentParenthesisPositions);
        open(in, me.farParentParenthesisPositions);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE, typename TStream>
    inline void
    save(TStream &out, const ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me)
    {
        save(out, me.parentheses);

        save(out, me.rankSamples);

        save(out, me.values);

        out.write((const char*) &me.MIN_NEAR_VALUE, sizeof(TSize));
        out.write((const char*) &me.MIN_FAR_VALUE, sizeof(TSize));

        save(out, me.nearClosingParenthesisPositions);
        save(out, me.farClosingParenthesisPositions);
        save(out, me.nearParentParenthesisPositions);
        save(out, me.farParentParenthesisPositions);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t RANK_SAMPLE_RATE>
    inline size_t
    size(const ParenthesesTrie<TValue, TSize, TBitBin, RANK_SAMPLE_RATE> &me)
    {
        size_t _size = 0;

        _size += size(me.parentheses);

        _size += size(me.rankSamples);

        _size += size(me.values);

        _size += sizeof(TSize); // TSize MIN_NEAR_VALUE
        _size += sizeof(TSize); // TSize MIN_NEAR_VALUE

        _size += size(me.nearClosingParenthesisPositions);
        _size += size(me.farClosingParenthesisPositions);
        _size += size(me.nearParentParenthesisPositions);
        _size += size(me.farParentParenthesisPositions);

        return _size;
    }

    // ----------------------------------------------------------------------------
    // LZNode operations
    // ----------------------------------------------------------------------------

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    inline TChildrenSize
    reserveChildren(LZNode<TIDValue, TValue, TChildrenSize> &me, const TChildrenSize new_children_capacity)
    {
        if(childrenCapacity(me) < new_children_capacity) {
            TIDValue *new_children = new TIDValue[new_children_capacity];
            memcpy(new_children, me.children, me.children_length * (size_t) BytesPerValue<TIDValue>::VALUE);

            delete[] me.children;

            me.children = new_children;
            me.children_capacity = new_children_capacity;
        }

        return childrenCapacity(me);
    }

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    inline void
    deleteChildren(LZNode<TIDValue, TValue, TChildrenSize> &me)
    {
        delete[] me.children;
    }

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    inline TChildrenSize
    childrenCapacity(const LZNode<TIDValue, TValue, TChildrenSize> &me)
    {
        return me.children_capacity;
    }

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    inline TChildrenSize
    childrenLength(const LZNode<TIDValue, TValue, TChildrenSize> &me)
    {
        return me.children_length;
    }

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    inline TIDValue*
    childrenBegin(const LZNode<TIDValue, TValue, TChildrenSize> &me)
    {
        return me.children;
    }

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    inline TIDValue*
    childrenEnd(const LZNode<TIDValue, TValue, TChildrenSize> &me)
    {
        return me.children + childrenLength(me);
    }

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    inline void
    appendChild(LZNode<TIDValue, TValue, TChildrenSize> &me, const TIDValue childId)
    {
        if(childrenLength(me) == childrenCapacity(me)) {
            TChildrenSize new_children_capacity = me.children_capacity > 0 ? 2 * me.children_capacity : 2;

            // there can not be more than ValueSize<TValue>::VALUE children
            //if(ValueSize<TValue>::VALUE < new_children_capacity) {
            //    new_children_capacity = (TChildrenSize) ValueSize<TValue>::VALUE;
            //}
            new_children_capacity = std::min(new_children_capacity, (TChildrenSize) ValueSize<TValue>::VALUE);

            reserveChildren(me, new_children_capacity);
        }

        me.children[me.children_length++] = childId;
    }

    // ----------------------------------------------------------------------------
    // LZTrie operations
    // ----------------------------------------------------------------------------

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline void
    clear(LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me)
    {
        clear(me.trie);

        clear(me.ids);
        clear(me.rids);

        clear(me.endPositionSamples);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline TSize
    length(const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me)
    {
        return length(me.trie);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline TSize
    getNodeId(const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me, const TSize preorderId)
    {
        return getValue(me.ids, preorderId);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline TSize
    getNodePreorderId(const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me, const TSize id)
    {
        return getValue(me.rids, id);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline TValue
    getNodeValue(const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me, const TSize preorderId)
    {
        return getNodeValue(me.trie, preorderId);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline TSize
    getNodeParent(const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me, const TSize preorderId)
    {
        return getNodeParent(me.trie, preorderId);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline TSize
    getNodeFirstChild(const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me, const TSize preorderId, TSize &_position)
    {
        return getNodeFirstChild(me.trie, preorderId, _position);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline TSize
    getNodeNextChild(const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me, TSize &_position)
    {
        return getNodeNextChild(me.trie, _position);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline TSize
    getNodeChild(const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me, const TSize preorderId, const TValue &_value)
    {
        TSize _position;
        TSize childId = getNodeFirstChild(me, preorderId, _position);
        while(childId > 0) {
            if(getNodeValue(me, childId) == _value) {
                return childId;
            }

            childId = getNodeNextChild(me, _position);
        }

        return 0;
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline TSize
    getNodeDepth(const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me, TSize preorderId)
    {
        TSize depth = 0;

        while(preorderId != 0) {
            ++depth;

            preorderId = getNodeParent(me, preorderId);
        }

        return depth;
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline TSize
    getNodeSubtreeSize(const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me, const TSize preorderId)
    {
        return getNodeSubtreeSize(me.trie, preorderId);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline TSize
    getNodeEndTextPosition(const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me, const TSize id)
    {
        SEQAN_ASSERT_LT_MSG(id, length(me), "");

        if(id % NODE_END_TEXT_POSTION_SAMPLE_RATE == 0) {
            return me.endPositionSamples[id / NODE_END_TEXT_POSTION_SAMPLE_RATE];
        }

        TSize _position = me.endPositionSamples[id / NODE_END_TEXT_POSTION_SAMPLE_RATE + 1];
        TSize nodeId = ::std::min((TSize) ((id / NODE_END_TEXT_POSTION_SAMPLE_RATE + 1) * NODE_END_TEXT_POSTION_SAMPLE_RATE), length(me) - 1);
        while(nodeId > id) {
            TSize preorderId = getNodePreorderId(me, nodeId);
            do {
                --_position;

                preorderId = getNodeParent(me, preorderId);
            } while(preorderId != 0);

            --nodeId;
        }

        return _position;
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline TSize
    getTextLength(const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me)
    {
        return getNodeEndTextPosition(me, length(me) - 1) + (!me.hasPseudoEndingCharacter ? 1 : 0);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE, typename TStream>
    inline void
    open(TStream &in, LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me)
    {
        open(in, me.trie);

        open(in, me.ids);
        open(in, me.rids);

        open(in, me.endPositionSamples);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE, typename TStream>
    inline void
    save(TStream &out, const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me)
    {
        save(out, me.trie);

        save(out, me.ids);
        save(out, me.rids);

        save(out, me.endPositionSamples);
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline size_t
    size(const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me)
    {
        return size(me.trie) + size(me.ids) + size(me.rids) + size(me.endPositionSamples);
    }

    // ----------------------------------------------------------------------------
    // LZRevNode operations
    // ----------------------------------------------------------------------------

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    inline TChildrenSize
    reserveChildren(LZRevNode<TIDValue, TValue, TChildrenSize> &me, const TChildrenSize new_children_capacity)
    {
        if(childrenCapacity(me) < new_children_capacity) {
            TIDValue *new_children = new TIDValue[new_children_capacity];
            memcpy(new_children, me.children, me.children_length * (size_t) BytesPerValue<TIDValue>::VALUE);

            delete[] me.children;

            me.children = new_children;
            me.children_capacity = new_children_capacity;
        }

        return childrenCapacity(me);
    }

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    inline void
    deleteChildren(LZRevNode<TIDValue, TValue, TChildrenSize> &me)
    {
        delete[] me.children;
    }

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    inline TChildrenSize
    childrenCapacity(const LZRevNode<TIDValue, TValue, TChildrenSize> &me)
    {
        return me.children_capacity;
    }

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    inline TChildrenSize
    childrenLength(const LZRevNode<TIDValue, TValue, TChildrenSize> &me)
    {
        return me.children_length;
    }

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    inline TIDValue*
    childrenBegin(const LZRevNode<TIDValue, TValue, TChildrenSize> &me)
    {
        return me.children;
    }

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    inline TIDValue*
    childrenEnd(const LZRevNode<TIDValue, TValue, TChildrenSize> &me)
    {
        return me.children + childrenLength(me);
    }

    template<typename TIDValue, typename TValue, typename TChildrenSize>
    inline void
    appendChild(LZRevNode<TIDValue, TValue, TChildrenSize> &me, TIDValue id)
    {
        if(childrenLength(me) == childrenCapacity(me)) {
            TChildrenSize new_children_capacity = me.children_capacity > 0 ? 2 * me.children_capacity : 2;

            // there can not be more than ValueSize<TValue>::VALUE children
            //if(ValueSize<TValue>::VALUE < new_children_capacity) {
            //    new_children_capacity = (TChildrenSize) ValueSize<TValue>::VALUE;
            //}
            new_children_capacity = std::min(new_children_capacity, (TChildrenSize) ValueSize<TValue>::VALUE);

            reserveChildren(me, new_children_capacity);
        }

        me.children[me.children_length++] = id;
    }

    // ----------------------------------------------------------------------------
    // LZRevTrie operations
    // ----------------------------------------------------------------------------

    template<typename TValue, typename TSize, typename TBitBin>
    inline void
    clear(LZRevTrie<TValue, TSize, TBitBin> &me)
    {
        clear(me.trie);
        clear(me.ids);
    }

    template<typename TValue, typename TSize, typename TBitBin>
    inline TSize
    length(const LZRevTrie<TValue, TSize, TBitBin> &me)
    {
        return length(me.trie);
    }

    template<typename TValue, typename TSize, typename TBitBin>
    inline TValue
    getNodeValue(const LZRevTrie<TValue, TSize, TBitBin> &me, const TSize preorderId)
    {
        return getNodeValue(me.trie, preorderId);
    }

    template<typename TValue, typename TSize, typename TBitBin>
    inline TSize
    getNodeFirstChild(const LZRevTrie<TValue, TSize, TBitBin> &me, const TSize preorderId, TSize &_position)
    {
        return getNodeFirstChild(me.trie, preorderId, _position);
    }

    template<typename TValue, typename TSize, typename TBitBin>
    inline TSize
    getNodeNextChild(const LZRevTrie<TValue, TSize, TBitBin> &me, TSize &_position)
    {
        return getNodeNextChild(me.trie, _position);
    }

    template<typename TValue, typename TSize, typename TBitBin>
    inline TSize
    getNodeChild(const LZRevTrie<TValue, TSize, TBitBin> &me, const TSize preorderId, const TValue &_value)
    {
        TSize _position;
        TSize childId = getNodeFirstChild(me, preorderId, _position);
        while(childId > 0) {
            if(getNodeValue(me, childId) == _value) {
                return childId;
            }

            childId = getNodeNextChild(me, _position);
        }

        return 0;
    }

    template<typename TValue, typename TSize, typename TBitBin>
    inline TSize
    getNodeSubtreeSize(const LZRevTrie<TValue, TSize, TBitBin> &me, const TSize preorderId)
    {
        return getNodeSubtreeSize(me.trie, preorderId);
    }

    template<typename TValue, typename TSize, typename TBitBin>
    inline TSize
    getForwardTriePreorderId(const LZRevTrie<TValue, TSize, TBitBin> &me, const TSize preorderId)
    {
        return getValue(me.ids, preorderId);
    }

    template<typename TValue, typename TSize, typename TBitBin, typename TStream>
    inline void
    open(TStream &in, LZRevTrie<TValue, TSize, TBitBin> &me)
    {
        open(in, me.trie);

        open(in, me.ids);
    }

    template<typename TValue, typename TSize, typename TBitBin, typename TStream>
    inline void
    save(TStream &out, const LZRevTrie<TValue, TSize, TBitBin> &me)
    {
        save(out, me.trie);

        save(out, me.ids);
    }

    template<typename TValue, typename TSize, typename TBitBin>
    inline size_t
    size(const LZRevTrie<TValue, TSize, TBitBin> &me)
    {
        return size(me.trie) + size(me.ids);
    }

    template<typename TValue, typename TSize, typename TBitBin, typename TIDValue, typename TChildrenSize>
    inline void
    __traverse(ParenthesesTrie<TValue, TSize, TBitBin> &me, const String<LZNode<TIDValue, TValue, TChildrenSize> > &trieNodes, const LZNode<TIDValue, TValue, TChildrenSize> &trieNode, TSize &numberOfNearClosingPositions, TSize &numberOfFarClosingPositions, TSize &numberOfNearParentPositions, TSize &numberOfFarParentPositions)
    {
        const TSize positionOfOpeningParenthesis = length(me.parentheses);

        appendValue(me.parentheses, false, Insist());
        appendValue(me.values, trieNode.character, Insist());

        for(TIDValue *it = childrenBegin(trieNode); it != childrenEnd(trieNode); ++it) {
            const TSize distanceToParentParenthesis = length(me.parentheses) - positionOfOpeningParenthesis;
            if(me.MIN_NEAR_VALUE <= distanceToParentParenthesis) {
                if(me.MIN_FAR_VALUE <= distanceToParentParenthesis) {
                    ++numberOfFarParentPositions;
                } else {
                    ++numberOfNearParentPositions;
                }
            }

            __traverse(me, trieNodes, trieNodes[*it], numberOfNearClosingPositions, numberOfFarClosingPositions, numberOfNearParentPositions, numberOfFarParentPositions);
        }

        appendValue(me.parentheses, true, Insist());

        const TSize distanceToClosingParenthesis = (length(me.parentheses) - 1) - positionOfOpeningParenthesis;
        if(me.MIN_NEAR_VALUE <= distanceToClosingParenthesis) {
            if(me.MIN_FAR_VALUE <= distanceToClosingParenthesis) {
                ++numberOfFarClosingPositions;
            } else {
                ++numberOfNearClosingPositions;
            }
        }
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE, typename TIDValue, typename TChildrenSize>
    inline void
    __traverse(LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me, String<LZNode<TIDValue, TValue, TChildrenSize> > &trieNodes, LZNode<TIDValue, TValue, TChildrenSize> &trieNode)
    {
        for(TIDValue *it = childrenBegin(trieNode); it != childrenEnd(trieNode); ++it) {
            appendValue(me.ids, (TSize) *it, Insist());

            __traverse(me, trieNodes, trieNodes[*it]);
        }

        deleteChildren(trieNode);
    }

    template<typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE, typename TIDValue, typename TChildrenSize, typename TText>
    inline void
    createTrie(LZTrie<typename Value<TText>::Type, typename SAValue<TText>::Type, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me, const TText &text, typename SAValue<TText>::Type &numberOfNearClosingPositions, typename SAValue<TText>::Type &numberOfFarClosingPositions, typename SAValue<TText>::Type &numberOfNearParentPositions, typename SAValue<TText>::Type &numberOfFarParentPositions)
    {
        typedef typename Value<TText>::Type TValue;
        typedef typename SAValue<TText>::Type TSAValue;

        String<LZNode<TIDValue, TValue, TChildrenSize> > trieNodes;
        reserve(trieNodes, length(text) / 10);

        typename Iterator<const TText>::Type textIt = begin(text);

        appendValue(trieNodes, LZNode<TIDValue, TValue, TChildrenSize>());
        while(textIt != end(text)) {
            TIDValue trieNodeId = 0;
            TValue character;

            // get the next phrase
            bool foundCharacter;
            do {
                character = *textIt;
                ++textIt;

                foundCharacter = false;
                for(TIDValue *it = childrenBegin(trieNodes[trieNodeId]); it != childrenEnd(trieNodes[trieNodeId]); ++it) {
                    if(character == trieNodes[*it].character) {
                        trieNodeId = *it;

                        if(textIt != end(text)) {
                            foundCharacter = true;
                        } else {
                            me.hasPseudoEndingCharacter = true;

                            character = convert<TValue>((size_t) ValueSize<TValue>::VALUE);

                            reserveChildren(trieNodes[trieNodeId], (TChildrenSize) (childrenLength(trieNodes[trieNodeId]) + 1));
                        }

                        break;
                    }
                }
            } while(foundCharacter);

            // insert the phrase into the LZTrie
            appendValue(trieNodes, LZNode<TIDValue, TValue, TChildrenSize>(character));
            appendChild(trieNodes[trieNodeId], (TIDValue) (length(trieNodes) - 1));
        }

        const TSAValue numberOfNodes = (TSAValue) length(trieNodes);

        me.trie = ParenthesesTrie<TValue, TSAValue, TBitBin>(((TSAValue) 1) << (bitsPerValue(2 * numberOfNodes) / 4), ((TSAValue) 1) << (3 * bitsPerValue(2 * numberOfNodes) / 4));
        reserve(me.trie, numberOfNodes, Exact());

        numberOfNearClosingPositions = 0;
        numberOfFarClosingPositions = 0;
        numberOfNearParentPositions = 0;
        numberOfFarParentPositions = 0;
        __traverse(me.trie, trieNodes, trieNodes[0], numberOfNearClosingPositions, numberOfFarClosingPositions, numberOfNearParentPositions, numberOfFarParentPositions);

        me.ids = __String<TSAValue, TSAValue, TBitBin>(bitsPerValue(numberOfNodes));
        reserve(me.ids, numberOfNodes, Exact());

        appendValue(me.ids, (TSAValue) 0, Insist());
        __traverse(me, trieNodes, trieNodes[0]);
    }

    template<typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE, typename TIDValue, typename TText>
    inline void
    createForwardTrie(LZTrie<typename Value<TText>::Type, typename SAValue<TText>::Type, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &me, const TText &text)
    {
        typedef typename Value<TText>::Type TValue;
        typedef typename SAValue<TText>::Type TSAValue;

        TSAValue numberOfNearClosingPositions, numberOfFarClosingPositions, numberOfNearParentPositions, numberOfFarParentPositions;

        if(ValueSize<TValue>::VALUE >= ULONG_MAX) {
            createTrie<TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE, TIDValue, unsigned long long, TText>(me, text, numberOfNearClosingPositions, numberOfFarClosingPositions, numberOfNearParentPositions, numberOfFarParentPositions);
        } else if(ValueSize<TValue>::VALUE >= UINT_MAX) {
            createTrie<TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE, TIDValue, unsigned long, TText>(me, text, numberOfNearClosingPositions, numberOfFarClosingPositions, numberOfNearParentPositions, numberOfFarParentPositions);
        } else if(ValueSize<TValue>::VALUE >= USHRT_MAX) {
            createTrie<TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE, TIDValue, unsigned int, TText>(me, text, numberOfNearClosingPositions, numberOfFarClosingPositions, numberOfNearParentPositions, numberOfFarParentPositions);
        } else if(ValueSize<TValue>::VALUE >= UCHAR_MAX) {
            createTrie<TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE, TIDValue, unsigned short, TText>(me, text, numberOfNearClosingPositions, numberOfFarClosingPositions, numberOfNearParentPositions, numberOfFarParentPositions);
        } else {
            createTrie<TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE, TIDValue, unsigned char, TText>(me, text, numberOfNearClosingPositions, numberOfFarClosingPositions, numberOfNearParentPositions, numberOfFarParentPositions);
        }

        createNodeChildrenQueryStructure(me.trie, numberOfNearClosingPositions, numberOfFarClosingPositions);
        createNodeParentQueryStructure(me.trie, numberOfNearParentPositions, numberOfFarParentPositions);
        createRankQueryStructure(me.trie);

        me.rids = __String<TSAValue, TSAValue, TBitBin>(bitsPerValue(length(me.ids)));
        resize(me.rids, length(me.ids), Exact());

        for(TSAValue i = 0; i < length(me.ids); ++i) {
            assignValue(me.rids, getValue(me.ids, i), i);
        }

        reserve(me.endPositionSamples, (typename Size<String<TSAValue> >::Type) ceil(length(me) / (double) NODE_END_TEXT_POSTION_SAMPLE_RATE) + 1, Exact());

        String<TSAValue> nodeDepths;
        resize(nodeDepths, length(me), Exact());

        nodeDepths[0] = 0;

        TSAValue endPosition = ~(TSAValue) 0;

        appendValue(me.endPositionSamples, endPosition, Insist());
        for(TSAValue i = 1; i < length(me); ++i) {
            TSAValue nodePreorderId = getNodePreorderId(me, i);
            const TSAValue nodeDepth = nodeDepths[nodePreorderId] = nodeDepths[getNodeParent(me, nodePreorderId)] + 1;

            endPosition += nodeDepth;

            if(i % NODE_END_TEXT_POSTION_SAMPLE_RATE == 0) {
                appendValue(me.endPositionSamples, endPosition, Insist());
            }
        }

        if((length(me) - 1) % NODE_END_TEXT_POSTION_SAMPLE_RATE > 0) {
            appendValue(me.endPositionSamples, endPosition, Insist());
        }
    }

    template<typename TValue, typename TSize, typename TBitBin, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE, typename TIDValue, typename TChildrenSize>
    inline void
    __traverse(const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &forwardTrie, TSize &parenthesisPosition, TSize &preorderId, String<TValue> &block, String<LZRevNode<TIDValue, TValue, TChildrenSize> > &trieNodes)
    {
        appendValue(block, getNodeValue(forwardTrie, preorderId));
        typename Position<String<TValue> >::Type blockIt = length(block);

        TIDValue trieNodeId = 0;
        TValue character;

        bool foundCharacter;
        do {
            character = block[blockIt - 1];

            foundCharacter = false;
            for(TIDValue *it = childrenBegin(trieNodes[trieNodeId]); it != childrenEnd(trieNodes[trieNodeId]); ++it) {
                if(trieNodes[*it].character == character) {
                    if(1 < blockIt) {
                        foundCharacter = true;
                    }

                    trieNodeId = *it;

                    --blockIt;

                    break;
                }
            }
        } while(foundCharacter);

        while(0 < blockIt) {
            appendValue(trieNodes, LZRevNode<TIDValue, TValue, TChildrenSize>(~(TIDValue) 0, block[blockIt - 1]));

            appendChild(trieNodes[trieNodeId], (TIDValue) (length(trieNodes) - 1));

            trieNodeId = (TIDValue) length(trieNodes) - 1;

            --blockIt;
        }

        trieNodes[trieNodeId].forwardTriePreorderId = (TIDValue) preorderId++;

        ++parenthesisPosition;
        while(getValue(forwardTrie.trie.parentheses, parenthesisPosition) == false) {
            __traverse(forwardTrie, parenthesisPosition, preorderId, block, trieNodes);
        }
        ++parenthesisPosition;

        resize(block, length(block) - 1, Exact());
    }

    template<typename TValue, typename TSize, typename TBitBin, typename TIDValue, typename TChildrenSize>
    inline void
    __traverse(ParenthesesTrie<TValue, TSize, TBitBin> &me, const String<LZRevNode<TIDValue, TValue, TChildrenSize> > &trieNodes, const LZRevNode<TIDValue, TValue, TChildrenSize> &trieNode, TSize &numberOfNearClosingPositions, TSize &numberOfFarClosingPositions)
    {
        const TSize positionOfOpeningParenthesis = length(me.parentheses);

        appendValue(me.parentheses, false, Insist());
        appendValue(me.values, trieNode.character, Insist());

        for(TIDValue *it = childrenBegin(trieNode); it != childrenEnd(trieNode); ++it) {
            __traverse(me, trieNodes, trieNodes[*it], numberOfNearClosingPositions, numberOfFarClosingPositions);
        }

        appendValue(me.parentheses, true, Insist());

        const TSize distanceToClosingParenthesis = length(me.parentheses) - 1 - positionOfOpeningParenthesis;
        if(me.MIN_NEAR_VALUE <= distanceToClosingParenthesis) {
            if(me.MIN_FAR_VALUE <= distanceToClosingParenthesis) {
                ++numberOfFarClosingPositions;
            } else {
                ++numberOfNearClosingPositions;
            }
        }
    }

    template<typename TValue, typename TSize, typename TBitBin, typename TIDValue, typename TChildrenSize>
    inline void
    __traverse(LZRevTrie<TValue, TSize, TBitBin> &me, String<LZRevNode<TIDValue, TValue, TChildrenSize> > &trieNodes, LZRevNode<TIDValue, TValue, TChildrenSize> &trieNode)
    {
        for(TIDValue *it = childrenBegin(trieNode); it != childrenEnd(trieNode); ++it) {
            appendValue(me.ids, (TSize) trieNodes[*it].forwardTriePreorderId, Insist());

            __traverse(me, trieNodes, trieNodes[*it]);
        }

        deleteChildren(trieNode);
    }

    template<typename TValue, typename TSize, typename TBitBin, typename TIDValue, typename TChildrenSize, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline void
    createTrie(LZRevTrie<TValue, TSize, TBitBin> &me, const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &forwardTrie, TSize &numberOfNearClosingPositions, TSize &numberOfFarClosingPositions)
    {
        String<LZRevNode<TIDValue, TValue, TChildrenSize> > trieNodes;
        reserve(trieNodes, (typename Size<String<LZRevNode<TIDValue, TValue, TChildrenSize> > >::Type) (1.33 * length(forwardTrie)));

        appendValue(trieNodes, LZRevNode<TIDValue, TValue, TChildrenSize>());

        reserveChildren(trieNodes[0], (TChildrenSize) (ValueSize<TValue>::VALUE + 1));

        String<TValue> block;

        TSize parenthesisPosition = 1, preorderId = 1;
        while(getValue(forwardTrie.trie.parentheses, parenthesisPosition) == false) {
            __traverse(forwardTrie, parenthesisPosition, preorderId, block, trieNodes);
        }

        const TSize numberOfNodes = (TSize) length(trieNodes);

        me.trie = ParenthesesTrie<TValue, TSize, TBitBin>(((TSize) 1) << (bitsPerValue(2 * numberOfNodes) / 4), ((TSize) 1) << (3 * bitsPerValue(2 * numberOfNodes) / 4));
        reserve(me.trie, numberOfNodes, Exact());

        numberOfNearClosingPositions = 0;
        numberOfFarClosingPositions = 0;
        __traverse(me.trie, trieNodes, trieNodes[0], numberOfNearClosingPositions, numberOfFarClosingPositions);

        me.ids = __String<TSize, TSize, TBitBin>(bitsPerValue(length(forwardTrie)));
        reserve(me.ids, numberOfNodes, Exact());

        appendValue(me.ids, (TSize) 0, Insist());
        __traverse(me, trieNodes, trieNodes[0]);
    }

    template<typename TValue, typename TSize, typename TBitBin, typename TIDValue, size_t NODE_END_TEXT_POSTION_SAMPLE_RATE>
    inline void
    createReverseTrie(LZRevTrie<TValue, TSize, TBitBin> &me, const LZTrie<TValue, TSize, TBitBin, NODE_END_TEXT_POSTION_SAMPLE_RATE> &forwardTrie)
    {
        TSize numberOfNearClosingPositions, numberOfFarClosingPositions;

        if(ValueSize<TValue>::VALUE >= ULONG_MAX) {
            impl::createTrie<TValue, TSize, TBitBin, TIDValue, unsigned long long, NODE_END_TEXT_POSTION_SAMPLE_RATE>(me, forwardTrie, numberOfNearClosingPositions, numberOfFarClosingPositions);
        } else if(ValueSize<TValue>::VALUE >= UINT_MAX) {
            impl::createTrie<TValue, TSize, TBitBin, TIDValue, unsigned long, NODE_END_TEXT_POSTION_SAMPLE_RATE>(me, forwardTrie, numberOfNearClosingPositions, numberOfFarClosingPositions);
        } else if(ValueSize<TValue>::VALUE >= USHRT_MAX) {
            impl::createTrie<TValue, TSize, TBitBin, TIDValue, unsigned int, NODE_END_TEXT_POSTION_SAMPLE_RATE>(me, forwardTrie, numberOfNearClosingPositions, numberOfFarClosingPositions);
        } else if(ValueSize<TValue>::VALUE >= UCHAR_MAX) {
            impl::createTrie<TValue, TSize, TBitBin, TIDValue, unsigned short, NODE_END_TEXT_POSTION_SAMPLE_RATE>(me, forwardTrie, numberOfNearClosingPositions, numberOfFarClosingPositions);
        } else {
            impl::createTrie<TValue, TSize, TBitBin, TIDValue, unsigned char, NODE_END_TEXT_POSTION_SAMPLE_RATE>(me, forwardTrie, numberOfNearClosingPositions, numberOfFarClosingPositions);
        }

        impl::createNodeChildrenQueryStructure(me.trie, numberOfNearClosingPositions, numberOfFarClosingPositions);
        impl::createRankQueryStructure(me.trie);
    }
}

// ============================================================================
// LZData: Tags, Classes, Enums
// ============================================================================

template<typename TValue, typename TSize, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>
class LZData {
public:
    TSize textLength;

    impl::LZTrie<TValue, TSize, TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> forwardTrie;
    impl::LZRevTrie<TValue, TSize, TBitBin> reverseTrie;

    LZData():
    textLength(0),
    forwardTrie(),
    reverseTrie()
    {
    }
};

// ============================================================================
// LZData: Functions
// ============================================================================

template<typename TValue, typename TSize, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>
inline void
clear(LZData<TValue, TSize, TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> &me)
{
    me.textLength = 0;

    clear(me.forwardTrie);
    clear(me.reverseTrie);
}

template<typename TValue, typename TSize, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE, typename TStream>
inline void
open(TStream &in, LZData<TValue, TSize, TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> &me)
{
    in.read((char*) &me.textLength, sizeof(TSize));

    impl::open(in, me.forwardTrie);
    impl::open(in, me.reverseTrie);
}

template<typename TValue, typename TSize, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE, typename TStream>
inline void
save(TStream &out, const LZData<TValue, TSize, TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> &me)
{
    out.write((const char*) &me.textLength, sizeof(TSize));

    impl::save(out, me.forwardTrie);
    impl::save(out, me.reverseTrie);
}

template<typename TValue, typename TSize, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>
inline size_t
size(const LZData<TValue, TSize, TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> &me)
{
    size_t _size = 0;

    _size += sizeof(TSize); // TSize textLength

    _size += size(me.forwardTrie);
    _size += size(me.reverseTrie);

    return _size;
}


// ============================================================================
// IndexLZ: Tags, Classes, Enums
// ============================================================================

template<typename TBitBin = size_t, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE = 128>
struct IndexLZ {};

template<typename TText, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>
class Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> > {
public:
    Holder<TText> text;

    LZData<typename Value<TText>::Type, typename SAValue<TText>::Type, TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> data;

    Index():
    text(),
    data()
    {
    }

    template<typename TOtherText>
    Index(TOtherText &_text):
    text(_text),
    data()
    {
    }
};

template<typename TText, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>
struct Fibre<Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> >, FibreSA> {
    typedef String<typename SAValue<TText>::Type> Type;
};

// ============================================================================
// IndexLZ: Functions
// ============================================================================

template<typename TText, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>
inline void
clear(Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> > &me)
{
    clear(me.text);
    clear(me.data);
}

template <typename TText, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>
inline bool
indexSupplied(Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> > &me, FibreSA const)
{
    return me.data.textLength > 0;
}

template <typename TText, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>
inline bool
indexSolveDependencies(Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> > &me, FibreSA const)
{
    return indexSupplied(me, FibreText());
}

template<typename TText, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>
inline bool
indexCreate(Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> > &me, FibreSA const)
{
    typedef typename Value<TText>::Type TValue;
    typedef typename SAValue<TText>::Type TSAValue;

    const TText &text = getFibre(me, FibreText());

    LZData<TValue, TSAValue, TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> &data = me.data;

    data.textLength = length(text);

    if(data.textLength >= ULONG_MAX) {
        impl::createForwardTrie<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE, unsigned long long, TText>(data.forwardTrie, text);
        impl::createReverseTrie<TValue, TSAValue, TBitBin, unsigned long long, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>(data.reverseTrie, data.forwardTrie);
    } else if(data.textLength >= UINT_MAX) {
        impl::createForwardTrie<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE, unsigned long, TText>(data.forwardTrie, text);
        impl::createReverseTrie<TValue, TSAValue, TBitBin, unsigned long, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>(data.reverseTrie, data.forwardTrie);
    } else if(data.textLength >= USHRT_MAX) {
        impl::createForwardTrie<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE, unsigned int, TText>(data.forwardTrie, text);
        impl::createReverseTrie<TValue, TSAValue, TBitBin, unsigned int, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>(data.reverseTrie, data.forwardTrie);
    } else if(data.textLength >= UCHAR_MAX) {
        impl::createForwardTrie<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE, unsigned short, TText>(data.forwardTrie, text);
        impl::createReverseTrie<TValue, TSAValue, TBitBin, unsigned short, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>(data.reverseTrie, data.forwardTrie);
    } else {
        impl::createForwardTrie<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE, unsigned char, const TText>(data.forwardTrie, text);
        impl::createReverseTrie<TValue, TSAValue, TBitBin, unsigned char, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>(data.reverseTrie, data.forwardTrie);
    }

    return true;
}

template<typename TReturnText, typename TText, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>
inline TReturnText
extract(const Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> > &me, const typename Position<TText>::Type _begin, const typename Position<TText>::Type _end)
{
    typedef typename Value<TText>::Type TValue;
    typedef typename SAValue<TText>::Type TSAValue;

    const LZData<TValue, TSAValue, TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> &data = me.data;

    if(data.textLength < _end || _end <= _begin) {
        return TReturnText();
    }

    const typename Size<TText>::Type textLength = _end - _begin;

    TReturnText text;
    resize(text, textLength);

    const impl::LZTrie<TValue, TSAValue, TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> &forwardTrie = data.forwardTrie;

    typename Position<String<TSAValue> >::Type l = 1, r = length(forwardTrie.endPositionSamples) - 1;
    while(l < r) {
        const typename Position<String<TSAValue> >::Type m = l + (r - l) / 2;

        if(forwardTrie.endPositionSamples[m] < _end) {
            l = m + 1;
        } else {
            r = m;
        }
    }

    TSAValue id = r < length(forwardTrie.endPositionSamples) - 1 ? (TSAValue) (r * LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE) : (impl::length(forwardTrie) - 1);
    typename Position<TText>::Type textPosition = forwardTrie.endPositionSamples[r];

    TSAValue preorderId = impl::getNodePreorderId(forwardTrie, id--);
    while(textPosition >= _end) {
        --textPosition;

        preorderId = impl::getNodeParent(forwardTrie, preorderId);
        if(preorderId == 0) {
            preorderId = impl::getNodePreorderId(forwardTrie, id--);
        }
    }

    for(typename Position<TText>::Type i = 0; i < (_end - _begin); ++i) {
        text[textLength - i - 1] = impl::getNodeValue(forwardTrie, preorderId);

        preorderId = impl::getNodeParent(forwardTrie, preorderId);
        if(preorderId == 0) {
            preorderId = impl::getNodePreorderId(forwardTrie, id--);
        }
    }

    return text;
}

template<typename TText, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>
inline typename Value<TText>::Type
extractValue(const Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> > &me, typename Position<TText>::Type _position)
{
    typedef typename Value<TText>::Type TValue;
    typedef typename SAValue<TText>::Type TSAValue;

    const LZData<TValue, TSAValue, TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> &data = me.data;

    if(data.textLength <= _position) {
        return convert<TValue>(0);
    }

    const impl::LZTrie<TValue, TSAValue, TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> &forwardTrie = data.forwardTrie;

    typename Position<String<TSAValue> >::Type l = 1, r = length(forwardTrie.endPositionSamples) - 1;
    while(l < r) {
        const typename Position<String<TSAValue> >::Type m = l + (r - l) / 2;

        if(forwardTrie.endPositionSamples[m] < _position) {
            l = m + 1;
        } else if(_position < forwardTrie.endPositionSamples[m]) {
            r = m;
        } else {
            return impl::getNodeValue(forwardTrie, impl::getNodePreorderId(forwardTrie, (TSAValue) (m * LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE)));
        }
    }

    TSAValue id = r < length(forwardTrie.endPositionSamples) - 1 ? (TSAValue) (r * LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE) : (impl::length(forwardTrie) - 1);
    typename Position<TText>::Type textPosition = forwardTrie.endPositionSamples[r];

    TSAValue preorderId = impl::getNodePreorderId(forwardTrie, id--);
    while(textPosition > _position) {
        --textPosition;

        preorderId = impl::getNodeParent(forwardTrie, preorderId);
        if(preorderId == 0) {
            preorderId = impl::getNodePreorderId(forwardTrie, id--);
        }
    }

    return impl::getNodeValue(forwardTrie, preorderId);
}

template <typename TText, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>
inline bool
open(Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> > &me, char const *filename) {
    std::ifstream in;

    in.open(filename, std::ios::in | std::ios::binary);
    if(!in.is_open()) {
        return false;
    }

    open(in, me.data);

    in.close();

    return true;
}

template <typename TText, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>
inline bool
save(const Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> > &me, char const *filename) {
    std::ofstream out;

    out.open(filename, std::ios::out | std::ios::trunc | std::ios::binary);
    if(!out.is_open()) {
        return false;
    }

    save(out, me.data);

    out.close();

    return true;
}

template<typename TText, typename TBitBin, size_t LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE>
inline size_t
size(const Index<TText, IndexLZ<TBitBin, LZTRIE_NODE_END_TEXT_POSTION_SAMPLE_RATE> > &me)
{
    return sizeof(Holder<TText>) + size(me.data);
}

}  // namespace seqan

#endif // SEQAN_INDEX_INDEX_LZ_H_
