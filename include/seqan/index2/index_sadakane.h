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
// Data structure and construction of the compressed suffix array of Sadakane.
// 
// [SS2001] K. Sadakane and T. Shibuya
//   "Indexing Huge Genome Sequences for Solving Various Problems"
//   Genome Informatics Series, 2001, 12, 175-183
// 
// [HLSS2003] W.-K. Hon, T.-W. Lam, K. Sadakane und W.-K. Sung
//   "Constructing Compressed Suffix Arrays with Large Alphabets",
//   Algorithms and Computation, 2003, 2906, 240-249,
//   http://dx.doi.org/10.1007/978-3-540-24587-2_26
// ==========================================================================

#ifndef SEQAN_INDEX_INDEX_SADAKANE_H_
#define SEQAN_INDEX_INDEX_SADAKANE_H_

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

    const unsigned char HSBPerByteValue[] =
    {
        0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8
    };

    // ============================================================================
    // Private: Functions
    // ============================================================================

    template<typename TValue, typename TSize, typename TBitBin>
    inline void
    encode(const String<TValue> &input, const TSize inputLength, TBitBin **output)
    {
        SEQAN_CHECKPOINT

        String<unsigned char> hsbPositions;
        resize(hsbPositions, inputLength, Exact());

        for(TSize i = 0; i < inputLength; ++i) {
            const unsigned char *valueByteString = (unsigned char*) &input[i];

            unsigned char j = BytesPerValue<TValue>::VALUE - 1;
            while(HSBPerByteValue[valueByteString[j]] == 0) {
                --j;
            }

            hsbPositions[i] = (unsigned char) (BitsPerValue<unsigned char>::VALUE * j) + HSBPerByteValue[valueByteString[j]];
        }

        size_t bitSize = 0;
        for(TSize i = 0; i < inputLength; ++i) {
            bitSize += 2 * hsbPositions[i] - 1;
        }

        const size_t bitBinSize = (size_t) ceil(bitSize / (double) BitsPerValue<TBitBin>::VALUE);

        *output = new TBitBin[bitBinSize];
        memset(*output, 0, bitBinSize * BytesPerValue<unsigned char>::VALUE);

        size_t bitPosition = 0;
        for(TSize i = 0; i < inputLength; ++i) {
            bitPosition += hsbPositions[i] - 1;

            setBit(*output, bitPosition++);
            for(unsigned char j = hsbPositions[i] - 1; j > 0; --j) {
                if (isBitSet(input[i], j - 1)) {
                    setBit(*output, bitPosition);
                }

                ++bitPosition;
            }
        }
    }

    template<typename TValue, typename TSize, typename TBitBin>
    inline void
    decode(const TBitBin *input, const TSize inputLength, String<TValue> &output, const TValue offsetValue)
    {
        SEQAN_CHECKPOINT

        TValue _value = offsetValue;

        size_t bitPosition = 0;
        for(TSize i = 0; i < inputLength; ++i) {
            output[i] = 0;

            const size_t beginBitPosition = bitPosition;
            while(!isBitSet(input, bitPosition)) {
                ++bitPosition;
            }

            setBit((TBitBin*) &output[i], bitPosition++ - beginBitPosition);
            for(unsigned char j = (unsigned char) (bitPosition - beginBitPosition - 1); j > 0; --j) {
                if(isBitSet(input, bitPosition++)) {
                    setBit((TBitBin*) &output[i], j - 1);
                }
            }

            _value = (output[i] += _value);
        }
    }

    template<typename TValue, typename TSize, typename TBitBin>
    inline TValue
    decodeValue(const TBitBin *input, const TSize _position)
    {
        SEQAN_CHECKPOINT

        TValue output = 0;

        size_t bitPosition = 0;
        for(TSize i = 0; i < _position; ++i) {
            TValue _value = 0;

            const size_t beginBitPosition = bitPosition;
            while(!isBitSet(input, bitPosition)) {
                ++bitPosition;
            }

            setBit((TBitBin*) &_value, bitPosition++ - beginBitPosition);
            for(unsigned char j = (unsigned char) (bitPosition - beginBitPosition - 1); j > 0; --j) {
                if(isBitSet(input, bitPosition++)) {
                    setBit((TBitBin*) &_value, j - 1);
                }
            }

            output += _value;
        }

        return output;
    }

    template<typename TTextSAValue, typename TBlockSAValue, typename TText>
    struct __PairSuffixLess :
    public ::std::binary_function < Pair<TBlockSAValue, TTextSAValue>, Pair<TBlockSAValue, TTextSAValue>, bool >
    {
        typedef typename Iterator<const TText>::Type TIterator;

        const TIterator blockBegin, blockEnd, textEnd;

        const String<TBlockSAValue> &previousInverseSuffixArray;

        __PairSuffixLess(const TText &text, const TBlockSAValue textLength, const String<TBlockSAValue> &_previousInverseSuffixArray):
        blockBegin(begin(text)),
        blockEnd(begin(text) + textLength),
        textEnd(end(text)),
        previousInverseSuffixArray(_previousInverseSuffixArray)
        {
            SEQAN_CHECKPOINT
        };

        inline bool operator() (Pair<TBlockSAValue, TTextSAValue> const &a, Pair<TBlockSAValue, TTextSAValue> const &b) const
        {
            SEQAN_CHECKPOINT

            if (a.i1 == b.i1) return false;

            if(a.i2 < b.i2) return true;
            if(b.i2 < a.i2) return false;

            TIterator itA = blockBegin + a.i1;
            TIterator itB = blockBegin + b.i1;
            if (a.i1 <= b.i1) {
                for(; itA != blockEnd && itB != textEnd; ++itA, ++itB) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }

                return itB != textEnd ? previousInverseSuffixArray[0] < previousInverseSuffixArray[itB - blockEnd] : false;
            } else {
                for(; itB != blockEnd && itA != textEnd; ++itA, ++itB) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }

                return itA != textEnd ? previousInverseSuffixArray[itA - blockEnd] < previousInverseSuffixArray[0] : true;
            }
        }
    };

    template<typename TTextSAValue, typename TBlockSAValue, typename TText>
    struct __PairSuffixLess2 :
    public ::std::binary_function < Pair<TBlockSAValue, TTextSAValue>, Pair<TBlockSAValue, TTextSAValue>, bool >
    {
        typedef typename Iterator<const TText>::Type TIterator;

        const TIterator blockBegin, blockEnd;

        const String<TBlockSAValue> &previousInverseSuffixArray;

        __PairSuffixLess2(const TText &text, const TBlockSAValue textLength, const String<TBlockSAValue> &_previousInverseSuffixArray):
        blockBegin(begin(text)),
        blockEnd(begin(text) + textLength),
        previousInverseSuffixArray(_previousInverseSuffixArray)
        {
            SEQAN_CHECKPOINT
        };

        inline bool operator() (Pair<TBlockSAValue, TTextSAValue> const &a, Pair<TBlockSAValue, TTextSAValue> const &b) const
        {
            SEQAN_CHECKPOINT

            if (a.i1 == b.i1) return false;

            if(a.i2 < b.i2) return true;
            if(b.i2 < a.i2) return false;

            TIterator itA = blockBegin + a.i1;
            TIterator itB = blockBegin + b.i1;
            if (a.i1 <= b.i1) {
                for(; itA != blockEnd; ++itA, ++itB) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }

                return previousInverseSuffixArray[0] < previousInverseSuffixArray[itB - blockEnd];
            } else {
                for(; itB != blockEnd; ++itA, ++itB) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }

                return previousInverseSuffixArray[itA - blockEnd] < previousInverseSuffixArray[0];
            }
        }
    };

    template<typename TSize, typename TBlockSAValue, typename TText>
    inline void
    increaseCharacterCount(TSize me[], const TText &text, const String<TBlockSAValue> &suffixArray)
    {
        SEQAN_CHECKPOINT

        typedef typename Value<TText>::Type TValue;

        TValue lastCharacter = convert<TValue>(0);
        for(typename Position<TText>::Type i = 0; i < length(suffixArray); ++i) {
            const TValue character = text[(typename Position<TText>::Type) suffixArray[i]];

            if(character != lastCharacter) {
                for(size_t j = ordValue(lastCharacter) + 1; j <= ordValue(character); ++j) {
                    me[j] += (TSize) i;
                }

                lastCharacter = character;
            }
        }

        for(size_t i = ordValue(lastCharacter) + 1; i < ValueSize<TValue>::VALUE + 1; ++i) {
            me[i] += (TSize) length(suffixArray);
        }
    }

    template<typename TSize, typename TTextSAValue, typename TBlockSAValue, typename TText>
    inline void
    increaseCharacterCount(TSize me[], const TText &text, const TBlockSAValue blockLength, const String<Pair<TBlockSAValue, TTextSAValue> > &suffixArray)
    {
        SEQAN_CHECKPOINT

        typedef typename Value<TText>::Type TValue;

        TValue lastCharacter = convert<TValue>(0);
        for(typename Position<TText>::Type i = 0; i < blockLength; ++i) {
            const TValue character = text[(typename Position<TText>::Type) suffixArray[i].i1];

            if (character != lastCharacter) {
                for(size_t j = ordValue(lastCharacter) + 1; j <= ordValue(character); ++j) {
                    me[j] += (TSize) i;
                }

                lastCharacter = character;
            }
        }

        for(size_t i = ordValue(lastCharacter) + 1; i < ValueSize<TValue>::VALUE + 1; ++i) {
            me[i] += (TSize) blockLength;
        }
    }

    template<typename TSAValue, typename TTextSAValue, typename TBlockSAValue>
    inline TBlockSAValue
    pMapping(const TSAValue _position, const String<Pair<TBlockSAValue, TTextSAValue> > &suffixArray, const TBlockSAValue n = 0)
    {
        SEQAN_CHECKPOINT

        typedef typename Position<const String<Pair<TBlockSAValue, TTextSAValue> > >::Type TPosition;

        TPosition l = n, r = length(suffixArray) - 1;
        while(l < r) {
            const TPosition m = l + (r - l) / 2;

            if(suffixArray[m].i2 <= _position) {
                l = m + 1;
            } else {
                r = m;
            }
        }

        return (TBlockSAValue) r;
    }
}

// ============================================================================
// SadakaneData: Tags, Classes, Enums
// ============================================================================

template<typename TValue, typename TSize, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
class SadakaneData {
public:
    TSize textLength;

    TSize *psiSamples;
    TBitBin **compressedPsi;

    TSize *I;
    TSize *J;

    TSize C[ValueSize<TValue>::VALUE + 1];

    SadakaneData():
    textLength(0),
    psiSamples(0),
    compressedPsi(0),
    I(0),
    J(0)
    {
        SEQAN_CHECKPOINT

        arrayFill(C, C + ValueSize<TValue>::VALUE + 1, 0);
    }
};

// ============================================================================
// SadakaneData: Functions
// ============================================================================

template<typename TValue, typename TSize, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
inline void
clear(SadakaneData<TValue, TSize, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> &me)
{
    SEQAN_CHECKPOINT

    delete[] me.psiSamples;
    me.psiSamples = 0;

    if(PSI_SAMPLE_RATE > 1) {
        for(TSize i = 0; i < (me.textLength - 1) / PSI_SAMPLE_RATE + 1; ++i) {
            delete[] me.compressedPsi[i];
        }

        delete[] me.compressedPsi;
        me.compressedPsi = 0;
    }


    delete[] me.I;
    me.I = 0;

    delete[] me.J;
    me.J = 0;

    arrayFill(me.C, me.C + ValueSize<TValue>::VALUE + 1, 0);

    me.textLength = 0;
}

template<typename TValue, typename TSize, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
inline TSize
__psiAt(const SadakaneData<TValue, TSize, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> &me, const TSize _position)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT_LEQ_MSG(_position, me.textLength, "");

    TSize psi = me.psiSamples[_position / PSI_SAMPLE_RATE];

    if (_position % PSI_SAMPLE_RATE > 0) {
        psi += impl::decodeValue<TSize, TSize, TBitBin>(me.compressedPsi[_position / PSI_SAMPLE_RATE], _position % PSI_SAMPLE_RATE);
    }

    return psi;
}

template<typename TValue, typename TSize, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
inline TSize
__saAt(const SadakaneData<TValue, TSize, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> &me, TSize _position)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT_LEQ_MSG(_position, me.textLength, "");

    TSize i = 0;
    while(_position % SA_SAMPLE_RATE > 0) {
        _position = __psiAt(me, _position);

        ++i;
    }

    return me.I[_position / SA_SAMPLE_RATE] - i;
}

template<typename TValue, typename TSize, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
inline TSize
__isaAt(const SadakaneData<TValue, TSize, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> &me, const TSize _position)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT_LEQ_MSG(_position, me.textLength, "");

    TSize i = me.J[_position / ISA_SAMPLE_RATE];
    for(TSize j = 0; j < _position % ISA_SAMPLE_RATE; ++j) {
        i = __psiAt(me, i);
    }

    return i;
}

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE, typename TTextSAValue>
inline void
create(SadakaneData<typename Value<TText>::Type, typename SAValue<TText>::Type, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> &me, const TText &text)
{
    SEQAN_CHECKPOINT

    typedef typename Value<TText>::Type TValue;
    typedef typename SAValue<TText>::Type TSAValue;

    arrayFill(me.C, me.C + ValueSize<TValue>::VALUE + 1, 1);

    me.textLength = length(text);

    String<TSAValue> suffixArray;
    resize(suffixArray, me.textLength, Exact());
    // create the suffix array
    if(me.textLength > 6) {
        createSuffixArray(suffixArray, text, Skew7());
    } else if(me.textLength > 1) {
        createSuffixArray(suffixArray, text, Skew3());
    } else {
        suffixArray[0] = 0;
    }

    me.I = new TSAValue[(me.textLength - 1) / SA_SAMPLE_RATE + 1];

    // sample the suffix array, sample length is SA_SAMPLE_RATE
    me.I[0] = me.textLength;
    for(TSAValue i = 1; i < (me.textLength - 1) / SA_SAMPLE_RATE + 1; ++i) {
        me.I[i] = suffixArray[i * SA_SAMPLE_RATE - 1];
    }

    String<TTextSAValue> inverseSuffixArray;
    resize(inverseSuffixArray, me.textLength, Exact());

    // create the inverse suffix array
    for(TSAValue i = 0; i < me.textLength; ++i) {
        inverseSuffixArray[suffixArray[i]] = (TTextSAValue) i;
    }

    me.J = new TSAValue[(me.textLength - 1) / ISA_SAMPLE_RATE + 1];

    // sample the inverse suffix array, sample length is ISA_SAMPLE_RATE
    for(TSAValue i = 0; i < (me.textLength - 1) / ISA_SAMPLE_RATE + 1; ++i) {
        me.J[i] = (TSAValue) inverseSuffixArray[i * ISA_SAMPLE_RATE] + 1;
    }

    me.psiSamples = new TSAValue[(me.textLength - 1) / PSI_SAMPLE_RATE + 1];
    if(PSI_SAMPLE_RATE > 1) {
        me.compressedPsi = new TBitBin*[(size_t) ceil(me.textLength / (double) PSI_SAMPLE_RATE)];
    }

    String<TSAValue> psiDiffs;
    resize(psiDiffs, PSI_SAMPLE_RATE - 1, Exact());

    // calculate psi and sample it, sample length is PSI_SAMPLE_RATE
    TSAValue previousPsiValue = me.psiSamples[0] = inverseSuffixArray[0] + 1;
    for(TSAValue i = 1; i <= me.textLength; ++i) {
        const TSAValue psiValue = suffixArray[i - 1] < me.textLength - 1 ? inverseSuffixArray[suffixArray[i - 1] + 1] + 1 : 0;

        if(i % PSI_SAMPLE_RATE == 0) {
            me.psiSamples[i / PSI_SAMPLE_RATE] = psiValue;

            // compress the values between two psi samples only if there are such values
            if(PSI_SAMPLE_RATE > 1) {
                impl::encode(psiDiffs, PSI_SAMPLE_RATE - 1, &me.compressedPsi[i / PSI_SAMPLE_RATE - 1]);
            }
        } else {
            psiDiffs[i % PSI_SAMPLE_RATE - 1] = (TSAValue) (psiValue - previousPsiValue);
        }

        previousPsiValue = psiValue;
    }

    // compress the values after the last psi sample
    if(me.textLength % PSI_SAMPLE_RATE > 0) {
        impl::encode(psiDiffs, me.textLength % PSI_SAMPLE_RATE, &me.compressedPsi[(me.textLength - 1) / PSI_SAMPLE_RATE]);
    }

    // adjust C
    impl::increaseCharacterCount(me.C, text, suffixArray);
}

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE, typename TTextSAValue, typename TBlockSAValue>
inline void
create(SadakaneData<typename Value<TText>::Type, typename SAValue<TText>::Type, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> &me, const TText &text, TBlockSAValue blockLength)
{
    SEQAN_CHECKPOINT

    typedef typename Value<TText>::Type TValue;
    typedef typename SAValue<TText>::Type TSAValue;

    arrayFill(me.C, me.C + ValueSize<TValue>::VALUE + 1, 1);

    me.textLength = length(text);

    const TTextSAValue numberOfBlocks = (TTextSAValue) ceil(me.textLength / (double) blockLength);

    // calculate the length of the last block
    TTextSAValue currentSuffixLength = (TTextSAValue) me.textLength % blockLength;
    if(currentSuffixLength == 0) {
        currentSuffixLength = (TTextSAValue) blockLength;
    }

    typename Suffix<const TText>::Type currentSuffix = suffix(text, me.textLength - currentSuffixLength);

    String<TSAValue> initialSuffixArray;
    resize(initialSuffixArray, currentSuffixLength, Exact());

    // create the suffix array for the last block
    if(currentSuffixLength > 6) {
        createSuffixArray(initialSuffixArray, currentSuffix, Skew7());
    } else if(currentSuffixLength > 1) {
        createSuffixArray(initialSuffixArray, currentSuffix, Skew3());
    } else {
            initialSuffixArray[0] = 0;
    }

    String<TBlockSAValue> inverseSuffixArray;
    resize(inverseSuffixArray, blockLength + 1, Exact());

    // create the inverse suffix arary for the last block
    for(TBlockSAValue i = 0; i < currentSuffixLength; ++i) {
        inverseSuffixArray[initialSuffixArray[i]] = i;
    }

    me.psiSamples = new TSAValue[(me.textLength - 1) / PSI_SAMPLE_RATE + 1];

    if(PSI_SAMPLE_RATE > 1) {
        me.compressedPsi = new TBitBin*[(size_t) ceil(me.textLength / (double) PSI_SAMPLE_RATE)];
    }

    String<TSAValue> psiDiffs;
    resize(psiDiffs, PSI_SAMPLE_RATE - 1, Exact());

    // calculate psi for currentSuffix and sample it, sample length is PSI_SAMPLE_RATE
    TSAValue previousPsiValue = me.psiSamples[0] = (TSAValue) inverseSuffixArray[0] + 1;
    for(TTextSAValue i = 1; i <= currentSuffixLength; ++i) {
        const TSAValue psiValue = (initialSuffixArray[i - 1] < (TSAValue) (currentSuffixLength - 1)) ? (TSAValue) inverseSuffixArray[initialSuffixArray[i - 1] + 1] + 1 : 0;

        if(i % PSI_SAMPLE_RATE == 0) {
            me.psiSamples[i / PSI_SAMPLE_RATE] = psiValue;

            // compress the values between two psi samples only if there are such values
            if(PSI_SAMPLE_RATE > 1) {
                impl::encode(psiDiffs, PSI_SAMPLE_RATE - 1, &me.compressedPsi[i / PSI_SAMPLE_RATE - 1]);
            }
        } else {
            psiDiffs[i % PSI_SAMPLE_RATE - 1] = (TSAValue) (psiValue - previousPsiValue);
        }

        previousPsiValue = psiValue;
    }

    // compress the values after the last psi sample
    if(currentSuffixLength % PSI_SAMPLE_RATE > 0) {
        impl::encode(psiDiffs, currentSuffixLength % PSI_SAMPLE_RATE, &me.compressedPsi[currentSuffixLength / PSI_SAMPLE_RATE]);
    }

    // adjust C
    impl::increaseCharacterCount(me.C, currentSuffix, initialSuffixArray);

    String<TSAValue> psiValues;
    resize(psiValues, PSI_SAMPLE_RATE - 1, Exact());

    String<Pair<TBlockSAValue, TTextSAValue> > suffixArray;
    resize(suffixArray, blockLength + 1, Exact());

    String<Pair<TTextSAValue, TTextSAValue> > psiBuffer;
    resize(psiBuffer, (typename Size<String<Pair<TTextSAValue, TTextSAValue> > >::Type) ceil((blockLength + 1) / (double) PSI_SAMPLE_RATE) * PSI_SAMPLE_RATE - 1, Exact());

    for(TTextSAValue i = 1; i < numberOfBlocks; ++i) {
        currentSuffixLength += (TTextSAValue) blockLength;

        currentSuffix = suffix(text, me.textLength - currentSuffixLength);

        // calculate the suffix array position for each of the blockLength longest suffixes and store it in .i2
        suffixArray[blockLength] = Pair<TBlockSAValue, TTextSAValue>(blockLength, (TTextSAValue) me.psiSamples[0]);
        for(TBlockSAValue j = blockLength; j > 0; --j) {
            const size_t character = ordValue(currentSuffix[j - 1]);

            suffixArray[j - 1] = Pair<TBlockSAValue, TTextSAValue>(j - 1, (TTextSAValue) me.C[character]);

            if(me.C[character] < me.C[character + 1]) {
                TSAValue l = me.C[character] / PSI_SAMPLE_RATE + 1, r = (me.C[character + 1] - 1) / PSI_SAMPLE_RATE;
                if(l <= r) {
                    if(suffixArray[j].i2 <= me.psiSamples[l]) {
                        if(me.C[character] % PSI_SAMPLE_RATE > 0 || me.psiSamples[l - 1] < suffixArray[j].i2) {
                            const size_t decodingStart = me.C[character] % PSI_SAMPLE_RATE - ((me.C[character] % PSI_SAMPLE_RATE > 0) ? 1 : 0);

                            impl::decode(me.compressedPsi[l - 1], PSI_SAMPLE_RATE - 1, psiValues, me.psiSamples[l - 1]);

                            if(psiValues[decodingStart] < suffixArray[j].i2) {
                                suffixArray[j - 1].i2 = (l - 1) * PSI_SAMPLE_RATE + (std::lower_bound(begin(psiValues) + decodingStart, begin(psiValues) + (PSI_SAMPLE_RATE - 1), suffixArray[j].i2) - begin(psiValues)) + 1;
                            } else if(me.C[character] % PSI_SAMPLE_RATE == 0) {
                                ++suffixArray[j - 1].i2;
                            }
                        }
                    } else if(me.psiSamples[r] < suffixArray[j].i2) {
                        const size_t decodingStop = (me.C[character + 1] - 1) % PSI_SAMPLE_RATE;

                        impl::decode(me.compressedPsi[r], decodingStop, psiValues, me.psiSamples[r]);

                        suffixArray[j - 1].i2 = r * PSI_SAMPLE_RATE + (std::lower_bound(begin(psiValues), begin(psiValues) + decodingStop, suffixArray[j].i2) - begin(psiValues)) + 1;
                    } else {
                        const TSAValue *p = std::lower_bound(me.psiSamples + l, me.psiSamples + r, suffixArray[j].i2);
                        if(suffixArray[j].i2 < *p) {
                            impl::decode(me.compressedPsi[p - me.psiSamples - 1], PSI_SAMPLE_RATE - 1, psiValues, me.psiSamples[p - me.psiSamples - 1]);

                            suffixArray[j - 1].i2 = (p - me.psiSamples - 1) * PSI_SAMPLE_RATE + (std::lower_bound(begin(psiValues), begin(psiValues) + PSI_SAMPLE_RATE - 1, suffixArray[j].i2) - begin(psiValues)) + 1;
                        } else {
                            suffixArray[j - 1].i2 = (p - me.psiSamples) * PSI_SAMPLE_RATE;
                        }
                    }
                } else {
                    if(me.C[character] % PSI_SAMPLE_RATE > 0 || me.psiSamples[l - 1] < suffixArray[j].i2) {
                        const size_t decodingStart = me.C[character] % PSI_SAMPLE_RATE - ((me.C[character] % PSI_SAMPLE_RATE > 0) ? 1 : 0);
                        const size_t decodingStop = (me.C[character + 1] - 1) % PSI_SAMPLE_RATE;

                        impl::decode(me.compressedPsi[l - 1], decodingStop, psiValues, me.psiSamples[l - 1]);

                        if(psiValues[decodingStart] < suffixArray[j].i2) {
                            suffixArray[j - 1].i2 = (l - 1) * PSI_SAMPLE_RATE + (std::lower_bound(begin(psiValues) + decodingStart, begin(psiValues) + decodingStop, suffixArray[j].i2) - begin(psiValues)) + 1;
                        } else if(me.C[character] % PSI_SAMPLE_RATE == 0) {
                            ++suffixArray[j - 1].i2;
                        }
                    }
                }
            }
        }

        // sort the blockLength longest suffixes
        if(currentSuffixLength < 2 * blockLength) {
            std::sort(begin(suffixArray), begin(suffixArray) + blockLength, impl::__PairSuffixLess<TTextSAValue, TBlockSAValue, typename Suffix<const TText>::Type const>(currentSuffix, blockLength, inverseSuffixArray));
        } else {
            std::sort(begin(suffixArray), begin(suffixArray) + blockLength, impl::__PairSuffixLess2<TTextSAValue, TBlockSAValue, typename Suffix<const TText>::Type const>(currentSuffix, blockLength, inverseSuffixArray));
        }

        suffixArray[blockLength].i2 += impl::pMapping(suffixArray[blockLength].i2, suffixArray);

        for(TBlockSAValue j = 0; j <= blockLength; ++j) {
            inverseSuffixArray[suffixArray[j].i1] = j;
        }

        TTextSAValue jf = 1;
        TBlockSAValue n = impl::pMapping(jf, suffixArray);

        resize(psiBuffer, ::std::min(capacity(psiBuffer), (typename Size<String<Pair<TTextSAValue, TTextSAValue> > >::Type) (currentSuffixLength - blockLength + 1)), Insist());

        // in order to avoid having the psi from the old iteration and psi from current iteration in memory, buffer the old psi values and reuse psiSamples and compressedPsi
        // thus fill the the buffer first
        typename Size<String<Pair<TTextSAValue, TTextSAValue> > >::Type k = 0;
        for(typename Size<String<Pair<TTextSAValue, TTextSAValue> > >::Type j = 0; j < length(psiBuffer) / PSI_SAMPLE_RATE; ++j) {
            impl::decode(me.compressedPsi[j], PSI_SAMPLE_RATE - 1, psiValues, me.psiSamples[j]);

            for(size_t l = 0; l < PSI_SAMPLE_RATE - 1; ++l){
                psiBuffer[k++] = Pair<TTextSAValue, TTextSAValue>(jf + n, (TTextSAValue) (psiValues[l] + impl::pMapping(psiValues[l], suffixArray)));

                n = impl::pMapping(++jf, suffixArray, n);
            }

            psiBuffer[k++] = Pair<TTextSAValue, TTextSAValue>(jf + n, (TTextSAValue) (me.psiSamples[j + 1] + impl::pMapping(me.psiSamples[j + 1], suffixArray)));

            n = impl::pMapping(++jf, suffixArray, n);
        }

        if(jf <= currentSuffixLength - blockLength) {
            impl::decode(me.compressedPsi[length(psiBuffer) / PSI_SAMPLE_RATE], ::std::min(PSI_SAMPLE_RATE - 1,
              (size_t) (currentSuffixLength - blockLength - jf + 1)), psiValues, me.psiSamples[length(psiBuffer) / PSI_SAMPLE_RATE]);

            for(size_t l = 0; l < length(psiBuffer) % PSI_SAMPLE_RATE; ++l){
                psiBuffer[k++] = Pair<TTextSAValue, TTextSAValue>(jf + n, (TTextSAValue) (psiValues[l] + impl::pMapping(psiValues[l], suffixArray)));

                n = impl::pMapping(++jf, suffixArray, n);
            }
        }

        const size_t maxCompressedPsi = (size_t) ceil((currentSuffixLength - blockLength) / (double) PSI_SAMPLE_RATE);

        typename Position<String<Pair<TTextSAValue, TTextSAValue> > >::Type psiBufferIt = 0;
        typename Iterator<const String<Pair<TBlockSAValue, TTextSAValue> > >::Type saIt = begin(suffixArray);

        // create psi for currentSuffix and sample it, sample length is PSI_SAMPLE_RATE
        TSAValue previousPsiValue = me.psiSamples[0] = (TSAValue) (suffixArray[inverseSuffixArray[0]].i2 + inverseSuffixArray[0]);
        for(TTextSAValue j = 1; j <= currentSuffixLength; ++j) {
            TSAValue psiValue;
            if(j == psiBuffer[psiBufferIt].i1) {
                psiValue = (TSAValue) psiBuffer[psiBufferIt].i2;

                // refill the buffer, only if there are any psi values left
                if(jf <= currentSuffixLength - blockLength) {
                    if(jf % PSI_SAMPLE_RATE == 0) {
                        psiBuffer[psiBufferIt] = Pair<TTextSAValue, TTextSAValue>(jf + n, (TTextSAValue) (me.psiSamples[jf / PSI_SAMPLE_RATE] + impl::pMapping(me.psiSamples[jf / PSI_SAMPLE_RATE], suffixArray)));

                        const size_t decodingStop = ::std::min(PSI_SAMPLE_RATE - 1, (size_t) (currentSuffixLength - blockLength - jf));

                        impl::decode(me.compressedPsi[jf / PSI_SAMPLE_RATE], decodingStop, psiValues, me.psiSamples[jf / PSI_SAMPLE_RATE]);
                    } else {
                        psiBuffer[psiBufferIt] = Pair<TTextSAValue, TTextSAValue>(jf + n, (TTextSAValue) (psiValues[jf % PSI_SAMPLE_RATE - 1] + impl::pMapping(psiValues[jf % PSI_SAMPLE_RATE - 1], suffixArray)));
                    }

                    n = impl::pMapping(++jf, suffixArray, n);
                }

                psiBufferIt = (psiBufferIt + 1) % length(psiBuffer);
            } else {
                psiValue = suffixArray[inverseSuffixArray[saIt->i1 + 1]].i2 + (inverseSuffixArray[saIt->i1 + 1] < blockLength ? (TSAValue) inverseSuffixArray[saIt->i1 + 1] : 0);

                ++saIt;
            }

            if(j % PSI_SAMPLE_RATE == 0) {
                me.psiSamples[j / PSI_SAMPLE_RATE] = psiValue;

                // compress the values between two psi samples only if there are such values
                if(PSI_SAMPLE_RATE > 1) {
                    // delete the old values, if ther are one
                    if(j / PSI_SAMPLE_RATE - 1 < maxCompressedPsi) {
                        delete[] me.compressedPsi[j / PSI_SAMPLE_RATE - 1];
                    }

                    impl::encode(psiDiffs, PSI_SAMPLE_RATE - 1, &me.compressedPsi[j / PSI_SAMPLE_RATE - 1]);
                }
            } else {
                psiDiffs[j % PSI_SAMPLE_RATE - 1] = (TSAValue) (psiValue - previousPsiValue);
            }

            previousPsiValue = psiValue;
        }

        // compress the values after the last psi sample
        if(currentSuffixLength % PSI_SAMPLE_RATE > 0) {
            // delete the old values, if ther are one
            if(currentSuffixLength / PSI_SAMPLE_RATE < maxCompressedPsi) {
                delete[] me.compressedPsi[currentSuffixLength / PSI_SAMPLE_RATE];
            }

            impl::encode(psiDiffs, currentSuffixLength % PSI_SAMPLE_RATE, &me.compressedPsi[currentSuffixLength / PSI_SAMPLE_RATE]);
        }

        // adjust C
        impl::increaseCharacterCount(me.C, currentSuffix, blockLength, suffixArray);
    }

    me.I = new TSAValue[(me.textLength - 1) / SA_SAMPLE_RATE + 1];
    me.J = new TSAValue[(me.textLength - 1) / ISA_SAMPLE_RATE + 1];

    TSAValue saValue = me.J[0] = me.psiSamples[0];
    if(saValue % SA_SAMPLE_RATE == 0) {
        me.I[saValue / SA_SAMPLE_RATE] = 0;
    }

    for(TTextSAValue i = 1; i <= me.textLength; ++i) {
        saValue = __psiAt(me, saValue);

        if(saValue % SA_SAMPLE_RATE == 0) {
            me.I[saValue / SA_SAMPLE_RATE] = (TSAValue) i;
        }

        if(i % ISA_SAMPLE_RATE == 0) {
            me.J[i / ISA_SAMPLE_RATE] = saValue;
        }
    }
}

template<typename TValue, typename TSize, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE, typename TStream>
inline void
open(TStream &in, SadakaneData<TValue, TSize, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> &me)
{
    SEQAN_CHECKPOINT

    in.read((char*) &me.textLength, sizeof(TSize));

    me.psiSamples = new TSize[(me.textLength - 1) / PSI_SAMPLE_RATE + 1];
    in.read((char*) me.psiSamples, ((me.textLength - 1) / PSI_SAMPLE_RATE + 1) * sizeof(TSize));

    if(PSI_SAMPLE_RATE > 1) {
        me.compressedPsi = new TBitBin*[(me.textLength - 1) / PSI_SAMPLE_RATE + 1];
        for(TSize i = 0; i < (me.textLength - 1) / PSI_SAMPLE_RATE + 1; ++i) {
            TSize row_length;
            in.read((char*) &row_length, sizeof(TSize));

            me.compressedPsi[i] = new TBitBin[row_length];
            in.read((char*) me.compressedPsi[i], row_length * sizeof(TBitBin));
        }
    }

    me.I = new TSize[(me.textLength - 1) / SA_SAMPLE_RATE + 1];
    in.read((char*) me.I, ((me.textLength - 1)/ SA_SAMPLE_RATE + 1) * sizeof(TSize));

    me.J = new TSize[(me.textLength - 1) / ISA_SAMPLE_RATE + 1];
    in.read((char*) me.J, ((me.textLength - 1) / ISA_SAMPLE_RATE + 1) * sizeof(TSize));

    in.read((char*) me.C, (ValueSize<TValue>::VALUE + 1) * sizeof(TSize));
}

template<typename TValue, typename TSize, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE, typename TStream>
inline void
save(TStream &out, const SadakaneData<TValue, TSize, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> &me)
{
    SEQAN_CHECKPOINT

    out.write((const char*) &me.textLength, sizeof(TSize));

    out.write((const char*) me.psiSamples, ((me.textLength - 1) / PSI_SAMPLE_RATE + 1) * sizeof(TSize));

    if(PSI_SAMPLE_RATE > 1) {
        for(TSize i = 0; i < (me.textLength - 1) / PSI_SAMPLE_RATE; ++i) {
            const TBitBin *row = me.compressedPsi[i];

            size_t position = 0;
            for(TSize j = 0; j < PSI_SAMPLE_RATE - 1; ++j) {
                const size_t startPosition = position;
                while(!impl::isBitSet(row, position)) {
                    ++position;
                }

                position += position - startPosition + 1;
            }

            const TSize row_length = (TSize) ceil(position / (double) BitsPerValue<TBitBin>::VALUE);
            out.write((const char*) &row_length, sizeof(TSize));

            out.write((const char*) row, row_length * sizeof(TBitBin));
        }

        const TBitBin *row = me.compressedPsi[(me.textLength - 1) / PSI_SAMPLE_RATE];

        size_t position = 0;
        for(TSize j = 0; j < ::std::min((TSize) PSI_SAMPLE_RATE, (TSize) (me.textLength % PSI_SAMPLE_RATE) - 1); ++j) {
            const size_t startPosition = position;
            while(!impl::isBitSet(row, position)) {
                ++position;
            }

            position += position - startPosition + 1;
        }

        const TSize row_length = (TSize) ceil(position / (double) BitsPerValue<TBitBin>::VALUE);
        out.write((const char*) &row_length, sizeof(TSize));

        out.write((const char*) row, row_length * sizeof(TBitBin));
    }

    out.write((const char*) me.I, ((me.textLength - 1) / SA_SAMPLE_RATE + 1) * sizeof(TSize));

    out.write((const char*) me.J, ((me.textLength - 1) / ISA_SAMPLE_RATE + 1) * sizeof(TSize));

    out.write((const char*) me.C, (ValueSize<TValue>::VALUE + 1) * sizeof(TSize));
}

template<typename TValue, typename TSize, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
inline size_t
size(const SadakaneData<TValue, TSize, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> &me)
{
    SEQAN_CHECKPOINT

    size_t compressedPsiSize = 0;
    if(PSI_SAMPLE_RATE > 1) {
        for(TSize i = 0; i < (me.textLength - 1) / PSI_SAMPLE_RATE; ++i) {
            const TBitBin *row = me.compressedPsi[i];

            size_t position = 0;
            for(TSize j = 0; j < PSI_SAMPLE_RATE - 1; ++j) {
                const size_t startPosition = position;
                while(!impl::isBitSet(row, position)) {
                    ++position;
                }

                position += position - startPosition + 1;
            }

            compressedPsiSize += sizeof(TBitBin*) + (size_t) ceil(position / (double) BitsPerValue<TBitBin>::VALUE);
        }

        const TBitBin *row = me.compressedPsi[(me.textLength - 1) / PSI_SAMPLE_RATE];

        size_t position = 0;
        for(TSize j = 0; j < ::std::min((TSize) PSI_SAMPLE_RATE, (TSize) (me.textLength % PSI_SAMPLE_RATE) - 1); ++j) {
            const size_t startPosition = position;
            while(!impl::isBitSet(row, position)) {
                ++position;
            }

            position += position - startPosition + 1;
        }

        compressedPsiSize += sizeof(TBitBin*) + (size_t) ceil(position / (double) BitsPerValue<TBitBin>::VALUE);
    }

    size_t _size = 0;

    _size += sizeof(TSize); // me.textLength

    _size += sizeof(TSize*) + ((me.textLength - 1) / PSI_SAMPLE_RATE + 1) * sizeof(TSize); //me.psiSamples;
    _size += sizeof(TBitBin**) + compressedPsiSize; // me.compressedPsi

    _size += sizeof(TSize*) + ((me.textLength - 1) / SA_SAMPLE_RATE + 1) * sizeof(TSize); // me.I

    _size += sizeof(TSize*) + ((me.textLength - 1) / ISA_SAMPLE_RATE + 1) * sizeof(TSize); // me.J

    _size += sizeof(me.C); // me.C

    return _size;
}

// ============================================================================
// IndexSadakane: Tags, Classes, Enums
// ============================================================================

template<typename TBitBin = unsigned char, size_t PSI_SAMPLE_RATE = 128, size_t SA_SAMPLE_RATE = 16, size_t ISA_SAMPLE_RATE = 16>
struct IndexSadakane {};

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
class Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> > {
protected:
    typedef typename Value<TText>::Type TValue;
    typedef typename SAValue<TText>::Type TSAValue;

public:
    Holder<TText> text;

    SadakaneData<TValue, TSAValue, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> data;

    Index():
    text(),
    data()
    {
        SEQAN_CHECKPOINT
    }

    template<typename TOtherText>
    Index(TOtherText &_text):
    text(_text),
    data()
    {
        SEQAN_CHECKPOINT
    }
};

// ============================================================================
// IndexSadakane: Functions
// ============================================================================

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
inline void
clear(Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> > &me)
{
    SEQAN_CHECKPOINT

    clear(me.text);
    clear(me.data);
}

template <typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
inline bool
indexSupplied(Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> > &me, FibreSA const)
{
    SEQAN_CHECKPOINT

    return me.data.textLength > 0;
}

template <typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
inline bool
indexSolveDependencies(Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> >& me, FibreSA const)
{
    SEQAN_CHECKPOINT

    return indexSupplied(me, FibreText());
}

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
inline bool
indexCreate(Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> > &me, FibreSA const)
{
    SEQAN_CHECKPOINT

    typedef typename Size<TText>::Type TTextSize;

    const TText &text = getFibre(me, FibreText());
    const TTextSize textLength = length(text);

    if(textLength > ULONG_MAX) {
        create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned long long>(me.data, text);
    } else if(textLength > UINT_MAX) {
        create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned long>(me.data, text);
    } else if(textLength > USHRT_MAX) {
        create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned int>(me.data, text);
    } else if(textLength > UCHAR_MAX) {
        create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned short>(me.data, text);
    } else {
        create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned char>(me.data, text);
    }

    return true;
}

// Incremental construction 
template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
inline bool
indexCreate(Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> > &me, FibreSA const, const typename Size<TText>::Type blockLength)
{
    SEQAN_CHECKPOINT
    
    typedef typename Size<TText>::Type TTextSize;
    const TText &text = getFibre(me, FibreText());
    const TTextSize textLength = length(text);
    
    if (blockLength == 0u || blockLength >= textLength) {
        // Fall back to the direct construction (non-incremental)
        return indexCreate(me, FibreSA());
    }
    
    SEQAN_ASSERT_GT_MSG(blockLength, 0u, "");

    if(textLength > ULONG_MAX) {
        if(blockLength >= ULONG_MAX) {
            create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned long long, unsigned long long>(me.data, text, (unsigned long long) blockLength);
        } else if(blockLength >= UINT_MAX) {
            create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned long long, unsigned long>(me.data, text, (unsigned long) blockLength);
        } else if(blockLength >= USHRT_MAX) {
            create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned long long, unsigned int>(me.data, text, (unsigned int) blockLength);
        } else if(blockLength >= UCHAR_MAX) {
            create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned long long, unsigned short>(me.data, text, (unsigned short) blockLength);
        } else {
            create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned long long, unsigned char>(me.data, text, (unsigned char) blockLength);
        }
    } else if(textLength >= UINT_MAX) {
        if(blockLength >= UINT_MAX) {
            create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned long, unsigned long>(me.data, text, (unsigned long) blockLength);
        } else if(blockLength >= USHRT_MAX) {
            create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned long, unsigned int>(me.data, text, (unsigned int) blockLength);
        } else if(blockLength >= UCHAR_MAX) {
            create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned long, unsigned short>(me.data, text, (unsigned short) blockLength);
        } else {
            create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned long, unsigned char>(me.data, text, (unsigned char) blockLength);
        }
    } else if(textLength >= USHRT_MAX) {
        if(blockLength >= USHRT_MAX) {
            create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned int, unsigned int>(me.data, text, (unsigned int) blockLength);
        } else if(blockLength > UCHAR_MAX) {
            create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned int, unsigned short>(me.data, text, (unsigned short) blockLength);
        } else {
            create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned int, unsigned char>(me.data, text, (unsigned char) blockLength);
        }
    } else if(textLength >= UCHAR_MAX) {
        if(blockLength >= UCHAR_MAX) {
            create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned short, unsigned short>(me.data, text, (unsigned short) blockLength);
        } else {
            create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned short, unsigned char>(me.data, text, (unsigned char) blockLength);
        }
    } else {
        create<TText, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE, unsigned char, unsigned char>(me.data, text, (unsigned char) blockLength);
    }

    return true;
}

template<typename TReturnText, typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
inline TReturnText
extract(const Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> > &me, const typename Position<TText>::Type _begin, const typename Position<TText>::Type _end)
{
    SEQAN_CHECKPOINT

    typedef typename Value<TText>::Type TValue;
    typedef typename SAValue<TText>::Type TSAValue;

    const SadakaneData<TValue, TSAValue, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> &data = me.data;

    if(data.textLength < _end || _end <= _begin) {
        return TReturnText();
    }

    const typename Size<TText>::Type textLength = _end - _begin;

    TReturnText text;
    resize(text, textLength);

    TSAValue pos = __isaAt(data, _begin);
    for (typename Position<TText>::Type i = 0; i < textLength; ++i) {
        size_t l = 1, r = ValueSize<TValue>::VALUE;
        while(l <= r) {
            const size_t m = l + (r - l) / 2;

            if(data.C[m - 1] < pos + 1 && pos + 1 <= data.C[m]) {
                text[i] = convert<TValue>(m - 1);

                break;
            } else if(data.C[m] < pos + 1) {
                l = m + 1;
            } else {
                r = m - 1;
            }
        }

        pos = __psiAt(data, pos);
    }

    return text;
}

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
inline typename Value<TText>::Type
extractValue(const Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> > &me, const typename Position<TText>::Type _position)
{
    SEQAN_CHECKPOINT

    typedef typename Value<TText>::Type TValue;
    typedef typename SAValue<TText>::Type TSAValue;

    const SadakaneData<TValue, TSAValue, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> &data = me.data;

    const TSAValue isaPosition = __isaAt(data, _position);

    size_t l = 1, r = ValueSize<TValue>::VALUE;
    while(l <= r) {
        const size_t m = l + (r - l) / 2;

        if(data.C[m - 1] < isaPosition + 1 && isaPosition + 1 <= data.C[m]) {
            return convert<TValue>(m - 1);
        } else if(data.C[m] < isaPosition + 1) {
            l = m + 1;
        } else {
            r = m - 1;
        }
    }

    return TValue();
}

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
inline bool
open(Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> > &me, const char *filename)
{
    SEQAN_CHECKPOINT

    std::ifstream in;

    in.open(filename);
    if(!in.is_open()) {
        return false;
    }

    open(in, me.data);

    in.close();

    return true;
}

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
inline bool
save(const Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> > &me, const char *filename)
{
    SEQAN_CHECKPOINT

    std::ofstream out;

    out.open(filename);
    if(!out.is_open()) {
        return false;
    }

    save(out, me.data);

    out.close();

    return true;
}

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
inline size_t
size(const Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> > &me)
{
    SEQAN_CHECKPOINT

    return sizeof(Holder<TText>) + size(me.data);
}

}  // namespace seqan

#endif // SEQAN_INDEX_INDEX_SADAKANE_H_
