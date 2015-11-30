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
// Pattern for IndexDigest.
// The pattern object is mainly responsible for converting the needle into a 
// binary representation and to iterate over the binary string.
// ==========================================================================

#ifndef SEQAN_INDEX_FIND_INDEX_DIGEST_PATTERN_H_
#define SEQAN_INDEX_FIND_INDEX_DIGEST_PATTERN_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// DigestBitIterator
//
// Iterate over the bits of a string.
// Therefore we iterate over the characters of the string (using currentPos)
// and over the bits of the current character (using bitPos and mask).
// ----------------------------------------------------------------------------

template <typename TString>
struct DigestBitIterator {
    typedef typename Value<TString>::Type                   TValue;
    typedef typename MakeUnsigned_<TValue>::Type            TMask;      // Yields e.g. unsigned char or SimpleType<unsigned char, Dna_>.
    typedef typename Value<TMask>::Type                     TMaskValue; // Yields in both above cases unsigned char.
    typedef typename Size<TString>::Type                    TSize;
    	
    static const unsigned BITS_PER_VALUE = BitsPerValue<TValue>::VALUE;
    static const unsigned MASK_HIGHEST_BIT = 1u << (BITS_PER_VALUE - 1u);
    
    Holder<TString> string;
    TSize stringLength;
    TSize currentPos;
    unsigned bitPos;
    TMaskValue mask;

    DigestBitIterator() {}
    
    DigestBitIterator(TString & string2):
        string(string2),
        stringLength(length(string2)),
        currentPos(0u),
        bitPos(BITS_PER_VALUE),
        mask(MASK_HIGHEST_BIT) { }

    bool operator *() {
        SEQAN_CHECKPOINT;
        return (mask & (TMaskValue) value(string)[currentPos]) != 0u;
    }

    void operator ++() {
        SEQAN_CHECKPOINT;
        mask >>= 1u;
        --bitPos;
        if (mask == 0u) {
            bitPos = BITS_PER_VALUE;
            mask = MASK_HIGHEST_BIT;
            ++currentPos;
        }
    }

    void operator +=(TSize i) {
        SEQAN_CHECKPOINT;
        if (i < bitPos) {
			bitPos -= i;
            mask >>= i;
            return;
        }
        i -= bitPos;
        
        TSize diffIndex = (i / BITS_PER_VALUE);
        currentPos += diffIndex + 1u;
        i -= diffIndex * BITS_PER_VALUE;
        if (currentPos < stringLength) {
            mask = 1u << (BITS_PER_VALUE - 1u - i);
            bitPos = BITS_PER_VALUE - i;
        }
    }
};

// ----------------------------------------------------------------------------
// Pattern
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec = Standard>
class PatternDigest {};

template <typename TNeedle, typename TIndex, typename TSpec>
struct Pattern<TNeedle, PatternDigest<TIndex, TSpec> > {
    typedef typename TIndex::TDigestSuffix                  TDigestSuffix;
    typedef DigestBitIterator<TNeedle>                      TBitIterator;

    Holder<TNeedle> _host;      // The needle we work on.
    TDigestSuffix digestSuffix; // The preprocessed pattern in the same format as the dividers.
    TBitIterator itBit;         // To iterate over the individual bits of the pattern.

    Pattern() {}

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle)
        //: _host(needle)
    {
        SEQAN_CHECKPOINT;
        setHost(*this, needle);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// We need to preprocess the pattern as well, depending on the config of the index.
template <typename TNeedle, typename THaystack, typename TSpec>
struct DefaultIndexPattern<TNeedle, Index<THaystack, IndexDigest<TSpec> > > {
    typedef Index<THaystack, IndexDigest<TSpec> >           TIndex;    
    typedef PatternDigest<TIndex>                           Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// atEnd
// ----------------------------------------------------------------------------

template <typename TString>
inline bool
atEnd(DigestBitIterator<TString> const & itBit) {
    return itBit.currentPos >= itBit.stringLength;
}

template <typename TString>
inline bool
atEnd(DigestBitIterator<TString> & itBit) {
    return itBit.currentPos >= itBit.stringLength;
}

// ----------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------

template <typename TString>
inline void
clear(DigestBitIterator<TString> & itBit) {
    typedef DigestBitIterator<TString>                      TBitIterator;
    
    SEQAN_CHECKPOINT;
    itBit.currentPos = 0u,
    itBit.bitPos = TBitIterator::BITS_PER_VALUE;
    itBit.mask = TBitIterator::MASK_HIGHEST_BIT;
}

template <typename TNeedle, typename TIndex, typename TSpec>
inline void
clear(Pattern<TNeedle, PatternDigest<TIndex, TSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    clear(pattern.itBit);
}

// ----------------------------------------------------------------------------
// host, needle
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TIndex, typename TSpec>
TNeedle const &
host(Pattern<TNeedle, PatternDigest<TIndex, TSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}

template <typename TNeedle, typename TIndex, typename TSpec>
TNeedle &
host(Pattern<TNeedle, PatternDigest<TIndex, TSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}

template <typename TNeedle, typename TIndex, typename TSpec>
TNeedle const &
needle(Pattern<TNeedle, PatternDigest<TIndex, TSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    return host(pattern);
}

template <typename TNeedle, typename TIndex, typename TSpec>
TNeedle &
needle(Pattern<TNeedle, PatternDigest<TIndex, TSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    return host(pattern);
}

// ----------------------------------------------------------------------------
// length
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TIndex, typename TSpec>
typename Size<TNeedle>::Type
length(Pattern<TNeedle, PatternDigest<TIndex, TSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    return length(needle(pattern));
}

// ----------------------------------------------------------------------------
// setHost
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TIndex, typename TSpec, typename TNeedle2>
inline void
setHost(Pattern<TNeedle, PatternDigest<TIndex, TSpec> > & pattern, TNeedle2 & needle2) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, PatternDigest<TIndex, TSpec> > TPattern;
    typedef typename TIndex::TPrefix                        TPrefix;
    typedef typename TPattern::TDigestSuffix                TDigestSuffix;
    typedef typename TPattern::TBitIterator                 TBitIterator;
    
    clear(pattern);
    setValue(pattern._host, needle2);
    TPrefix patternPrefix;
    _getPrefix<TIndex::PREFIX_LENGTH>(patternPrefix, needle(pattern), 0u);
    pattern.digestSuffix = TDigestSuffix(0u, patternPrefix);
    // TODO(krugel) Implement setHost for BitIterator
    pattern.itBit = TBitIterator(needle(pattern));
}

//template <typename TNeedle, typename TIndex, typename TSpec, typename TNeedle2>
//inline void
//setHost(Pattern<TNeedle, PatternDigest<TIndex, TSpec> > & pattern, TNeedle2 & needle2) {
//    setHost(pattern, const_cast<TNeedle2 const &>(needle2));
//}

//template <typename TNeedle, typename TIndex, typename TSpec, typename TNeedle2>
//inline void
//setNeedle(Pattern<TNeedle, PatternDigest<TIndex, TSpec> > & pattern, TNeedle2 const & needle2) {
//    setHost(pattern, needle2);
//}

template <typename TNeedle, typename TIndex, typename TSpec, typename TNeedle2>
inline void
setNeedle(Pattern<TNeedle, PatternDigest<TIndex, TSpec> > & pattern, TNeedle2 & needle2) {
    //setHost(pattern, const_cast<TNeedle2 const &>(needle2));
    setHost(pattern, needle2);
}

// ----------------------------------------------------------------------------
// begin, end, beginPosition, endPosition
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TIndex, typename TSpec, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type
begin(Pattern<TNeedle, PatternDigest<TIndex, TSpec> > const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    return begin(needle(pattern), spec);
}

template <typename TNeedle, typename TIndex, typename TSpec, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type
end(Pattern<TNeedle, PatternDigest<TIndex, TSpec> > const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    return end(needle(pattern), spec);
}

template <typename TNeedle, typename TIndex, typename TSpec>
typename Position<TNeedle>::Type
beginPosition(Pattern<TNeedle, PatternDigest<TIndex, TSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    return 0u;
}

template <typename TNeedle, typename TIndex, typename TSpec>
typename Position<TNeedle>::Type
endPosition(Pattern<TNeedle, PatternDigest<TIndex, TSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    return length(needle(pattern));
}

}  // namespace seqan

#endif  // SEQAN_INDEX_FIND_INDEX_DIGEST_PATTERN_H_
