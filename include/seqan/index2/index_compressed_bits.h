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
// Bit-packed string.
// ==========================================================================

#ifndef SEQAN_INDEX_INDEX_COMPRESSED_BITS_H_
#define SEQAN_INDEX_INDEX_COMPRESSED_BITS_H_

#include <cmath>

namespace seqan {

namespace impl
{
    template<typename TSize>
    inline unsigned char
    bitsPerValue(TSize value)
    {
        SEQAN_CHECKPOINT
        
        return (unsigned char) (log((double) value) / log(2.0)) + 1;
    }
    
    template<typename TValue>
    inline bool
    isBitSet(const TValue &me, const unsigned char position)
    {
        SEQAN_CHECKPOINT
        
        return me & (TValue) ((TValue) 1 << position);
    }
    
    template<typename TValue>
    inline void
    setBit(TValue &me, const unsigned char position)
    {
        SEQAN_CHECKPOINT
        
        me |= (TValue) ((TValue) 1 << position);
    }
    
    template<typename TValue>
    inline void
    clearBit(TValue &me, const unsigned char position)
    {
        SEQAN_CHECKPOINT
        
        me &= ~((TValue) 1 << position);
    }
    
    template<typename TBitBin, typename TSize>
    inline bool
    isBitSet(const TBitBin *me, const TSize position)
    {
        SEQAN_CHECKPOINT
        
        return isBitSet(me[position / BitsPerValue<TBitBin>::VALUE], position % BitsPerValue<TBitBin>::VALUE);
    }
    
    template<typename TBitBin, typename TSize>
    inline void
    clearBit(TBitBin *me, const TSize position)
    {
        SEQAN_CHECKPOINT
        
        clearBit(me[position / BitsPerValue<TBitBin>::VALUE], position % BitsPerValue<TBitBin>::VALUE);
    }
    
    template<typename TBitBin, typename TSize>
    inline void
    setBit(TBitBin *me, const TSize position)
    {
        SEQAN_CHECKPOINT
        
        setBit(me[position / BitsPerValue<TBitBin>::VALUE], position % BitsPerValue<TBitBin>::VALUE);
    }
}

}  // namespace seqan

#endif // SEQAN_INDEX_INDEX_COMPRESSED_BITS_H_
