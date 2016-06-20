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
// Author: Tobias Stadler <tobiasstadler@mytum.de>
// ==========================================================================
// Finder class for the compressed suffix array of Sadakane.
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

#ifndef SEQAN_INDEX_FIND_INDEX_SADAKANE_H_
#define SEQAN_INDEX_FIND_INDEX_SADAKANE_H_

namespace seqan {

template<typename TSpec = void>
struct FinderSadakane {};

// ============================================================================
// Metafunctions
// ============================================================================

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE>
struct DefaultFinder<Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> > > {
    typedef FinderSadakane<> Type;
};

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE, typename TSpecFinder >
class Finder< Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> >, FinderSadakane<TSpecFinder>  >
{
protected:
    typedef Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> >    TIndex;
    typedef typename Value<TIndex>::Type                                                               TValue;
    typedef typename SAValue<TIndex>::Type                                                             TSAValue;
    typedef typename Size<TIndex>::Type                                                                TSize;

public:
    Holder<TIndex> index;
    Pair<TSAValue> range;
    TSAValue data_iterator;
    TSAValue data_length;

    Finder()
    {
        clear(*this);
    }

    Finder(TIndex &_index):
    index()
    {
        clear(*this);

        index = _index;
    }
};

// ============================================================================
// Functions
// ============================================================================

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE, typename TSpecFinder>
inline void
clear(Finder< Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> >, FinderSadakane<TSpecFinder>  > &me)
{
    me.range.i1 = me.range.i2 = me.data_iterator = 0;
    me.data_length = 0;
}

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE, typename TSpecFinder>
inline typename Position< Finder< Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> >, FinderSadakane<TSpecFinder>  > >::Type &
hostIterator(Finder< Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> >, FinderSadakane<TSpecFinder> > &me)
{
    return me.data_iterator;
}

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE, typename TSpecFinder>
inline typename Position< Finder< Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> >, FinderSadakane<TSpecFinder>  > >::Type const &
hostIterator(Finder< Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> >, FinderSadakane<TSpecFinder> > const &me)
{
    return me.data_iterator;
}

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE, typename TSpecFinder>
inline typename Position< Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> > >::Type
beginPosition(Finder< Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> >, FinderSadakane<TSpecFinder> > &me)
{
    SEQAN_ASSERT(!empty(me));

    return __saAt(host(me).data, me.data_iterator);
}

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE, typename TSpecFinder>
inline typename Position<Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> > >::Type
beginPosition(Finder< Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> >, FinderSadakane<TSpecFinder> > const &me)
{
    SEQAN_ASSERT(!empty(me));

    return __saAt(host(me).data, me.data_iterator);
}

template<typename TText, typename TBitBin, size_t PSI_SAMPLE_RATE, size_t SA_SAMPLE_RATE, size_t ISA_SAMPLE_RATE, typename TSpecFinder, typename TPatternText, typename TFinderSpec2>
inline void _findFirstIndex(Finder<Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> >, FinderSadakane<TSpecFinder> > &me, const TPatternText &pattern, const TFinderSpec2)
{
    typedef typename Value<TText>::Type              TValue;
    typedef typename SAValue<TText>::Type            TSAValue;

    Index<TText, IndexSadakane<TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> > &index = haystack(me);
    const SadakaneData<TValue, TSAValue, TBitBin, PSI_SAMPLE_RATE, SA_SAMPLE_RATE, ISA_SAMPLE_RATE> &data = index.data;

    indexRequire(index, FibreSA());

    // if the pattern has zero length or is greater the the text length, there is no occurence of pattern in the text
    typename Size<TPatternText>::Type patternLength = length(pattern);
    if(patternLength == 0 || data.textLength < patternLength) {
        me.range = Pair<TSAValue, TSAValue>(0, 0);

        return;
    }

    // use backward search to find the pattern
    TSAValue _begin = 0, _end = data.textLength;
    for (typename Position<TPatternText>::Type i = patternLength; i > 0; --i) {
        const size_t character = ordValue(convert<TValue>(pattern[i - 1]));

        const TSAValue prevBegin = _begin, prevEnd = _end;

        _begin = data.textLength + 1;
        _end = 0;

        TSAValue l = data.C[character], r = data.C[character + 1] - 1;
        while (l <= r) {
            const TSAValue m = l + (r - l) / 2;

            if(__psiAt(data, m) < prevBegin) {
                l = m + 1;
            } else {
                r = m - 1;

                if(__psiAt(data, m) <= prevEnd) {
                    _begin = m;
                }
            }
        }

        l = _begin, r = data.C[character + 1] - 1;
        while (l <= r) {
            const TSAValue m = l + (r - l) / 2;

            if (__psiAt(data, m) <= prevEnd) {
                l = m + 1;

                _end = m;
            } else {
                r = m - 1;
            }
        }

        // now, all occurences of suffix(pattern, i - 1) are SA[_begin .. end]

        if (_begin > _end) {
            me.range = Pair<TSAValue, TSAValue>(0, 0);

            return;
        }
    }

    me.range = Pair<TSAValue, TSAValue>(_begin, _end + 1);
}

}  // namespace seqan

#endif // SEQAN_INDEX_FIND_INDEX_SADAKANE_H_
