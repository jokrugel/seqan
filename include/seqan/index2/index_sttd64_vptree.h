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
// Author: Alexander Aumann <a.aumann@web.de>
// ==========================================================================
// STTD64 suffix tree data structure and construction algorithm.
// 
// [HST2007] M. Halachev, N. Shiri, A. Thamildurai
//   "Efficient and Scalable Indexing Techniques for Biological Sequence Data",
//   First International Conference on Bioinformatics Research and Development
//   (BIRD'07), Springer, 2007, 4414, 464-479
//   http://dx.doi.org/10.1007/978-3-540-71233-6_36
// ==========================================================================

#ifndef SEQAN_INDEX_INDEX_STTD64_VPTREE_H_
#define SEQAN_INDEX_INDEX_STTD64_VPTREE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

//this tree is used to be able to iterate over the whole tree. It especially also contains
//the last prefixLen suffices (otherwise omitted) and is able to virtually build nodes
//for strings containing less than sequence_length characters
template<typename TextValue>
class _VirtualPrefixTree
{
public:
    //unsigned * prefixesContained;   //stores a flag for each possible prefix, enumerated from left to right
    unsigned alphabetSize;
    unsigned prefixLen;
    std::vector<bool> prefixesContained;
    String<TextValue> prefLenSuffix;    //the last prefLen characters of the input text
    //unsigned repStartSuffix;    //if the iterator is still possible within the suffix
    //unsigned repEndSuffix;      //these contain the start and end positions within the suffix
    //unsigned const sizeOfOneArrayMem;
    _VirtualPrefixTree(unsigned alphabetSize, unsigned prefixLen)
    :
            alphabetSize(alphabetSize), prefixLen(prefixLen),
            prefixesContained(static_cast<unsigned>(pow(static_cast<double>(alphabetSize), static_cast<double>(prefixLen))), false)
            /*, sizeOfOneArrayMem(sizeof prefixesContained[0] * 8)*/
    {
        /* use std::vector<bool> instead, as it does the same space saving*/
        //unsigned arraySize = pow(static_cast<double>(alphabetSize), static_cast<double>(prefixLen));
        //if (arraySize % sizeOfOneArrayMem == 0)
        //    arraySize = arraySize / sizeOfOneArrayMem;
        //else
        //    arraySize = arraySize / sizeOfOneArrayMem + 1;
        //prefixesContained = new unsigned[arraySize];
        //memset(prefixesContained, 0, arraySize * (sizeof prefixesContained[0]));

    }
    ~_VirtualPrefixTree()
    {
        //delete[] prefixesContained;
        //prefixesContained = 0;
    }

    template <typename TSeq>
    inline void setPrefLenSuffix(TSeq const & suf)
    {
        SEQAN_ASSERT_EQ(length(suf), prefixLen - 1);
        //typename Iterator<TSeq const>::Type inputIt; const it comparison prob!#280
        unsigned i;
        resize(prefLenSuffix, prefixLen - 1, Exact());
        for (i = 0; i < length(suf); ++i)
        {
            assignValue(prefLenSuffix, i, getValue(suf, i));
        }

    }

    inline void setPartitionHasValues(unsigned partitionNum)
    {
        //unsigned arrayPos = partitionNum / sizeOfOneArrayMem;
        //unsigned bit = partitionNum % sizeOfOneArrayMem;
        //prefixesContained[arrayPos] |= (1 << (sizeOfOneArrayMem - bit));
        prefixesContained[partitionNum] = true;
    }

    inline unsigned numPartitions() const
    {
        return prefixesContained.size();
    }

    inline bool partitionHasValues(unsigned const partitionNum) const
    {
        //unsigned arrayPos = partitionNum / sizeOfOneArrayMem;
        //unsigned bit = partitionNum % sizeOfOneArrayMem;
        //return ((prefixesContained[arrayPos] & (1 << (sizeOfOneArrayMem - bit))) != 0);
        return prefixesContained[partitionNum];
    }

    template <typename TSeq>
    inline bool children (unsigned & firstPartitionChild,
                          unsigned & lastPartitionChild,
                          bool & hasPartitionChildren,
                              std::vector<unsigned> & prefLenSuffices,
                              TSeq const & seq,
                              unsigned checkInternalSufficesFromPos) const
    {
        SEQAN_CHECKPOINT;
        prefLenSuffices = prefLenSufficesStartingWith(seq, prefLenSuffices,
                                                      checkInternalSufficesFromPos);
        unsigned partChildCount = partitionChildren(firstPartitionChild, lastPartitionChild, seq);
        if (partChildCount >= 1)
        {
            hasPartitionChildren = true;
        }
        else
        {
            hasPartitionChildren = false;
        }
        return (hasPartitionChildren || !prefLenSuffices.empty());
    }

    //determines the right sibling of a prefix
    //goes strictly lexicographical
    template <typename TSeq>
    inline bool
    rightSibling (TSeq & prefix, bool & node, unsigned & partitionNum) const
    {
        SEQAN_CHECKPOINT;
        unsigned rightSuffSiblingCount;
        unsigned leftMostSuffSibling = ordValueOfLeftMostSuffSibling(rightSuffSiblingCount, prefix, false);
        unsigned firstPartSibling, lastPartSibling;
        unsigned partSiblingCount;
        unsigned currentSiblingRightMostOrdVal;
        partitionNum = std::numeric_limits<unsigned>::max();
        TSeq currentSibling = prefix;
        SEQAN_ASSERT_EQ(length(currentSibling), length(prefix));

        for (currentSiblingRightMostOrdVal = ordValue(getValue(prefix, length(prefix) - 1)) + 1;
                currentSiblingRightMostOrdVal < alphabetSize;
                ++currentSiblingRightMostOrdVal)
        {
            assignValue(currentSibling, length(currentSibling) - 1, currentSiblingRightMostOrdVal);
            partSiblingCount = partitionChildren(firstPartSibling, lastPartSibling, currentSibling, 2);
            if (partSiblingCount >= 1 || leftMostSuffSibling == currentSiblingRightMostOrdVal)
            {
                //some sibling has been found

                if (partSiblingCount > 1
                        || (leftMostSuffSibling == currentSiblingRightMostOrdVal && partSiblingCount >= 1))
                {
                    node = true;
                }
                else if (leftMostSuffSibling == currentSiblingRightMostOrdVal/* && partSiblingCount == 0*/)
                {
                    SEQAN_ASSERT_EQ(partSiblingCount, 0u);
                    //check if more than one of the same "suff siblings"
                    std::vector<unsigned> prefSuffixChildrenOfCurSibling;
                    prefSuffixChildrenOfCurSibling = prefLenSufficesStartingWith(currentSibling,
                                                     prefSuffixChildrenOfCurSibling, 0);
                    if (prefSuffixChildrenOfCurSibling.size() >= 2)
                        node = true;
                    else
                        node = false;

                }
                else
                {
                    SEQAN_ASSERT(partSiblingCount == 1 && leftMostSuffSibling != currentSiblingRightMostOrdVal);
                    partitionNum = firstPartSibling;
                    node = false;
                }
                prefix = currentSibling;
                return true;
            }
        }
        return false;
    }

    template <typename TSeq>
    inline unsigned
    ordValueOfLeftMostSuffSibling(unsigned & rightSuffSiblingCount, TSeq const & prefix, bool includeSelfAsSibling) const
    {
        SEQAN_CHECKPOINT;
        std::vector<unsigned> potentialSuffices;
        std::vector<unsigned> sufficesMatching;
        potentialSuffices = prefLenSufficesStartingWith(infix(prefix, 0, length(prefix) - 1),
                                                        sufficesMatching,
                                                        0);
        rightSuffSiblingCount = 0;
        unsigned ordValueOfLeftMostSibling = alphabetSize;
        TSeq curLeftMostRightSibling = prefix;
        SEQAN_ASSERT_EQ(length(curLeftMostRightSibling), length(prefix));
        //potentialSuffices now contains all suffices starting with prefix-without-last-character
        unsigned lastSeqCharOrdVal = ordValue(getValue(prefix, length(prefix) - 1));

        //unsigned selfExclude;
        //if (includeSelfAsSibling)
        //    selfExclude = 0;   //called by leftMostChild
        //else
        //    selfExclude = 1;   //called by rightSibling

        String<TextValue> currentPotSuffix;
        unsigned curPotSuffixOrdValLastChar;
        for (unsigned i = 0; i < potentialSuffices.size(); i++)
        {
            if (potentialSuffices[i] + length(prefix) > prefixLen - 1) continue;
            if (length(prefLenSuffix != 1)) //workaround segment bug
            {
                assign(currentPotSuffix, infix(prefLenSuffix, potentialSuffices[i],
                                               potentialSuffices[i] + length(prefix)));
            }
            else
            {
                //workaround
                resize(currentPotSuffix, 1);
                assignValue(currentPotSuffix, 0, getValue(prefLenSuffix, potentialSuffices[i]));
            }
            SEQAN_ASSERT_EQ(length(currentPotSuffix), length(prefix));
            curPotSuffixOrdValLastChar
                = ordValue(getValue(currentPotSuffix, length(currentPotSuffix) - 1));
            if (curPotSuffixOrdValLastChar >= lastSeqCharOrdVal && includeSelfAsSibling)
            {
                ++rightSuffSiblingCount;
                if (curPotSuffixOrdValLastChar < ordValueOfLeftMostSibling)
                {
                    ordValueOfLeftMostSibling = curPotSuffixOrdValLastChar;
                }

            }
            else if (curPotSuffixOrdValLastChar > lastSeqCharOrdVal)
            {
                ++rightSuffSiblingCount;
                if (curPotSuffixOrdValLastChar < ordValueOfLeftMostSibling)
                {
                    ordValueOfLeftMostSibling = curPotSuffixOrdValLastChar;
                }

            }

        }

        //for (unsigned i = lastSeqCharOrdVal + selfExclude; i < alphabetSize; i++)
        //{
        //    assignValue(curLeftMostRightSibling, length(curLeftMostRightSibling) - 1, i);
        //    sufficesMatching = prefLenSufficesStartingWith(curLeftMostRightSibling,
        //                                                potentialSuffices,
        //                                                length(prefix) - 1);
        //    if (sufficesMatching.size() >= 1)
        //    {
        //        TSeq matchingSuffix;
        //        getSuffixAtRelativePos(matchingSuffix, sufficesMatching[0], length(prefix));
        //        SEQAN_ASSERT_EQ(length(prefix), length(matchingSuffix));
        //        if (ordValueOfLeftMostSibling > ordValue(getValue(matchingSuffix, length(prefix) - 1)))
        //        {
        //            ordValueOfLeftMostSibling = ordValue(getValue(matchingSuffix, length(prefix) - 1));
        //        }
        //        rightSuffSiblingCount += sufficesMatching.size();
        //    }
        //}
        return ordValueOfLeftMostSibling;
    }

    template <typename TSeq>
    inline bool leftMostChild(TSeq & leftMostChild, bool & node, unsigned & partChildNum,
                              TSeq const & seq) const
    {
        SEQAN_CHECKPOINT;
        //TODO(aumann): assume no leaf, optimize?
        partChildNum = std::numeric_limits<unsigned>::max();
        unsigned firstChild, lastChild;
        unsigned partChildrenCountFather = partitionChildren(firstChild, lastChild, seq, 2);
        SEQAN_ASSERT(partChildrenCountFather <= 2);
        TSeq leftmostPossibleChildSequence = seq;
        appendValue(leftmostPossibleChildSequence, 0);
        unsigned suffChildrenFatherCount, ordValueLeftMostSuffChild;
        ordValueLeftMostSuffChild = ordValueOfLeftMostSuffSibling(suffChildrenFatherCount,
                                                                  leftmostPossibleChildSequence,
                                                                  true);
        //std::vector<unsigned> prefSufficesStartingWithSeq
        //        = prefLenSufficesStartingWith(seq, std::vector<unsigned>(), 0);
        TSeq leftmostPrefSuffix;
        TSeq curPrefSuffix;
        bool prefSuffixIsLeftMost = false;
        if (suffChildrenFatherCount > 0)
        {

            //replace last character of leftMostPossibleChild by the leftMost ordValue
            //resulting in the leftmost prefSuffix matching
            leftmostPrefSuffix = leftmostPossibleChildSequence;
            assignValue(leftmostPrefSuffix, length(leftmostPrefSuffix) - 1, ordValueLeftMostSuffChild);
            SEQAN_ASSERT_EQ(length(leftmostPrefSuffix), length(seq) + 1);
        }
        
        //if (prefSufficesStartingWithSeq.size() > 0)
        //{
        //    getSuffixAtRelativePos(leftmostPrefSuffix, prefSufficesStartingWithSeq[0], length(seq));
        //    //find the leftmost of them
        //    for (unsigned i = 1; i < prefSufficesStartingWithSeq.size(); i++)
        //    {
        //        getSuffixAtRelativePos(curPrefSuffix, prefSufficesStartingWithSeq[i], length(seq));
        //        if (curPrefSuffix < leftmostPrefSuffix)
        //        {
        //            leftmostPrefSuffix = curPrefSuffix;
        //        }
        //    }
        //}

        if (partChildrenCountFather >= 1 && suffChildrenFatherCount == 0)
        {
            _getPrefixForPartitionNumber(leftMostChild, firstChild, alphabetSize, prefixLen, length(seq) + 1);

            //count partition children to check if "virtual leaf" (directly go down to partition) or not
            //signal that the iterator now points to a node: node = false, partChildNum != limits::max()
            unsigned partChildrenChild;
            if (length(leftMostChild) < prefixLen - 1)
            {
                unsigned dumpFirst, dumpLast;
                partChildrenChild = partitionChildren(dumpFirst, dumpLast,
                                                               leftMostChild, 2);
                SEQAN_ASSERT_EQ(dumpFirst, firstChild);    //if not use dumpFirst below (if at all sensible that assertion doesn't hold)
            }
            else
            {
                partChildrenChild = 1;
            }
            if (partChildrenChild >= 2)
            {
                node = true;
            }
            else
            {
                partChildNum = firstChild;  //signal that iterator points at a partition now
                node = false;
            }
            return true;
        }
        if (partChildrenCountFather == 0 && suffChildrenFatherCount >= 1)   //else if
        {
            leftMostChild = leftmostPrefSuffix;
            //count how many prefSuffices start with "leftmostPrefSuffix" in order to find
            //out if node or leaf
            std::vector<unsigned> prefSufficesStartingWithLM
                    = prefLenSufficesStartingWith(leftmostPrefSuffix, std::vector<unsigned>(), 0);
            SEQAN_ASSERT_NOT(prefSufficesStartingWithLM.empty());
            if (prefSufficesStartingWithLM.size() > 1)
            {
                node = true;
            }
            else
            {
                node = false;
            }
            return true;
        }

        if (partChildrenCountFather >= 1 && suffChildrenFatherCount >= 1)
        {
            //determine which of these is the leftmost
            _getPrefixForPartitionNumber(leftMostChild, firstChild, alphabetSize, prefixLen, length(seq) + 1);
            if (leftmostPrefSuffix < leftMostChild)
            {
                leftMostChild = leftmostPrefSuffix;
                prefSuffixIsLeftMost = true;
            }
            else
            {
                prefSuffixIsLeftMost = false;
            }
            //now count how many suff and partition children the leftMostChild has
            std::vector<unsigned> prefSufficesStartingWithLM
                    = prefLenSufficesStartingWith(leftMostChild, std::vector<unsigned>(), 0);
            unsigned firstPartChildLM, lastPartChildLM, partChildrenLMChildCount;
            if (length(leftMostChild) < prefixLen)
            {
                partChildrenLMChildCount = partitionChildren(firstPartChildLM, lastPartChildLM,
                                                                      leftMostChild, 2);
                SEQAN_ASSERT_EQ(firstPartChildLM, firstChild);
            }
            else
            {
                partChildrenLMChildCount = 1;
            }
            SEQAN_ASSERT(prefSufficesStartingWithLM.size() > 0 || partChildrenLMChildCount > 0);
            SEQAN_ASSERT(prefSuffixIsLeftMost || partChildrenLMChildCount >= 1);
            if (!prefSuffixIsLeftMost)
            {
                if (prefSufficesStartingWithLM.empty() && partChildrenLMChildCount == 1)
                {
                    //found a partition, signal that we need to go down to it
                    node = false;
                    partChildNum = firstChild;
                }
                else if (!prefSufficesStartingWithLM.empty() || partChildrenLMChildCount >= 2)
                {
                    node = true;
                }
                else
                {
                    //only one partition child
                    node = false;
                }
            }
            else
            {
                SEQAN_ASSERT_GEQ(prefSufficesStartingWithLM.size(), 1u);
                if (partChildrenLMChildCount >= 1 || prefSufficesStartingWithLM.size() >= 2)
                {
                    node = true;
                }
                else
                {
                    node = false;
                }
            }
            return true;
        }   //partChildrenCountFather >= 1 && suffChildrenFatherCount >= 1

        return false;

//            if (partChildrenCountFather >= 1)                                   //else if
//            {
//                _getPrefixForPartitionNumber(leftMostChild, firstChild, alphabetSize, prefixLen, length(seq) + 1);
//                SEQAN_ASSERT_EQ(length(leftMostChild), length(seq) + 1);
//                SEQAN_ASSERT_EQ(length(leftmostPrefSuffix), length(leftMostChild));
//                if (suffChildrenFatherCount > 0)
//                {
//                    if (leftmostPrefSuffix == leftMostChild)
//                    {
//                        node = true;
//                    }
//                    if (leftmostPrefSuffix < leftMostChild)
//                    {
//                        leftMostChild = leftmostPrefSuffix;
//                    }
//                }
//                else
//                {
//                    //check whether the partition child has a unique prefix (=> "virtual leaf")
//                    unsigned dumpFirst, dumpLast;
//                    unsigned partChildrenChild = partitionChildren(dumpFirst, dumpLast,
//                                                                   leftMostChild, 2);
//                    SEQAN_ASSERT_EQ(dumpFirst, firstChild);    //if not use dumpFirst below (if at all sensible that assertion doesn't hold)
//                    if (partChildrenChild >= 2)
//                    {
//                        node = true;
//                    }
//                    else
//                    {
//                        partChildNum = firstChild;  //signal that iterator points at a partition now
//                        node = false;
//                    }
//
//                }
//                return true;
//            }
//            else
//            {
//                if (suffChildrenFatherCount > 0)
//                {
//                    leftMostChild = leftmostPrefSuffix;
//                    if (suffChildrenFatherCount >= 2)
//                    {
//                        node = true;
//                    }
//                    else
//                    {
//                        node = false;
//                    }
//                    return true;
//                }
//                return false;
//            }
    }

    template <typename TSeq>
    inline void
    getSuffixAtRelativePos(TSeq & seq, unsigned relStartPos,
                                       typename Size<TSeq>::Type length) const
    {
        SEQAN_CHECKPOINT;
        resize(seq, length);
        SEQAN_ASSERT_LEQ(relStartPos + length, prefixLen - 1);
        for (unsigned i = 0; i < length; i++) {
            assignValue(seq, i, getValue(prefLenSuffix, relStartPos + i));
        }

    }

    template <typename TSeq>
    inline std::vector<unsigned>
    prefLenSufficesStartingWith (TSeq const & seq,
                                 std::vector<unsigned> const & potentialSuffices,
                                 unsigned checkInternalSufficesFromPos) const
    {
        SEQAN_CHECKPOINT;
        std::vector<unsigned> potentialSufficesTMP;
        std::vector<unsigned> sufficesStartingWithSeq;
        SEQAN_ASSERT_EQ(length(prefLenSuffix), prefixLen - 1);
        if (checkInternalSufficesFromPos == 0)
        {
            //all preflen tree suffices are potential matches
            //fill list with every starting position
            potentialSufficesTMP.reserve(length(prefLenSuffix));
            for (unsigned i = 0; i < length(prefLenSuffix); i++)
            {
                potentialSufficesTMP.push_back(i);
            }
        } else {
            potentialSufficesTMP.assign(potentialSuffices.begin(), potentialSuffices.end());
        }

        if (length(seq) == 0) return potentialSufficesTMP;  //all children

        //now check all prefLenSuffices for a match
        //simple quadratic approach, as the length is <= 7 this sould not matter
        std::vector<unsigned>::const_iterator potSufficesIt;
        typename Iterator<const TSeq>::Type seqIt;
        typename Iterator<const String<TextValue> >::Type prefLenSuffIt;

        bool misMatch;
        for (potSufficesIt = potentialSufficesTMP.begin();
                potSufficesIt < potentialSufficesTMP.end();
                potSufficesIt++)
        {
            if (length(seq) + 1 + *potSufficesIt > prefixLen) continue; //match impossible

            seqIt = iter(seq, checkInternalSufficesFromPos);
            prefLenSuffIt = iter(prefLenSuffix, *potSufficesIt + checkInternalSufficesFromPos);
            misMatch = false;
            for (; seqIt != end(seq); seqIt++, prefLenSuffIt++)
            {
                if (*seqIt != *prefLenSuffIt)
                {
                    misMatch = true;
                    break;
                }
            }
            if (!misMatch)
            {
                sufficesStartingWithSeq.push_back(*potSufficesIt);
            }
        }
        return sufficesStartingWithSeq;
    }

    template <typename TSeq>
    inline unsigned
    partitionChildren(unsigned & firstChild, unsigned & lastChild, TSeq const & sequence,
                      unsigned searchLimit) const
    {
        SEQAN_CHECKPOINT;
        unsigned startPartition;
        unsigned endPartition;
        unsigned count = 0;
        lastChild = std::numeric_limits<unsigned>::max();
        SEQAN_ASSERT_LEQ(length(sequence), prefixLen);
        SEQAN_ASSERT_GEQ(searchLimit, 1u);
        getStartAndEndOfPartition(startPartition, endPartition, sequence);

        if (startPartition == endPartition)
        {
            if (prefixesContained[startPartition] == true)
            {
                firstChild = lastChild = startPartition;
                return 1;
            }
            firstChild = lastChild = std::numeric_limits<unsigned>::max();
            return 0;
        }
        bool childFound = false;
        std::vector<bool>::const_iterator prefConIt;
        unsigned pos;
        for (prefConIt = prefixesContained.begin() + startPartition, pos = startPartition;
                prefConIt <= prefixesContained.begin() + endPartition && count < searchLimit;
                ++prefConIt, ++pos)
        {
            if (*prefConIt == true)
            {
                ++count;
                if (!childFound) {
                    childFound = true;
                    firstChild = pos;
                }
                else
                {
                    lastChild = pos;
                }
            }
        }
        return count;
    }

    //determines the start and end partition numbers for the current prefix prefix
    //e.g. prefixLen = 3 and sequence = "AG" for DNA alphabet
    //result would be: start= partition number "AGA"; end = partitionNumber("AGT")
    template <typename TSeq>
    inline void
    getStartAndEndOfPartition(unsigned & start, unsigned & end, TSeq const & sequence) const
    {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_LEQ(length(sequence), prefixLen);
        typedef String<typename Value<TSeq>::Type > TStr;
        if (length(sequence) < prefixLen) {
            TStr beginStr;
            TStr endStr;
            typename Value<TSeq>::Type firstOfAlphabet = 0;
            typename Value<TSeq>::Type lastOfAlphabet = alphabetSize - 1;
            beginStr = sequence;
            endStr = sequence;
            appendValue(beginStr, firstOfAlphabet);
            appendValue(endStr, lastOfAlphabet);
            _getPartitionNumber(start, beginStr, alphabetSize, prefixLen);
            _getPartitionNumber(end, endStr, alphabetSize, prefixLen);
        }
        else
        {
            _getPartitionNumber(start, sequence, alphabetSize, prefixLen);
            end = start;
        }
    }

}; // VirtualPrefixTree

} //namespace seqan

#endif // SEQAN_INDEX_INDEX_STTD64_VPTREE_H_
