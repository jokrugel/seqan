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

#ifndef SEQAN_INDEX_INDEX_STTD64_BASE_H_
#define SEQAN_INDEX_INDEX_STTD64_BASE_H_

#ifndef CONF_SW_EXT_MMAP
// TODO(krugel) External failed on Windows with MinGW and e.g. ./data/prefix/dna-human4-14.txt
#define CONF_SW_EXT_MMAP External
//#define CONF_SW_EXT_MMAP MMap
#endif

#include <limits>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template<typename TextValue>
class _VirtualPrefixTree;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

//STTD64 fibres
struct Sttd64SuffixTrees_;
struct Sttd64Partitions_;
struct Sttd64VirtualTree_;
//struct Sttd64Suffices_;

typedef FibreText Sttd64Text;
typedef FibreRawText Sttd64RawText;
typedef Tag<Sttd64SuffixTrees_> const Sttd64SuffixTrees;
typedef Tag<Sttd64Partitions_> const Sttd64Partitions;
typedef Tag<Sttd64VirtualTree_> const Sttd64VirtualTree;
//typedef Tag<Sttd64Suffices_> const Sttd64Suffices;

struct Sttd64Node_
{
    unsigned lp;
    unsigned pointerOrDepth;

    Sttd64Node_()
    : lp(0), pointerOrDepth(0)
    {
    }
    inline void setLeafBit(bool set) {
        if (set) {
            pointerOrDepth |= (1 << 31);
        } else {
            pointerOrDepth &= 0x7fffffff;
        }
    }
    inline void setRightMostBit(bool set) {
        if (set) {
            pointerOrDepth |= (1 << 30);
        } else {
            pointerOrDepth &= 0xbfffffff;
        }
    }
    inline unsigned cleanPointerOrDepth() const {
        return (pointerOrDepth & 0x3fffffff);
    }
    inline bool isLeaf() const {
        if ((pointerOrDepth & (1 << 31)) == 0) {
            return false;
        }
        return true;
    }
    inline bool isRightMost() const {
        if ((pointerOrDepth & (1 << 30)) == 0) {
            return false;
        }
        return true;
    }

    CharString toString() {
        CharString retVal;
        std::stringstream cStream;
        cStream << "[";
        if (isLeaf()) {
            cStream << "L";
        } else {
            cStream << "N";
        }

        cStream << "(" << lp << "," << cleanPointerOrDepth() << ")";

        if (isRightMost()) {
            cStream << "R";
        }
        cStream << "]";
        retVal = cStream.str();
        return retVal;
    }
};

template<unsigned PREFIX_LENGTH_ = 1,
        unsigned PAGESIZE_ = 8192,
        unsigned STRBUFF_FRAME_NUM_ = 5, //250 MB
        unsigned TREEBUFF_FRAME_NUM_ = 5,
        unsigned TEMPBUFF_FRAME_NUM_ = 5, ///1GB
        unsigned SUFBUFF_FRAME_NUM_ = 5, ///1GB
        unsigned TOTAL_MEM_SIZE_ = 2u * 1024 * 1024 * 1024,
        bool DYNAMIC_BUFFER_SIZE_CFG_ = true>
struct Sttd64Config
{
    //typedef unsigned Type;
    //typedef TPrefixType_ TPrefixType;

    static unsigned const PREFIX_LENGTH = PREFIX_LENGTH_;
    static unsigned const PAGESIZE = PAGESIZE_;
    static unsigned const STRBUFF_FRAME_NUM = STRBUFF_FRAME_NUM_;
    static unsigned const TREEBUFF_FRAME_NUM = TREEBUFF_FRAME_NUM_;
    static unsigned const TEMPBUFF_FRAME_NUM = TEMPBUFF_FRAME_NUM_;
    static unsigned const SUFBUFF_FRAME_NUM = SUFBUFF_FRAME_NUM_;
    static unsigned const TOTAL_MEM_SIZE = TOTAL_MEM_SIZE_;

    static bool const DYNAMIC_BUFFER_SIZE_CFG = DYNAMIC_BUFFER_SIZE_CFG_;
    static unsigned const SIZE_T_TEMP = 4;
    static unsigned const SIZE_T_SUF = 4;
    static unsigned const SIZE_T_TREE = 8;
    //optimum would be a dynamically calculated buffer size which depends
    //on the given alphabet size,
    //as described in Tian et al.'s work about TDD
    //thus the additional parameters dynamic_buffer_size_cfg and TOTAL_MEM_SIZE

    template<typename TFile = ExternalConfig<>::TFile, unsigned PAGESIZE = PAGESIZE, unsigned FRAMES = TREEBUFF_FRAME_NUM>
    struct TreeBuffConfig
    {
        typedef ExternalConfig<TFile, PAGESIZE / SIZE_T_TREE, FRAMES> Type;
    };

    template<typename TFile = ExternalConfig<>::TFile, unsigned PAGESIZE = PAGESIZE, unsigned FRAMES = SUFBUFF_FRAME_NUM>
    struct SufBuffConfig
    {
        typedef ExternalConfig<TFile, PAGESIZE / SIZE_T_SUF, FRAMES> Type;
    };

    template<typename TFile = ExternalConfig<>::TFile, unsigned PAGESIZE = PAGESIZE, unsigned FRAMES = 2>
    struct SufBuffPartitioningConfig
    {
        typedef ExternalConfig<TFile, PAGESIZE, FRAMES> Type;
    };

    template<typename TFile = ExternalConfig<>::TFile, unsigned PAGESIZE = PAGESIZE, unsigned FRAMES = STRBUFF_FRAME_NUM>
    struct StrBuffConfig
    {
        typedef ExternalConfig<TFile, PAGESIZE, FRAMES> Type;
    };

    template<typename TFile = ExternalConfig<>::TFile, unsigned PAGESIZE = PAGESIZE, unsigned FRAMES = TEMPBUFF_FRAME_NUM>
    struct TempBuffConfig
    {
        typedef ExternalConfig<TFile, PAGESIZE / SIZE_T_TEMP, FRAMES> Type;
    };
};

template<typename TConfig = Sttd64Config<> >
struct IndexSttd64 {};

template<typename _TText, typename TConfig>
class Index<_TText, IndexSttd64<TConfig> >
{
public:
    //typedef typename Fibre<Index, Sttd64Text>::Type             TText;
    //typedef String<typename Value<TText>::Type/*, CONF_SW_EXT_MMAP<_TStrConf> */> _TIntText;
    typedef _TText                                              TText;
    typedef Index<TText, IndexSttd64<TConfig> >                 _TIndex;
    typedef Holder<TText>                                       THolder;

    typedef typename TConfig::template SufBuffConfig<>::Type   _TSufConf;
    typedef typename TConfig::template SufBuffPartitioningConfig<>::Type _TSufPartiConf;
    typedef typename TConfig::template TreeBuffConfig<>::Type  _TTreeConf;
    typedef typename TConfig::template TempBuffConfig<>::Type  _TTempConf;
    typedef typename TConfig::template StrBuffConfig<>::Type   _TStrConf;
    // Defining as below might lead to size _uint64 for each partition, which is useless, as only _uint32 length strings are the input
    //typedef typename Position<_TIntText>::Type                  TPosition;
    typedef unsigned                                            TPosition;
    typedef String<TPosition, CONF_SW_EXT_MMAP<_TSufConf> >     TOneWork;
    typedef String<TPosition, CONF_SW_EXT_MMAP<_TTempConf> >    _TTemp;
    typedef String<TPosition, CONF_SW_EXT_MMAP<_TSufPartiConf> > TOnePartSuffices;
    typedef String<Sttd64Node_, CONF_SW_EXT_MMAP<_TTreeConf> >  TPartTree;
    typedef typename Fibre<_TIndex, Sttd64Partitions>::Type     TSufficesParts;
    typedef typename Fibre<_TIndex, Sttd64SuffixTrees>::Type    TTrees;
    typedef typename Fibre<_TIndex, Sttd64VirtualTree>::Type    VirtualTree;
    typedef String<TOneWork>                                    TWorks;
    static const unsigned PREFIX_LENGTH = TConfig::PREFIX_LENGTH;

    //_TIntText txt;
    THolder txt;
    TSufficesParts partitions;
    TWorks workPartitions;
    mutable TTrees trees;   //as this needs to be loaded and closed at any time while iterating...
    std::string dirName; //the directory where all externally mapped strings are stored
    VirtualTree virtualPrefixTree;
    bool leafDepthsFilledIn;
    bool inMem;

    char const * _treeNamePart;
    char const * _sufficesNamePart;
    char const * _stringNamePart;
    char const * _tempNamePart;

    Index()
        :
            virtualPrefixTree(0, 0)
    {
    }
    
    template <typename TText_>
    Index(TText_ & text):
            /*txt(text),*/dirName(""),
            virtualPrefixTree((unsigned) ValueSize<typename Value<TText>::Type>::VALUE, PREFIX_LENGTH),
            leafDepthsFilledIn(false), inMem(true)
    {
        ConstrInit();
        setValue(this->txt, text);
    }

    Index(Index & other) :
            txt(other.txt),
            partitions(other.partitions),
            trees(other.trees),
            virtualPrefixTree(other.virtualPrefixTree),
            leafDepthsFilledIn(other.leafDepthsFilledIn),
            inMem(other.inMem)
    {
        ConstrInit();
        std::cerr << "Unexpected.";
    }
    Index(Index const & other) :
            txt(other.txt),
            partitions(other.partitions),
            trees(other.trees),
            virtualPrefixTree(other.virtualPrefixTree),
            leafDepthsFilledIn(other.leafDepthsFilledIn),
            inMem(other.inMem)
    {
        ConstrInit();
    }
    
private:
    inline void ConstrInit()
    {
        _treeNamePart = ".suffix_tree.";
        _sufficesNamePart = ".suffix_part.";
        _stringNamePart = ".string_cpy.";
        _tempNamePart = ".temp_buff.";
    }

};

// TODO(krugel) save/open

// ----------------------------------------------------------------------------
// save/open
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
inline bool
configureSave(Index<TText, IndexSttd64<TSpec> > & index, const char * fileName) {
    index.dirName = fileName;
    index.inMem = false;
    return true;
}

template <typename TText, typename TSpec>
inline bool
save(Index<TText, IndexSttd64<TSpec> > & index, const char * fileName, int openMode) {
	String<char> name;

    name = fileName; append(name, ".txt");
    if (! save(getFibre(index, Sttd64Text()), toCString(name), openMode)) return false;

    index.dirName = fileName;
    index.inMem = false;
    return true;
}

template <typename TText, typename TSpec>
inline bool
save(Index<TText, IndexSttd64<TSpec> > & index, const char * fileName) {
    return save(index, fileName, OPEN_WRONLY | OPEN_CREATE);
}

template <typename TText, typename TSpec>
inline bool
open(Index<TText, IndexSttd64<TSpec> > & index, const char * fileName, int openMode) {
    clear(index);

    String<char> name;

    name = fileName; append(name, ".txt");
    if (! open(getFibre(index, Sttd64Text()), toCString(name), openMode)) return false;

    index.dirName = fileName;
    index.inMem = false;
    return true;
}

template <typename TText, typename TSpec>
inline bool
open(Index<TText, IndexSttd64<TSpec> > & index, const char * fileName) {
    return open(index, fileName, OPEN_RDONLY);
}

// ============================================================================
// Metafunctions
// ============================================================================

template<typename TText, typename TSpec>
struct Fibre<Index<TText, IndexSttd64<TSpec> >, Sttd64SuffixTrees>
{
    typedef Index<TText, IndexSttd64<TSpec> > _TIndex;
    typedef String<typename _TIndex::TPartTree> Type;
};

template<typename TText, typename TSpec>
struct Fibre<Index<TText, IndexSttd64<TSpec> >, Sttd64VirtualTree>
{
    typedef Index<TText, IndexSttd64<TSpec> > _TIndex;
    typedef _VirtualPrefixTree<typename Value<typename _TIndex::TText>::Type>  Type;
};

template<typename TText, typename TSpec>
struct Fibre<Index<TText, IndexSttd64<TSpec> >, Sttd64Partitions>
{
    typedef Index<TText, IndexSttd64<TSpec> > _TIndex;
    typedef String<typename _TIndex::TOnePartSuffices> Type;
};

template<typename TText, typename TSpec>
struct Fibre<Index<TText, IndexSttd64<TSpec> >, Sttd64Text>
{
    typedef Index<TText, IndexSttd64<TSpec> > _TIndex;
    typedef typename _TIndex::TText Type;
};

template <typename TText, typename TSpec>
struct SAValue<Index<TText, IndexSttd64<TSpec> > >
{
    typedef Index<TText, IndexSttd64<TSpec> > _TIndex;
    typedef typename Size<typename _TIndex::TText>::Type Type;
};

// TODO(krugel) Enable as soon as suffix tree iterators with parent links are implemented
template <typename TSpec>
struct SupportsSuffixTreeIteration<IndexSttd64<TSpec> > {
    typedef False Type;
};

// ============================================================================
// Functions
// ============================================================================

template<typename TText_, typename TConfig>
inline bool createPartitions(Index<TText_, IndexSttd64<TConfig> > & index)
{
    SEQAN_CHECKPOINT;
    typedef Index<TText_, IndexSttd64<TConfig> > TIndex;
    typedef typename Fibre<TIndex, Sttd64Text>::Type TText;

    unsigned const alphabetSize = (unsigned) ValueSize<typename Value<TText>::Type>::VALUE;
    unsigned const prefLen = TIndex::PREFIX_LENGTH;
    unsigned numPartitions = static_cast<unsigned int>(pow(static_cast<double>(alphabetSize), static_cast<double>(prefLen)));
    typename Iterator<TText const>::Type textit;
    unsigned partitionPosition;

    //create the necessary partitions
    if (!_createPartitionFiles(index, index.dirName, numPartitions))
    {
        std::cerr << "Could not create partition files!";
        return false;
    }

    //distribute suffices across partitions
    if (numPartitions > 1)
    {
        unsigned i = 0;
        for (i = 0, textit = begin(indexText(index)); i < length(indexText(index)) - prefLen + 1; ++textit, ++i)
        {
            SEQAN_ASSERT(textit != end(indexText(index)));
            _getPartitionNumber(partitionPosition, index, infix(indexText(index), i, i + prefLen));
            typename TIndex::TOnePartSuffices & partitionForPrefix = value(index.partitions, partitionPosition);
            appendValue(partitionForPrefix, i + prefLen);
        }

        typename Iterator<typename TIndex::TSufficesParts>::Type partIt;
        unsigned j;
        for (j = 0, partIt = begin(index.partitions); partIt != end(index.partitions);++j, ++partIt)
        {
            if (!empty(*partIt))
            {
                index.virtualPrefixTree.setPartitionHasValues(j);
            }

            _unloadFibreFile(*partIt, index.inMem, j, index.dirName, index._sufficesNamePart);
        }

        index.virtualPrefixTree.setPrefLenSuffix(suffix(indexText(index), i));
    }
    else
    {
        SEQAN_ASSERT_EQ(prefLen, 0u);
        typename TIndex::TOnePartSuffices & partitionForPrefix = value(index.partitions, 0);
        unsigned i;
        for (i = 0, textit = begin(indexText(index)); i < length(indexText(index)); ++textit, ++i)
        {
            SEQAN_ASSERT(textit != end(indexText(index)));
            appendValue(partitionForPrefix, i);
        }
        index.virtualPrefixTree.prefLenSuffix = "";
        index.virtualPrefixTree.setPartitionHasValues(0);
    }

    return true;
}

template<typename TConfig, typename TText>
inline bool
_createPartitionFiles(Index<TText, IndexSttd64<TConfig> > & index,
                      std::string dirName,
                      unsigned numPartitions)
{
    SEQAN_CHECKPOINT;
    typedef Index<TText, IndexSttd64<TConfig> > TIndex;

    resize(index.partitions, numPartitions);
    resize(index.workPartitions, numPartitions);

    for (unsigned i = 0; i < numPartitions; i++)
    {

        typename TIndex::TOnePartSuffices & curPartString = value(index.partitions, i);
        if (!_openFibreFile(curPartString, index.inMem, i, dirName, index._sufficesNamePart, false, false))
        {
            std::cout << "Could not open partition file for partition " << i << "!";
            return false;
        }
    }
    return true;
}

template<typename _TText, typename TConfig>
inline bool createSuffixTrees(Index<_TText, IndexSttd64<TConfig> > & index)
{
    SEQAN_CHECKPOINT;
    typedef Index<_TText, IndexSttd64<TConfig> > TIndex;
    typedef typename Fibre<TIndex, Sttd64Text>::Type TText;
    unsigned const prefLen = TIndex::PREFIX_LENGTH;
    unsigned const alphabetSize = (unsigned) ValueSize<typename Value<TText>::Type>::VALUE;
    unsigned const sentinelCount = 1;

    unsigned lcp;

    //typedef typename Iterator<typename TIndex::TOnePartSuffices, Rooted>::Type TSuffPartIt;

    typedef typename Iterator<typename TIndex::TWorks, Rooted>::Type TWorksIt;

    TWorksIt currentWorkPart;

    typedef typename Iterator<typename TIndex::TOneWork, Rooted>::Type TOneWorkIt;

    //arrays for count sort, defined and initialized here in order
    //to avoid unnecessary (de)allocations
    unsigned * counts = new unsigned[alphabetSize + sentinelCount];
    TOneWorkIt* startPForLetter = new TOneWorkIt[alphabetSize + sentinelCount];
    unsigned* evaluationOrder = new unsigned[alphabetSize + sentinelCount];
    bool* hasBeenFound = new bool[alphabetSize + sentinelCount];
    String<typename Value<TText>::Type > lcpBuffer;
    resize(lcpBuffer, 2048);

    typename TIndex::_TTemp TEMP;
    if (!_openFibreFile(TEMP, index.inMem, 0, index.dirName, index._tempNamePart, false, false))
    {
        std::cerr << "Could not open TEMP buffer file.";
        return false;
    }

    String<typename Position<typename TIndex::TPartTree>::Type, Block<512> > stack;

    if (index.inMem) {
        resize(index.trees, index.virtualPrefixTree.numPartitions());
    } else {
        resize(index.trees, 1);
        resize(index.partitions, 1);
        resize(index.workPartitions, 1);
    }

    typename Iterator<typename TIndex::TSufficesParts>::Type partitionsIt;
    typename Iterator<typename TIndex::TTrees>::Type treesIt;
    unsigned treeLen;

    unsigned currentPartitionNumber = 0;
    for (currentPartitionNumber = 0, partitionsIt = begin(index.partitions), treesIt = begin(index.trees), currentWorkPart = begin(index.workPartitions);
            currentPartitionNumber < index.virtualPrefixTree.numPartitions();
            currentPartitionNumber++, (index.inMem ? (void)(++partitionsIt, ++treesIt, ++currentWorkPart) : (void)(partitionsIt)))
    {
        typename TIndex::TOnePartSuffices & curPartSuffixParti = *partitionsIt;
        typename TIndex::TOneWork & curPartSuffices = *currentWorkPart;

        _initWorkPart(curPartSuffices, curPartSuffixParti, currentPartitionNumber, index);

        if (empty(curPartSuffices))
        {
            _unloadWorkPart(curPartSuffices, curPartSuffixParti, currentPartitionNumber, index);

            continue;
        }
        typename TIndex::TPartTree & curPartTree = *treesIt;
        if (!_openFibreFile(curPartTree, index.inMem, currentPartitionNumber, index.dirName, index._treeNamePart, false, false))
        {
            std::cerr << "Could not open tree file (Partition: " << currentPartitionNumber
                    << ")." << std::endl;
            return false;
        }
        //special case, where tree only contains one leaf
        typename Size<typename TIndex::TOnePartSuffices>::Type curPartSuffLen;
        curPartSuffLen = length(curPartSuffices);
        if (curPartSuffLen == 1)
        {
            Sttd64Node_ leaf;
            leaf.lp = getValue(curPartSuffices, 0) - prefLen; //subtract preflen to get real start position
            leaf.pointerOrDepth = 0;
            leaf.setRightMostBit(true);
            leaf.setLeafBit(true);
            appendValue(curPartTree, leaf);
            //std::cout << "Tree (Partition Number: "<< currentPartitionNumber << ")" << std::endl;
            //_printPartTree(curPartTree, false);
            //std::cout << std::endl;
            _unloadWorkPart(curPartSuffices, curPartSuffixParti, currentPartitionNumber, index);
            _unloadFibreFile(curPartTree, index.inMem, currentPartitionNumber, index.dirName, index._treeNamePart);
            continue;
        }

        //Iterator con- and destructing slow for external strings
        //so instead of creating an infix give iterators to the right positions
        //maybe const's for further optis?? would need to look up internals
        TOneWorkIt beginSuffInf;
        TOneWorkIt endSuffInf;
        beginSuffInf = begin(curPartSuffices);
        endSuffInf = end(curPartSuffices);
        SEQAN_ASSERT_EQ(curPartSuffLen, (typename Size<typename TIndex::TOnePartSuffices>::Type) difference(beginSuffInf, endSuffInf));
        lcp = _longestCommonPrefix(indexText(index), beginSuffInf, endSuffInf, 0u, lcpBuffer);
        Sttd64Node_ root;
        root.lp = *beginSuffInf - prefLen;
        root.pointerOrDepth = 1;
        root.setRightMostBit(true);
        appendValue(curPartTree, root);
        treeLen = 1;

        resize(TEMP, curPartSuffLen);

        _evaluateNode(indexText(index), curPartTree, treeLen, beginSuffInf, endSuffInf, TEMP, stack, lcp,
                      counts, startPForLetter, evaluationOrder, hasBeenFound, alphabetSize,
                      sentinelCount);
        bool nodeWasRM;
        typename Position<typename TIndex::TPartTree>::Type nodePosInTree;
        //pop a node from the stack
        while (!empty(stack))
        {
            nodePosInTree = back(stack);

            //the following line (and the according usage) leads to a seg fault
            //after evaluateNode -> string size may have changed and addresses might have moved
            //=> address of reference becomes invalid (e.g. for MMap files this is probable to occurr)
            Sttd64Node_ node = getValue(curPartTree, nodePosInTree);
            nodeWasRM = node.isRightMost();
            eraseBack(stack);
            beginSuffInf = iter(curPartSuffices, node.lp);
            endSuffInf = iter(curPartSuffices, node.cleanPointerOrDepth());
            //the following write should be reflected (value instead of getValue used)! make sure in asserts
            //!!!! cannot use old reference in asserts (see above), so don't use reference at all
            SEQAN_ASSERT_EQ(length(curPartTree), treeLen);
            SEQAN_ASSERT_EQ(treeLen & 3 << 30, 0u);
            node.pointerOrDepth = treeLen;
            node.lp = *beginSuffInf;
            node.setRightMostBit(nodeWasRM);
            //write changes back instead
            assignValue(curPartTree, nodePosInTree, node);

            lcp = _longestCommonPrefix(indexText(index), beginSuffInf, endSuffInf, 1u, lcpBuffer);

            _evaluateNode(indexText(index), curPartTree, treeLen, beginSuffInf, endSuffInf, TEMP, stack,
                          lcp, counts, startPForLetter, evaluationOrder, hasBeenFound, alphabetSize,
                          sentinelCount);
            SEQAN_ASSERT_EQ(node.lp, value(curPartTree, nodePosInTree).lp);
            SEQAN_ASSERT_EQ(node.pointerOrDepth, value(curPartTree, nodePosInTree).pointerOrDepth);
        }
        //std::cout << "Tree (Partition Number: " << currentPartitionNumber << ")" << std::endl;
        //_printPartTree(curPartTree, false);
        //std::cout << std::endl;
        _unloadFibreFile(curPartTree, index.inMem, currentPartitionNumber, index.dirName, index._treeNamePart);
        _unloadWorkPart(curPartSuffices, curPartSuffixParti, currentPartitionNumber, index);
    }   //for

    _unloadFibreFile(TEMP, index.inMem, 0, index.dirName, index._tempNamePart);

    //clean up the count arrays
    delete[] counts;
    counts = 0;
    delete[] startPForLetter;
    startPForLetter = 0;
    delete[] evaluationOrder;
    evaluationOrder = 0;
    delete[] hasBeenFound;
    hasBeenFound = 0;
    return true;
}

template <typename TWork, typename TParti, typename TIndex>
inline bool _initWorkPart(TWork & works, TParti & parti, unsigned currentPartitionNumber,
                          TIndex & index)
{
    //if persistent storage was enabled, the partition files have been closed after partitioning
    if (index.inMem)
    {
        //no works opened yet
        if (!_openFibreFile(works, index.inMem, currentPartitionNumber, index.dirName,
                            (std::string(index._sufficesNamePart) + "wx").c_str(),
                            false, false))
        {
            std::cerr << "Could not reopen partition file (Partition: " << currentPartitionNumber
                    << ")." << std::endl;
            return false;
        }
        assign(works, parti);
        return true;
    }

    _unloadFibreFile(parti,
                     index.inMem,
                     currentPartitionNumber,
                     index.dirName, index._sufficesNamePart);

    //if persistent storage was enabled, the partition files have been closed after partitioning
    if (!_openFibreFile(works, index.inMem, currentPartitionNumber, index.dirName,
                        index._sufficesNamePart,
                        true, false))
    {
        std::cerr << "Could not reopen partition file (Partition: " << currentPartitionNumber
                << ")." << std::endl;
        return false;
    }
    return true;
}

template <typename TWork,typename TParti, typename TIndex>
inline void _unloadWorkPart(TWork & works, TParti & parti, unsigned currentPartitionNumber,
                           TIndex & index)
{
    if (index.inMem) {
        _unloadFibreFile(works,
                         index.inMem,
                         currentPartitionNumber, index.dirName,
                         (std::string(index._sufficesNamePart) + "wx"));
        _unloadFibreFile(parti,
                         index.inMem,
                         currentPartitionNumber,
                         index.dirName, index._sufficesNamePart);
        return;
    }

    _unloadFibreFile(works,
                     index.inMem,
                     currentPartitionNumber,
                     index.dirName, index._sufficesNamePart);

}

template<typename _TText, typename TConfig>
inline bool _createLeafDepths(Index<_TText, IndexSttd64<TConfig> > & index)
{
    //TODO(aumann): use Iterator<TopDown<ParentLink> > as soon as implemented
    if (index.leafDepthsFilledIn) return true;

    typedef typename Index<_TText, IndexSttd64<TConfig> >::TTrees TTrees;
    typename Iterator<TTrees>::Type treesIt;
    unsigned i;
    for (i = 0, treesIt = begin(index.trees);
            i < index.virtualPrefixTree.numPartitions();
            ++i, index.inMem ? ++treesIt : treesIt)
    {
        SEQAN_ASSERT(treesIt != end(index.trees));
        if(!index.virtualPrefixTree.partitionHasValues(i)) continue;
        if (!_reopenTreeFileRW(index, i))
        {
            std::cerr << "Could not reopen tree file (Partition: " << i
                    << ")." << std::endl;
            return false;
        }

        if (!empty(*treesIt))
        {
            _tempFillPartTreeLeafDepths(*treesIt);
        }
        _unloadTreeFile(index, i);
    }
    index.leafDepthsFilledIn = true;
    return true;
}

//opens the tree or partition file
//determined by the name given (e.g. "suffix_tree", "suffix_part")
template<typename TFibre>
inline bool _openFibreFile(TFibre & openedFile,
                           bool inMem,
                           unsigned const partitionNumber,
                           std::string const dirName,
                           std::string const fibreNamePart,
                           bool reOpen,
                           bool readOnly)
{
    SEQAN_CHECKPOINT;
    std::ostringstream ostr;
    ostr << partitionNumber;
    std::string currentPartitionSuffix = ostr.str();

    SEQAN_ASSERT(!readOnly || reOpen);  //convention
    if (inMem)
    {
        if (!reOpen)
        {
            if (!openTemp(openedFile))
                return false;
        }
    }
    else
    {
        flush(openedFile);
        close(openedFile);
        if (reOpen)
        {
            if (!readOnly)
            {
                //directly mapping the string RDWR does not seem to work. bug??
                TFibre tmpFile;
                std::remove((dirName + fibreNamePart + currentPartitionSuffix + "RO").c_str());
                int renRes = std::rename( (dirName + fibreNamePart + currentPartitionSuffix).c_str(),
                                     (dirName + fibreNamePart + currentPartitionSuffix + "RO").c_str());
                if (renRes != 0)
                {
                    std::perror("Could not rename");
                    return false;
                }
                if (!open(tmpFile, (dirName + fibreNamePart + currentPartitionSuffix + "RO").c_str(), OPEN_RDONLY))
                {
                    return false;
                }
                if (!open(openedFile, (dirName + fibreNamePart + currentPartitionSuffix).c_str(),
                          OPEN_RDWR | OPEN_CREATE))
                {
                    return false;
                }
                resize(openedFile, length(tmpFile), Exact());
                assign(openedFile, tmpFile);
                close(tmpFile);
                std::remove((dirName + fibreNamePart + currentPartitionSuffix + "RO").c_str());
            }
            else
            {
                if(!open(openedFile, (dirName + fibreNamePart + currentPartitionSuffix).c_str(),
                    OPEN_RDONLY)){
                    return false;
                }
            }
        }
        else
        {
            std::remove((dirName + fibreNamePart + currentPartitionSuffix).c_str());
            if (!open(openedFile, (dirName + fibreNamePart + currentPartitionSuffix).c_str(),
                      OPEN_RDWR | OPEN_CREATE))
                return false;
        }
    }
    return true;

}

template<typename TFibre>
inline bool _unloadFibreFile(TFibre & closedFile,
                             bool inMem,
                             unsigned const partitionNumber,
                             std::string const dirName,
                             std::string const fibreNamePart)
{
    SEQAN_CHECKPOINT;
    if (!inMem)
    {
        //closes the fibre file only iff a dirName was given (i.e. not only temp files used
        //which would then be lost)
        std::ostringstream ostr;
        ostr << partitionNumber;
        std::string currentPartitionSuffix = ostr.str();

        //flush(closedFile);
        if (!close(closedFile)) {
            std::cerr << "Could not close file: " << dirName + fibreNamePart + currentPartitionSuffix << std::endl;
            return false;
        }
        clear(closedFile);
    }
    return true;
}

template <typename TIndex>
inline bool _reopenTreeFile(TIndex const & index, unsigned partitionNumber)
{
   return _openFibreFile(value(index.trees, 0),
                         index.inMem,
                         partitionNumber,
                         index.dirName,
                         index._treeNamePart,
                         true,
                         true);
}

template <typename TIndex>
inline bool _reopenTreeFileRW(TIndex const & index, unsigned partitionNumber)
{
    if (!index.inMem)
    {
        return _openFibreFile(value(index.trees, 0),
                                     index.inMem,
                                     partitionNumber,
                                     index.dirName,
                                     index._treeNamePart,
                                     true,
                                     false);
    }

   return _openFibreFile(value(index.trees, partitionNumber),
                         index.inMem,
                         partitionNumber,
                         index.dirName,
                         index._treeNamePart,
                         true,
                         false);
}

template <typename TIndex>
inline bool _unloadTreeFile(TIndex const & index, unsigned partitionNumber)
{
    if (!index.inMem)
    {
        return _unloadFibreFile(value(index.trees, 0), index.inMem, partitionNumber,
                                index.dirName, index._treeNamePart);
    }

    return _unloadFibreFile(value(index.trees, partitionNumber), index.inMem, partitionNumber,
                            index.dirName, index._treeNamePart);
}

template<typename TText, typename TTree, typename Stack, typename TPartSuffixIt,
            typename TTemp>
inline void _evaluateNode(TText const & text,
                          TTree & tree,
                          unsigned & treeLen,
                          TPartSuffixIt & beginSuffInf,
                          TPartSuffixIt & endSuffInf,
                          TTemp & TEMP,
                          Stack & stack,
                          unsigned lcp,
                          unsigned counts[],
                          TPartSuffixIt startPForLetter[],
                          unsigned evaluationOrder[],
                          bool hasBeenFound[],
                          unsigned const alphabetSize,
                          unsigned const sentinelCount)
{

    _countSort(text, beginSuffInf, endSuffInf, /*suffInfLen,*/ TEMP, lcp, counts, startPForLetter, evaluationOrder,
               hasBeenFound, alphabetSize, sentinelCount);
    unsigned i = 0;
    for (unsigned curEval = 0; curEval < alphabetSize + sentinelCount; ++curEval)
    {
        i = evaluationOrder[curEval];
        if (i == std::numeric_limits<unsigned>::max())
            continue; //has not been found (no children)

        SEQAN_ASSERT_GEQ(counts[i], 1u);
        if (counts[i] == 1)
        { //leaf
            Sttd64Node_ leaf;
            leaf.lp = *startPForLetter[i];
            leaf.pointerOrDepth = 0;
            leaf.setLeafBit(true);
            appendValue(tree, leaf);
            ++treeLen;
        }
        else
        {
            Sttd64Node_ nodeMarker;
            nodeMarker.lp = position(startPForLetter[i]);   //same as * above
            nodeMarker.pointerOrDepth = position(startPForLetter[i]) + counts[i];
            appendValue(tree, nodeMarker);
            appendValue(stack, treeLen++);  //new length - 1 = position of last node = post++
        }
    }
    //as all children of the node have now been added, the last one must be the rightmost
    Sttd64Node_ & lastNode = value(tree, treeLen - 1);
    lastNode.setRightMostBit(true);
}

// TODO(krugel) Move ... ?
template<typename TString>
inline typename Position<ExtStringIterator<TString> >::Type
position(ExtStringIterator<TString>  const & it)
{
    return it.pageOfs + it.pageNo * it.PAGESIZE;
}

template<typename TString>
inline typename Position<ExtStringConstIterator<TString> >::Type
position(ExtStringConstIterator<TString>  const & it)
{
    return it.pageOfs + it.pageNo * it.PAGESIZE;
}

template<typename TString>
inline typename Position<ExtStringFwdIterator<TString> >::Type
position(ExtStringFwdIterator<TString>  const & it)
{
    return it.pageOfs + it.pageNo * it.PAGESIZE;
}

template<typename TString>
inline typename Position<ExtStringFwdConstIterator<TString> >::Type
position(ExtStringFwdConstIterator<TString>  const & it)
{
    return it.pageOfs + it.pageNo * it.PAGESIZE;
}

template<typename TText, typename TSufficesIt, typename TTextBuffer>
inline unsigned _longestCommonPrefixBlock(TText const & text,
                                     TSufficesIt & beginSuffInf,
                                     TSufficesIt & endSuffInf,
                                     unsigned lcpOverall,
                                     TTextBuffer & curLcpTextCopy) //buffer allocated outside of loops
{
    SEQAN_CHECKPOINT;

    unsigned maxLcpCurrentBlock = 2;
    unsigned currLcpCurrentBlock;

    typename Difference<TSufficesIt>::Type suffInfLen = difference(beginSuffInf, endSuffInf);
    SEQAN_ASSERT_EQ(suffInfLen, (typename Difference<TSufficesIt>::Type) position(endSuffInf) - (typename Difference<TSufficesIt>::Type) position(beginSuffInf));
    SEQAN_ASSERT_GEQ(suffInfLen, 0);
    typename Iterator<TText const>::Type textIt;
    typename Iterator<TText const>::Type textEnd = end(text);

    unsigned posInText;
    typedef typename Position<TTextBuffer>::Type TBufferPos;
    TBufferPos bufferPos;

    bool repeat = false;
    do {
        currLcpCurrentBlock = maxLcpCurrentBlock;
        posInText = *beginSuffInf + lcpOverall;
        textIt = iter(text, posInText);
        for (bufferPos = 0;
                textIt != textEnd && bufferPos < maxLcpCurrentBlock;
                ++textIt, ++bufferPos)
        {
            assignValue(curLcpTextCopy, bufferPos, *textIt);
        }
        SEQAN_ASSERT_LEQ(bufferPos, maxLcpCurrentBlock);
        if (bufferPos != maxLcpCurrentBlock)
            currLcpCurrentBlock = bufferPos;    //if text not long enough

        ++beginSuffInf;
        while (beginSuffInf != endSuffInf && currLcpCurrentBlock > 0)
        {
            posInText = *beginSuffInf + lcpOverall;
            textIt = iter(text, posInText);

            for (bufferPos = 0;
                    textIt != textEnd && bufferPos < currLcpCurrentBlock;
                    ++bufferPos, ++textIt)
            {
                if (getValue(curLcpTextCopy, bufferPos) != *textIt) break;
            }

            SEQAN_ASSERT(lcpOverall > 0 || bufferPos >= 1u); //minimum of matches is the first character
            currLcpCurrentBlock = bufferPos;    //if text not long enough

            ++beginSuffInf;
        }

        SEQAN_ASSERT_GEQ(difference(beginSuffInf, endSuffInf), 0);
        SEQAN_ASSERT_GEQ(suffInfLen - difference(beginSuffInf, endSuffInf), 0);
        SEQAN_ASSERT_EQ(difference(beginSuffInf, endSuffInf), (typename Difference<TSufficesIt>::Type) position(endSuffInf) - (typename Difference<TSufficesIt>::Type) position(beginSuffInf));
        //typedef typename Size<typename Container<TSufficesIt>::Type >::Type TDiffUnsigned;
        beginSuffInf -= suffInfLen - difference(beginSuffInf, endSuffInf); //reset to first suffix
        lcpOverall += currLcpCurrentBlock;
        SEQAN_ASSERT_LEQ(currLcpCurrentBlock, maxLcpCurrentBlock);
        if (currLcpCurrentBlock == maxLcpCurrentBlock)
        {
            maxLcpCurrentBlock *= 2;    //TODO(aumann): check what would be good? 2 seems fine though
            //loadFactor *= 2;
            if (maxLcpCurrentBlock > 2048)
                maxLcpCurrentBlock = 2048;
            repeat = true;
        }
        else
        {
            repeat = false;
        }

    } while (repeat);

    SEQAN_ASSERT_GEQ(lcpOverall, 1u);
    return lcpOverall;
}

template<typename TText, typename TSufficesIt, typename TTextBuff>
inline unsigned _longestCommonPrefix(TText const & text,
                                     TSufficesIt & beginSuffInf,
                                     TSufficesIt & endSuffInf,
                                     unsigned minLCPGuaranteed,
                                     TTextBuff & textBuffer)
{
    unsigned startOfFirstSuffix;
    unsigned startOfCompSuffix;
    startOfFirstSuffix = *beginSuffInf;
    typename Value<TText>::Type charFirstSuffix;
    typename Value<TText>::Type compareChar;
    typename Size<TText>::Type textLen = length(text);
    typename Size<typename Container<TSufficesIt>::Type >::Type movedFromBeginSuffInf = 0;
    //as the lcp is at least minLCPGuaranteed, we start at start position + minLCPGuaranteed

    //in a one-string scenario, the text length can be used as sentinel
    SEQAN_ASSERT_LEQ(startOfFirstSuffix + minLCPGuaranteed, textLen);
    if (startOfFirstSuffix + minLCPGuaranteed == textLen)
    {
        return minLCPGuaranteed;    // = minLCPGuaranteed
    }

    charFirstSuffix = getValue(text, startOfFirstSuffix + minLCPGuaranteed);
    for (++beginSuffInf, ++movedFromBeginSuffInf; beginSuffInf != endSuffInf; ++beginSuffInf, ++movedFromBeginSuffInf)
    {
        startOfCompSuffix = (*beginSuffInf) + minLCPGuaranteed;
        SEQAN_ASSERT_LEQ(startOfCompSuffix, textLen);
        if (startOfCompSuffix == textLen)
        {
            beginSuffInf -= movedFromBeginSuffInf;
            return minLCPGuaranteed;
        }
        compareChar = getValue(text, startOfCompSuffix);
        if (charFirstSuffix != compareChar)
        {
            beginSuffInf -= movedFromBeginSuffInf;
            return minLCPGuaranteed;
        }
    }
    beginSuffInf -= movedFromBeginSuffInf;

    //if control reaches this area the block lcp ist started (lcp of at least minLcp + 1 = 2 (1 for 1st call overall, and even then should be rare for larger text)
    return _longestCommonPrefixBlock(text, beginSuffInf, endSuffInf, minLCPGuaranteed + 1, textBuffer);
}

template<typename TText, typename TPartSuffixIt, typename TTemp>
inline void _countSort(TText const & text,
                       TPartSuffixIt & beginSuffInf,
                       TPartSuffixIt & endSuffInf,
                       TTemp & TEMP,
                       const unsigned prefLenOrLCP,
                       unsigned counts[],
                       TPartSuffixIt startPForLetter[],
                       unsigned evaluationOrder[],
                       bool hasBeenFound[],
                       unsigned const alphabetSize,
                       unsigned const sentinelCount)
{
    typedef unsigned Position;
    typename Iterator<TTemp, Rooted>::Type tempIt;
    unsigned curLetterVal;
    typename Size<TText>::Type textLen = length(text);
    SEQAN_ASSERT_GEQ(difference(beginSuffInf, endSuffInf), 2);
    for (unsigned i = 0; i < alphabetSize + sentinelCount; ++i)
    {
        counts[i] = 0;
        hasBeenFound[i] = false;
    }

    Position curSuffixPos;
    typename Size<typename Container<TPartSuffixIt>::Type >::Type lenPartSuffices = difference(beginSuffInf, endSuffInf);
    SEQAN_ASSERT_GEQ(lenPartSuffices, 2u);
    SEQAN_ASSERT_EQ(lenPartSuffices, position(endSuffInf) - position(beginSuffInf));
    SEQAN_ASSERT_GEQ(length(TEMP), lenPartSuffices);
    tempIt = begin(TEMP);
    unsigned nextEvaluationPos = 0;
    //first step: count letter occurrences
    bool allHaveBeenFound = false;
    Position i;
    for (i = 0; i < lenPartSuffices - 1 && !allHaveBeenFound; ++i)
    {
        SEQAN_ASSERT(beginSuffInf != endSuffInf);
        curSuffixPos = *beginSuffInf + prefLenOrLCP;
        //if we are one position after the end of the text, a sentinel is found(one-sentinel scenario)
        SEQAN_ASSERT_LT(curSuffixPos, textLen); //equality only possible at the very end!
        curLetterVal = ordValue(getValue(text, curSuffixPos));

        ++counts[curLetterVal];
        if (!hasBeenFound[curLetterVal])
        {
            hasBeenFound[curLetterVal] = true;
            evaluationOrder[nextEvaluationPos] = curLetterVal;
            ++nextEvaluationPos;
            if (nextEvaluationPos == alphabetSize)
            {
                allHaveBeenFound = true;    //go to optimized loop the one without the comparison
            }
        }
        assignValue(tempIt, curSuffixPos);
        ++beginSuffInf;
        ++tempIt;
    }

    //optimized second loop which does not test the hasBeenFound criteria
    for (/*i from above*/; i < lenPartSuffices - 1; ++i)
    {
        SEQAN_ASSERT(beginSuffInf != endSuffInf);
        curSuffixPos = *beginSuffInf + prefLenOrLCP;
        //if we are one position after the end of the text, a sentinel is found(one-sentinel scenario)
        SEQAN_ASSERT_LT(curSuffixPos, textLen);    //can only be last suffix, if always increasing order
        curLetterVal = ordValue(getValue(text, curSuffixPos));

        ++counts[curLetterVal];
        assignValue(tempIt, curSuffixPos);
        ++beginSuffInf;
        ++tempIt;
    }

    curSuffixPos = *beginSuffInf + prefLenOrLCP;
    SEQAN_ASSERT_LEQ(curSuffixPos, textLen);
    if (curSuffixPos == textLen)
    {
        curLetterVal = alphabetSize;    //first value no real character has
        hasBeenFound[alphabetSize] = true;
        //do not set the evaluation order for the sentinel here, as its always the last one
    }
    else
    {
        curLetterVal = ordValue(getValue(text, curSuffixPos));
        if (!hasBeenFound[curLetterVal])
        {
            hasBeenFound[curLetterVal] = true;
            evaluationOrder[nextEvaluationPos] = curLetterVal;
            ++nextEvaluationPos;
        }
    }
    assignValue(tempIt, curSuffixPos);
    ++counts[curLetterVal];

    beginSuffInf -= (lenPartSuffices - 1);  //the -1 is as the last iteration was not executed

    for (unsigned j = nextEvaluationPos; j < alphabetSize + sentinelCount; ++j)
    {
        evaluationOrder[j] = std::numeric_limits<unsigned>::max(); //mark as invalid, if some letters were not found
    }

    //sort in at the evaluation order positions
    //not necessary for the algo to work correctly (? think about this again)
    //but may slightly increase disk locality - this it does, so this outweighs alphabetical ordering?
    //Alphabet ordering (except for first child) would cost +5% construction time but may well be worth that.
    unsigned sumBefore = 0;
    SEQAN_ASSERT_EQ(evaluationOrder[alphabetSize], std::numeric_limits<unsigned>::max());
    for (unsigned j = 0; j < alphabetSize; ++j)
    {
        if (evaluationOrder[j] == std::numeric_limits<unsigned>::max())
        {
            break; //no more suffices
        }
        SEQAN_ASSERT_LT(evaluationOrder[j], alphabetSize + sentinelCount);
        startPForLetter[evaluationOrder[j]] = beginSuffInf + sumBefore;
        sumBefore += counts[evaluationOrder[j]];
    }

    if (hasBeenFound[alphabetSize])
    {
        SEQAN_ASSERT_EQ(counts[alphabetSize], 1u);
        SEQAN_ASSERT_LT(nextEvaluationPos, alphabetSize + sentinelCount);
        //evaluationOrder[nextEvaluationPos] = alphabetSize; only at the very end
        startPForLetter[alphabetSize] = beginSuffInf + sumBefore;
    }

    for (tempIt = begin(TEMP), i = 0; i < lenPartSuffices; goNext(tempIt), ++i)
    {
        SEQAN_ASSERT(tempIt != end(TEMP));
        curSuffixPos = *tempIt;
        if (curSuffixPos == textLen)
        {
            SEQAN_ASSERT(hasBeenFound[alphabetSize]);
            curLetterVal = alphabetSize;
            //SEQAN_ASSERT_NEQ(startPForLetter[curLetterVal], std::numeric_limits<unsigned>::max());
            SEQAN_ASSERT(startPForLetter[curLetterVal] != endSuffInf);
            //assignValue(partSuffices, startPForLetter[curLetterVal]++, curSuffixPos);
            *(startPForLetter[curLetterVal]) = curSuffixPos;
        }
        else
        {
            curLetterVal = ordValue(getValue(text, curSuffixPos));
            //SEQAN_ASSERT_NEQ(startPForLetter[curLetterVal], std::numeric_limits<unsigned>::max());
            SEQAN_ASSERT(startPForLetter[curLetterVal] != endSuffInf);
            //assignValue(partSuffices, startPForLetter[curLetterVal]++, curSuffixPos);
            *(startPForLetter[curLetterVal]) = curSuffixPos;
            ++startPForLetter[curLetterVal];
        }

    }

    sumBefore = 0;
    unsigned j;
    for (j = 0; j < alphabetSize; ++j)
    {
        if (evaluationOrder[j] == std::numeric_limits<unsigned>::max())
        {
            break; //no more suffices
        }
        startPForLetter[evaluationOrder[j]] = beginSuffInf + sumBefore;
        sumBefore += counts[evaluationOrder[j]];
    }

    if (hasBeenFound[alphabetSize])
    {
        SEQAN_ASSERT_EQ(j, nextEvaluationPos);
        SEQAN_ASSERT_LEQ(nextEvaluationPos, alphabetSize + sentinelCount);
        evaluationOrder[nextEvaluationPos] = alphabetSize;
    }

}

template<typename TInfixPartSuffices>
inline void _printPositions(TInfixPartSuffices const & partSuffices)
{
    typename Iterator<TInfixPartSuffices const>::Type sufficesIt;
    std::cout << "Partition from: " << beginPosition(partSuffices) << " to " << endPosition(partSuffices) << std::endl;

    std::cout << "(";
    for (sufficesIt = begin(partSuffices); sufficesIt != end(partSuffices); ++sufficesIt)
    {
        std::cout << *sufficesIt << ",";
    }
    std::cout << ")" << std::endl;
}

template<typename TTree>
inline void _printRawTreeTable(TTree const & treeTable, std::ostream & out = std::cout)
{
    SEQAN_CHECKPOINT;
    typename Iterator<TTree const>::Type treeIt;
    Sttd64Node_ node;
    //out << "Tree table: ";
    for (treeIt = begin(treeTable); treeIt != end(treeTable); ++treeIt)
    {
        node = *treeIt;
        out << node.toString();
    }
    out << std::endl;
}

//recursive approach, iterators should be used instead when available, I guess?
//just for testing purposes
template<typename TTree>
inline void _tempFillPartTreeLeafDepths(TTree & tree)
{
    SEQAN_CHECKPOINT;
    Sttd64Node_ root;
    root = getValue(tree, 0);
    if (root.isLeaf())
    { //only a single leaf in the table
        //nothing to do, depth is zero
        SEQAN_ASSERT_EQ(root.cleanPointerOrDepth(), 0u);
    }
    else
    {
        //incomingEdgeLengthToRoot = getValue(tree, 1).lp - root.lp;
        //_printPartTreeRec(tree, 0, incomingEdgeLengthToRoot, printDepthsInLeaves);
        _tempFillPartTreeLeafDepths(tree, 0, 0);
    }
}

template<typename TTree>
inline void _tempFillPartTreeLeafDepths(TTree & tree,
                              typename Position<TTree>::Type currentFatherPos,
                              unsigned currentDepth)
{
    SEQAN_CHECKPOINT;
    typename Position<TTree>::Type childPos;
    Sttd64Node_ father;
    Sttd64Node_ child;

    father = getValue(tree, currentFatherPos);
    childPos = father.cleanPointerOrDepth();
    child = getValue(tree, childPos); //the length of the edge incoming to the
                                      //father is determined by the leftmost child node
    unsigned lengthIncomingEdgeFather = child.lp - father.lp;
    do
    {
        child = getValue(tree, childPos);
        if (child.isLeaf())
        {
            value(tree, childPos).pointerOrDepth += currentDepth + lengthIncomingEdgeFather;
            SEQAN_ASSERT_EQ(getValue(tree, childPos).cleanPointerOrDepth(), (currentDepth + lengthIncomingEdgeFather));
            //out << (child.lp - (currentDepth + lengthIncomingEdgeFather));
        }
        else
        {
            _tempFillPartTreeLeafDepths(tree, childPos, currentDepth + lengthIncomingEdgeFather);
        }
        ++childPos;
    } while (!child.isRightMost());

}
template<typename TTree>
inline void _printPartTree(TTree const & tree, bool printDepthsInLeaves, std::ostream & out = std::cout)
{
    SEQAN_CHECKPOINT;
    Sttd64Node_ root;
    root = getValue(tree, 0);
    if (root.isLeaf())
    { //only a single leaf in the table
        out << "(";
        out << root.lp;
        if (printDepthsInLeaves)
        { //depths have already been filled in, print to verify
            out << "-" << root.cleanPointerOrDepth();
        }
        out << ")";
    }
    else
    {
        _printPartTreeRec(tree, 0,  0, printDepthsInLeaves, out);
    }
    out << std::endl;
}

template<typename TTree>
inline void _printPartTreeRec(TTree const & tree,
                              typename Position<TTree>::Type currentFatherPos,
                              unsigned currentDepth,
                              bool printDepthsInLeaves,
                              std::ostream & out = std::cout)
{
    SEQAN_CHECKPOINT;
    //typedef typename Index<TText, IndexSttd64<TConfig> > TIndex;
    typename Position<TTree>::Type childPos;
    //typename Iterator<TIndex::TTree const> treeIt;
    Sttd64Node_ father;
    Sttd64Node_ child;

    father = getValue(tree, currentFatherPos);
    childPos = father.cleanPointerOrDepth();
    out << "(";
    child = getValue(tree, childPos); //the length of the edge incoming to the
                                      //father is determined by the leftmost child node
    unsigned lengthIncomingEdgeFather = child.lp - father.lp;
    do
    {
        child = getValue(tree, childPos);
        if (child.isLeaf())
        {
            out << "(";
            out << (child.lp - (currentDepth + lengthIncomingEdgeFather));
            if (printDepthsInLeaves)
            { //depths have already been filled in, print to verify
                out << "-" << child.cleanPointerOrDepth();
            }
            out << ")";
        }
        else
        {
            _printPartTreeRec(tree, childPos, currentDepth + lengthIncomingEdgeFather, printDepthsInLeaves, out);
        }
        ++childPos;
    } while (!child.isRightMost());
    out << ")";
}

template<typename TText, typename TConfig>
inline void clear(Index<TText, IndexSttd64<TConfig> > & index)
{
    SEQAN_CHECKPOINT;
    //typedef Index<TText, IndexSttd64<TConfig> > TIndex;
    //typename Iterator<typename TIndex::TSufficesParts>::Type partitionsIt;
    //typename Iterator<typename TIndex::TTrees>::Type treesIt;

    //TODO(aumann) FIXME: this does not work, how does deallocation work for strings of external strings
    //How do you destroy the iterators in a defined way?

    //first iterate over all partitions and trees and clear them
//        for (partitionsIt = begin(index.partitions); partitionsIt != end(index.partitions); ++partitionsIt)
//        {
//            clear(*partitionsIt);
//        }
    clear(index.partitions);

//        for (treesIt = begin(index.trees); treesIt != end(index.trees); ++treesIt)
//        {
//            clear(*treesIt);
//        }
    clear(index.trees);
    clear(index.txt);
}

template<typename TText_, typename TConfig, typename TSeq>
inline bool _getPartitionNumber(unsigned & partitionNum,
                                Index<TText_, IndexSttd64<TConfig> > const & index,
                                TSeq const & prefix)
{
    SEQAN_CHECKPOINT;
    typedef Index<TText_, IndexSttd64<TConfig> > TIndex;
    typedef typename Fibre<TIndex const, Sttd64Text>::Type TText;
    //typename Iterator<typename Infix<TText>::Type>::Type prefixIt;
    typename Iterator<typename TIndex::TTrees>::Type treesIt;
    //find first non-empty tree partition
    if (length(prefix) == 0)
    {
        SEQAN_ASSERT(false);
        //TODO(aumann): ASSERT tree has been built
        //WARNING: would not work anymore atm, as inMem != external trees
        for (treesIt = begin(index.trees), partitionNum = 0; treesIt != end(index.trees); ++treesIt, ++partitionNum)
        {
            if (!empty(*treesIt))
            {
                return true;
            }
        }
        return false; //no non-empty tree
    }
    else
    {
        //find correct partition for the given prefix

        unsigned const prefLen = TIndex::PREFIX_LENGTH;
        SEQAN_ASSERT_EQ(length(prefix), prefLen);
        unsigned const alphabetSize = (unsigned) ValueSize<typename Value<TText>::Type>::VALUE;
        //std::cerr << "alphabetSize=" << alphabetSize << std::endl;
        //unsigned const numPartitions = pow(static_cast<double>(alphabetSize), static_cast<double>(prefLen));
        _getPartitionNumber(partitionNum, prefix, alphabetSize, prefLen);
        return true;

    }
}
template <typename TSeq>
inline void _getPartitionNumber(unsigned & partitionNum, TSeq const & prefix,
                                unsigned alphabetSize, unsigned prefixLenTotal)
{
    unsigned numPartitions = static_cast<unsigned>(pow(static_cast<double>(alphabetSize), static_cast<double>(prefixLenTotal)));
    unsigned const prefixLen = length(prefix);
    partitionNum = 0;
    SEQAN_ASSERT_LEQ(prefixLen, prefixLenTotal);
    //typename Iterator<TSeq const>::Type prefixIt;
    //prefixIt = begin(prefix);
    for (unsigned j = 0; j < prefixLen; ++j/*, ++prefixIt*/)
    {
        SEQAN_ASSERT_EQ(numPartitions % alphabetSize, 0u);
        //SEQAN_ASSERT(prefixIt != end(prefix));
        partitionNum += (numPartitions / alphabetSize) * ordValue(getValue(prefix, j));
        numPartitions /= alphabetSize;
    }

    //unsigned const numPartitions = pow(static_cast<double>(alphabetSize), static_cast<double>(prefixLenTotal));
    //unsigned const prefixLen = length(prefix);
    //partitionNum = 0;
    //typename Iterator<TSeq>::Type prefixIt;
    //prefixIt = begin(prefix);
    //for (unsigned j = 0; j < prefixLen; ++j, ++prefixIt)
    //{
    //    SEQAN_ASSERT_EQ(numPartitions % (alphabetSize * (j + 1)), 0u);
    //    SEQAN_ASSERT(prefixIt != end(prefix));
    //    partitionNum += (numPartitions / (alphabetSize * (j + 1))) * ordValue(*prefixIt);
    //}
}

template <typename TSeq>
inline void _getPrefixForPartitionNumber(TSeq & prefix,
                                         unsigned partitionNumber,
                                         unsigned alphabetSize,
                                         unsigned prefixLenTotal,
                                         unsigned returnedPrefLen)
{
    SEQAN_CHECKPOINT;
    //unsigned const numPartitions = pow(static_cast<double>(alphabetSize), static_cast<double>(prefixLenTotal));
    SEQAN_ASSERT_LEQ(partitionNumber, static_cast<unsigned>(pow(static_cast<double>(alphabetSize), static_cast<double>(prefixLenTotal))));
    SEQAN_ASSERT_GEQ(returnedPrefLen, 0u);
    typename Value<TSeq>::Type letter;

    resize(prefix, returnedPrefLen);

    //first skip last letters, not relevant for requested returnedPrefLen long sequence
    for (unsigned i = prefixLenTotal - 1; i >= returnedPrefLen; i--) {
        partitionNumber = partitionNumber / alphabetSize;
    }

    for (int i = returnedPrefLen - 1; i >= 0; i--) {
        letter = partitionNumber % alphabetSize;
        assignValue(prefix, i, letter);
        partitionNumber = partitionNumber / alphabetSize;
    }
}

template <typename TText, typename TSpec>
void _indexRequireTopDownIteration(Index<TText, IndexSttd64<TSpec> > &index)
{
    indexRequire(index, Sttd64Partitions());
    indexRequire(index, Sttd64SuffixTrees());
    _createLeafDepths(index);
}

template <typename TText, typename TSpec, typename TSpecAlg>
inline bool indexCreate(Index<TText, TSpec> &index, Sttd64SuffixTrees, TSpecAlg const) {
    SEQAN_CHECKPOINT;
    //atm only one algo, so std creation routine is called
    return createSuffixTrees(index);
}

template <typename TText, typename TSpec, typename TSpecAlg>
inline bool indexCreate(Index<TText, TSpec> &index, Sttd64VirtualTree, TSpecAlg const) {
    SEQAN_CHECKPOINT;
    return indexRequire(index, Sttd64SuffixTrees());
}

template <typename TText, typename TSpec, typename TSpecAlg>
inline bool indexCreate(Index<TText, TSpec> &index, Sttd64Partitions, TSpecAlg const) {
    SEQAN_CHECKPOINT;
    //atm only one algo, so std creation routine is called
    return createPartitions(index);
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, TSpec> const, Sttd64VirtualTree>::Type &
getFibre(Index<TText, TSpec> const &index, Sttd64VirtualTree) {
    return index.virtualPrefixTree;
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, TSpec> const, Sttd64Partitions>::Type &
getFibre(Index<TText, TSpec> const &index, Sttd64Partitions) {
    return index.partitions;
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, TSpec> const, Sttd64SuffixTrees>::Type &
getFibre(Index<TText, TSpec> const &index, Sttd64SuffixTrees) {
    return index.trees;
}

//as soon as operating on StringSets, the text getFibres can be removed (index_shims.h should provide)
template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, IndexSttd64<TSpec> >, FibreText>::Type &
getFibre(Index<TText, IndexSttd64<TSpec> > &index, FibreText) {
    return value(index.txt);
}
template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, IndexSttd64<TSpec> > const, FibreText>::Type &
getFibre(Index<TText, IndexSttd64<TSpec> > const &index, FibreText) {
    return value(index.txt);
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, TSpec>, FibreRawText>::Type &
getFibre(Index<TText, IndexSttd64<TSpec> > &index, FibreRawText) {
    return value(index.txt);
}
template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, TSpec> const, FibreRawText>::Type &
getFibre(Index<TText, IndexSttd64<TSpec> > const &index, FibreRawText) {
    return value(index.txt);
}

} //namespace seqan

#endif // SEQAN_INDEX_INDEX_STTD64_BASE_H_
