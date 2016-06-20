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

#ifndef SEQAN_INDEX_INDEX_STTD64_STREE_H_
#define SEQAN_INDEX_INDEX_STTD64_STREE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template<typename TIndex, typename TPosition_>
struct VertexSttd64
{
    typedef TPosition_ TPosition;
    TPosition nodePos; // position of current node entry in directory
    Sttd64Node_ node;
    TPosition parentRepLen; // representative length of parent node
    TPosition edgeLen; // length of edge above current node
    unsigned partitionNum; // the number of the suffix tree containing this node
    // as there are AlphabetSize^PrefixLength partitions
    bool inVirtual; //if the incoming edge from the root to the iterator's node <= prefixLen
    bool isVirtualNode;
    String<typename Value<typename TIndex::TText>::Type > prefix;
    //std::vector<unsigned> prefLenSuffixMatches; //contains potential matches in the last prefLen suffices of the text
    //unsigned leftMostPrefPartition;
    //unsigned rightMostPrefPartition;
    VertexSttd64() :
        nodePos(0),
        parentRepLen(0),
        edgeLen(0),
        partitionNum(std::numeric_limits<unsigned>::max()),
        inVirtual(true),
        isVirtualNode(true)/*,
        leftMostPrefPartition(std::numeric_limits<unsigned>::max()),
        rightMostPrefPartition(std::numeric_limits<unsigned>::max())*/
    {
    }
    VertexSttd64(MinimalCtor) :
        nodePos(0),
        parentRepLen(0),
        edgeLen(0),
        partitionNum(std::numeric_limits<unsigned>::max()),
        inVirtual(true),
        isVirtualNode(true)/*,
        leftMostPrefPartition(std::numeric_limits<unsigned>::max()),
        rightMostPrefPartition(std::numeric_limits<unsigned>::max())*/
    {
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TText, typename TSpec>
struct Size<Index<TText, IndexSttd64<TSpec> > >
{
    typedef Index<TText, IndexSttd64<TSpec> > TIndex;
    typedef typename Size<typename TIndex::TText>::Type Type;
};

template<typename TText, typename TSpec>
struct VertexDescriptor<Index<TText, IndexSttd64<TSpec> > >
{
    typedef typename Size<Index<TText, IndexSttd64<TSpec> > >::Type TSize;
    typedef VertexSttd64<Index<TText, IndexSttd64<TSpec> >, TSize> Type;
};

template<typename TText, typename TSpec>
struct GetOccurrences<Index<TText, IndexSttd64<TSpec> > > {
    typedef Index<TText, IndexSttd64<TSpec> >               TIndex;
    typedef String<typename SAValue<TIndex>::Type>          Type;
};

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template<typename TText, typename TIndexSpec, typename TSpec>
class Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TSpec> >
{
public:
    typedef Index<TText, IndexSttd64<TIndexSpec> > TIndex;
    typedef typename VertexDescriptor<TIndex>::Type TVertexDesc;

    TIndex const *index; // container of all necessary tables
    TVertexDesc vDesc; // current interval in suffix array and

    TVertexDesc _parentDesc;        //to be able to goUp one node

    Iter(TIndex &_index) :
            index(&_index)
    {
        _indexRequireTopDownIteration(_index);
        goRoot(*this);
    }

    Iter(TIndex const &_index, MinimalCtor) :
            index(&_index), vDesc(MinimalCtor())
    {
    }

    Iter(Iter const &_origin) :
            index(&container(_origin)), vDesc(value(_origin))
    {
    }

};

template<typename TText, typename TIndexSpec, typename TSpec>
class Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TopDown<TSpec> > >
    : public Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TSpec> >
{
public:
    typedef Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TSpec> > TBase_;
    typedef typename TBase_::TIndex TIndex;
    typedef typename TBase_::TVertexDesc TVertexDesc;
    
    Iter(TIndex &_index) :
            TBase_(_index)
    {
        _indexRequireTopDownIteration(_index);
    }

    Iter(TIndex const &_index, MinimalCtor) :
            TBase_(_index, MinimalCtor())
    {
    }

    Iter(Iter const &_origin) :
        TBase_(_origin)
    {
    }

};

// TODO(krugel) So far only a stub, parent link functions not implemented yet
template<typename TText, typename TIndexSpec, typename TSpec>
class Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > >
    : public Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TopDown<TSpec> > >
{
public:
    typedef Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TopDown<TSpec> > > TBase_;
    typedef typename TBase_::TIndex TIndex;
    typedef typename TBase_::TVertexDesc TVertexDesc;
    
    //typedef Iter<TIndex, VSTree<TopDown<> > >		        TBase;
    //typedef	typename HistoryStackEntry_<Iter>::Type	        TStackEntry;
    //typedef String<TStackEntry, Block<> >			        TStack;

    //TStack			history;	// contains all previously visited nodes (allows to go up)

    Iter(TIndex &_index) :
            TBase_(_index)
    {
        _indexRequireTopDownIteration(_index);
    }

    Iter(TIndex const &_index, MinimalCtor) :
            TBase_(_index, MinimalCtor())
    {
    }

    Iter(Iter const &_origin) :
        TBase_(_origin)
    {
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

//template < typename TText, typename TIndexSpec, typename TSpec >
//struct HistoryStackEntry_<Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > >
//{
//    typedef Index<TText, IndexSttd64<TIndexSpec> >	        TIndex;
//    typedef typename Size<TIndex>::Type				        TSize;
//    typedef HistoryStackWotdModified_<TSize>		        Type;
//};

// ============================================================================
// Functions
// ============================================================================

template<typename TText, typename TIndexSpec, class TSpec>
inline void clear(Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TSpec> > &it)
{
    typedef typename VertexDescriptor<Index<TText, IndexSttd64<TIndexSpec> > >::Type VertexD;
    if (container(it).virtualPrefixTree.prefixLen > 0)
    {
        VertexD & valueIt = value(it);
        if (!value(it).inVirtual)
        {
            SEQAN_ASSERT_NEQ(value(it).partitionNum, std::numeric_limits<unsigned>::max());
            _unloadTreeFile(container(it), value(it).partitionNum);
        }
        valueIt = VertexD();
        SEQAN_ASSERT_EQ(value(it).parentRepLen, valueIt.parentRepLen);
        SEQAN_ASSERT_EQ(value(it).nodePos, valueIt.nodePos);
        SEQAN_ASSERT_EQ(value(it).edgeLen, valueIt.edgeLen);
        SEQAN_ASSERT_EQ(value(it).partitionNum, valueIt.partitionNum);
        //value(it) = valueIt;  //should not be necessary
    }
    else
    {
        VertexD childDesc;
        _getVertexOfPartitionFather(childDesc, 0, container(it));
        value(it) = childDesc;
    }
}

template<typename TText, typename TIndexSpec, class TSpec>
inline void goRoot(Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TSpec> > &it)
{
    clear(it);
}

template<typename TText, typename TIndexSpec, class TSpec>
inline bool isRoot(Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TSpec> > const &it)
{
    if (value(it).inVirtual && length(value(it).prefix) == 0)
        return true;
    if (Index<TText, IndexSttd64<TIndexSpec> >::PREFIX_LENGTH < 1 && value(it).nodePos == 0)
    {
        return true;
    }

    return false;
}

//    template <typename TIndex, class TSpec>
//    inline bool goDown(Iter<TIndex, VSTree<TopDown<TSpec> > > &it) {
//        if (_goDown(it, typename GetVSTreeIteratorTraits<Iter<TIndex, VSTree<TopDown<TSpec> > > >::Type())) {
//            _onGoDown(it);
//            return true;
//        } else
//            return false;
//    }

//TODO(aumann): _goDowns implementieren, dann sollte das funktionieren
//TODO(aumann): _goRights implementieren
//TODO(aumann): _getNodeByChar implementieren
//TODO(aumann): _historyClear/_historyPush implementieren für history TD-iterator
//TODO(aumann): _isLeafs implementieren
//TODO(aumann): lca, lcp
//TODO(aumann): parentEdge... -Sachen

// go down the leftmost edge (including empty $-edges)
/*
template<typename TText, class TIndexSpec, class TSpec, typename TDfsOrder>
inline bool _goDown(Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
    VSTreeIteratorTraits<TDfsOrder, False> const)
*/

template<typename TText, class TIndexSpec, class TSpec>
inline bool _goDown(Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
    DeepestSpec<VSTree<TopDown<> > >::Type const)
{
    //typedef Index<TText, IndexSttd64<TIndexSpec> > TIndex;
    if (_isLeaf(value(it)))
    {
        return false;
    }
    typename VertexDescriptor<Index<TText, IndexSttd64<TIndexSpec> > >::Type vdFather = value(it);
    _getLeftMostChild(value(it), vdFather, container(it));

    return true;

//        if (_isLeaf(it, EmptyEdges()))
//            return false;
//        _historyPush(it);
//
//        if (value(it).isVirtual) {
//            return _goDownVirtual(it);
//        }
//
//        TIndex const &index = container(it);
//        Sttd64Node_ curNode;
//        Sttd64Node_ leftMostChildOfIt, leftMostChildOfChild;
//        if (!_getInternalNodeOfIterator(it, curNode))
//        {
//            //look for the first prefix-partition available
//            if (!_getPartitionNumber(it.vd.partitionNum, index, ""))
//                return false; //no tree at all
//
//            _getLeftMostChildNode(it, curNode); //get root of leftmost partition
//            //partition number already set
//            it.vd.node = 0;
//            it.vd.parentRepLen = curNode.lp;
//            it.vd.edgeLen = curNode.lp;
//            return true;
//        }
//        _getLeftMostChildNode(it, leftMostChildOfIt);
//        leftMostChildOfChild = getValue(getValue(index.trees, it.vd.partitionNum), leftMostChildOfIt.cleanPointerOrDepth());
//
//        it.vd.node = curNode.cleanPointerOrDepth();
//        it.vd.parentRepLen = it.vd.parentRepLen + it.vd.edgeLen;
//        it.vd.edgeLen = (leftMostChildOfChild.lp - leftMostChildOfIt.lp);
//        return true;
}

template <typename TIndex, typename TSize >
inline bool _isLeaf(VertexSttd64<TIndex, TSize> const & vDesc)
{
    //this is a leaf if
    //1) virtual and only one suffix in the preflen last suffices of the text (no partition child)
    //2) non-virtual and node is a leaf

    if (vDesc.inVirtual) {
        if (vDesc.isVirtualNode == false)
            return true;

        return false;

    }

    if (vDesc.node.isLeaf())
        return true;

    return false;
}

// is this a leaf? (including empty $-edges)
template <typename TText, typename TIndexSpec, class TSpec, typename TDfsOrder >
inline bool _isLeaf(
    Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TSpec> > const &it,
    VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _isLeaf(value(it));
}

// is this a leaf? (hide empty $-edges)
template <typename TText, typename TIndexSpec, class TSpec, typename TDfsOrder >
inline bool _isLeaf(
    Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TSpec> > const &it,
    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    if (_isLeaf(value(it))) return true;
    if (isRoot(it)) return false;

    typedef Index<TText, IndexSttd64<TIndexSpec> > TIndex;
    
    TIndex const& index = container(it);
    const typename VertexDescriptor<TIndex>::Type & vd = value(it);
    SEQAN_ASSERT_NOT(vd.inVirtual); //TODO(aumann): implement as if as soon as working on string sets
//        if (vd.inVirtual)
//        {
//            //TOD: adjust to new virtual tree concept
//            unsigned firstPartitionChild;
//            unsigned lastPartitionChild;
//            bool hasPartitionChildren;
//            bool hasChildren = index.virtualPrefixTree.children(firstPartitionChild,
//                                                                lastPartitionChild,
//                                                                hasPartitionChildren,
//                                                                std::vector<unsigned>(),
//                                                                vd.prefix,
//                                                                0);
//            SEQAN_ASSERT(hasChildren);
//            SEQAN_ASSERT_NOT(vd.preflenSuffixMatches.empty());  //otherwise should not be virtual
//
//            if (vd.leftMostPrefPartition == std::numeric_limits<unsigned>::max()
//                    && vd.prefLenSuffixMatches.size() < 2)
//                return true;
//
//            //if the tree is virtual and has at least one partition as child, there are at least
//            //2 children in sum
//            return false;
//        }

    Sttd64Node_ currentChild;
    typedef typename VertexDescriptor<TIndex>::Type TVertex;
    typename TVertex::TPosition currentChildPos;

    _getLeftMostChildNode(currentChild, currentChildPos, value(it), container(it));

    //this is only necessary when operating on string sets (i.e. in the future)
    //a single string can only produce nodes with at most one sentinel edge
    bool foundNotSentinelLeafChild = false;
    do {
        if (!currentChild.isLeaf() || currentChild.lp < length(indexText(index)))
        {
            foundNotSentinelLeafChild = true;
        }
        ++currentChildPos;
        unsigned internalFilePos;
        if (index.inMem)
            internalFilePos = vd.partitionNum;
        else
            internalFilePos = 0;
        currentChild = getValue(getValue(index.trees, internalFilePos), currentChildPos);
    } while (!currentChild.isRightMost());
    if (!foundNotSentinelLeafChild)
        return true;    //all children are sentinel leaves

    return false;
}

template <typename TText, class TIndexSpec, class TSpec, typename TValue >
inline bool
_getNodeByChar(
    Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TSpec> > const &it,
    TValue c,
    typename VertexDescriptor<Index<TText, IndexSttd64<TIndexSpec> > >::Type &childDesc)
{
    typedef Index<TText, IndexSttd64<TIndexSpec> > TIndex;
    if (_isLeaf(value(it))) return false;

    _getLeftMostChild(childDesc, value(it), container(it));
    if (childDesc.inVirtual)
    {
        if (getValue(childDesc.prefix, length(childDesc.prefix) - 1) == c)
            return true;
    }
    else
    {
        unsigned positionOnEdge;
        if (childDesc.nodePos == 0 && TIndex::PREFIX_LENGTH > 0)
            positionOnEdge = length(childDesc.prefix) - 1;
        else
            positionOnEdge = 0;
        if (getValue(indexText(container(it)), childDesc.node.lp + positionOnEdge) == c)
            return true;
    }
    Iter<TIndex, VSTree<TSpec> > subTreeIt(container(it), MinimalCtor());
    value(subTreeIt) = childDesc;
    SEQAN_ASSERT_EQ(value(subTreeIt).prefix, childDesc.prefix);    //check if assignment works
    TValue compC;
    SEQAN_ASSERT_EQ(value(subTreeIt).prefix, childDesc.prefix);    //check if assignment works
    if (!goRight(subTreeIt)) return false;

    unsigned positionOnEdge;    //if the iterator came from the virtual tree, there might be an edge offset
    do
    {
        if (value(subTreeIt).inVirtual)
        {
            compC = getValue(value(subTreeIt).prefix, length(value(subTreeIt).prefix) - 1);
        }
        else
        {
            if (value(subTreeIt).nodePos == 0 && TIndex::PREFIX_LENGTH > 0)
                positionOnEdge = length(value(subTreeIt).prefix) - 1;
            else
                positionOnEdge = 0;
            compC = getValue(indexText(container(subTreeIt)), value(subTreeIt).node.lp + positionOnEdge);
        }
        if (compC == c)
        {
            childDesc = value(subTreeIt);
            return true;
        }
        //TODO(aumann): implement correctly, uncomment line: if (compC > c) return false;    // all children after first child lexicographically sorted
    } while(goRight(subTreeIt));

    return false;
}

template <typename TText, class TIndexSpec, class TSpec, typename TDfsOrder, typename THideEmptyEdges >
inline bool _goRight(
    Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
    VSTreeIteratorTraits<TDfsOrder, THideEmptyEdges> const)
{
    typedef Index<TText, IndexSttd64<TIndexSpec> > TIndex;
    TIndex const & index = container(it);
    if (isRoot(it)) return false;

    if (value(it).inVirtual || (value(it).nodePos == 0 && TIndex::PREFIX_LENGTH > 0))
    {
        if (!value(it).inVirtual)
        {
            value(it).inVirtual = true; //set to virtual, as virtual iteration is attempted
            //is set to false in _getVertexOfPartitionFather if necessary

            //close the current partition's tree file
            SEQAN_ASSERT_NEQ(value(it).partitionNum, std::numeric_limits<unsigned>::max());
            _unloadTreeFile(container(it), value(it).partitionNum);
        }

        bool canGoRight = index.virtualPrefixTree
                .rightSibling(value(it).prefix,
                              value(it).isVirtualNode,
                              value(it).partitionNum);
        if (canGoRight && value(it).isVirtualNode == false
                && value(it).partitionNum != std::numeric_limits<unsigned>::max())
        {
            //points at a partition, need to "go down"
            _getVertexOfPartitionFather(value(it), value(it).partitionNum, index);
        }
        return canGoRight;
    }

    if (value(it).node.isRightMost()) return false;

    ++value(it).nodePos;
    typename Iterator<typename TIndex::TPartTree const>::Type treeIt;
    unsigned internalFilePos;
    if (index.inMem)
        internalFilePos = value(it).partitionNum;
    else
        internalFilePos = 0;
    treeIt = iter(getValue(index.trees, internalFilePos), value(it).nodePos);
    Sttd64Node_ rbNode = *treeIt;
    value(it).node = rbNode;
    //determine incoming edge length
    SEQAN_ASSERT_NOT(value(it).inVirtual);
    _setVertexEdgeLen(value(it), index);
    while ((THideEmptyEdges::VALUE && emptyParentEdge(it)) || !nodeHullPredicate(it))
    {
        if (value(it).node.isRightMost()) return false;
        ++value(it).nodePos;
        ++treeIt;
        value(it).node = *treeIt;
        _setVertexEdgeLen(value(it), index);
    }

    return true;
}

template <typename TText, class TIndexSpec, class TSpec>
inline void
_historyPush(Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TopDown<TSpec> > > &it)
{
    it._parentDesc = value(it);
}
template <typename TText, class TIndexSpec, class TSpec>
inline void
_historyPush(Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > &it)
{
    //TODO(aumann): implement for HistoryTopDown Iterator
//        typedef Iter<Index<TText, IndexEsa<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > TIter;
//        typename HistoryStackEntry_<TIter>::Type h;
//        h.range = value(it).range;
//
//        value(it).parentRight = value(it).range.i2;
//        appendValue(it.history, h);
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline typename Size<Index<TText, IndexSttd64<TIndexSpec> > >::Type
repLength(Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TopDown<TSpec> > > const &it)
{
    if (value(it).inVirtual)
    {
        return length(value(it).prefix);
    }
    return value(it).parentRepLen + value(it).edgeLen;
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline typename Size<Index<TText, IndexSttd64<TIndexSpec> > >::Type
parentRepLength(Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TopDown<TSpec> > > const &it)
{
    if (value(it).inVirtual)
    {
        return length(value(it).prefix) - 1;
    }
    return value(it).parentRepLen;
}

template <typename TText, typename TIndexSpec,class TSpec>
inline typename SAValue<Index<TText, IndexSttd64<TIndexSpec> > >::Type
getOccurrence(Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TSpec> > const &it)
{
    Index<TText, IndexSttd64<TIndexSpec> > const & index = container(it);
    if (_isLeaf(value(it)))
    {
        if (value(it).inVirtual)
        {
            //if a virtual leaf is encountered this can only be a "suffix-leaf"
            //otherwise the iterator would already point to a partition
            unsigned const prefixLength = Index<TText, IndexSttd64<TIndexSpec> >::PREFIX_LENGTH;
            std::vector<unsigned> sufficesMatching;
            sufficesMatching = index.virtualPrefixTree
                    .prefLenSufficesStartingWith(value(it).prefix,
                                                 sufficesMatching,
                                                 0);
            SEQAN_ASSERT_NOT(sufficesMatching.empty());

            return length(indexText(index)) - prefixLength + sufficesMatching[0];

        }
        return value(it).node.lp - value(it).node.cleanPointerOrDepth();    //leaf node
    }

    //find leftmost leaf (only if in virtual)
    typename VertexDescriptor<Index<TText, IndexSttd64<TIndexSpec> > >::Type leftMostChild;
    typename VertexDescriptor<Index<TText, IndexSttd64<TIndexSpec> > >::Type curFather;
    leftMostChild = value(it);
    //_getLeftMostChild(leftMostChild, value(it), container(it));
    if (leftMostChild.inVirtual)
    {
        while (!_isLeaf(leftMostChild))
        {
            curFather = leftMostChild;
            _getLeftMostChild(leftMostChild, curFather, container(it));
        }
    }

    if (leftMostChild.inVirtual)
    {
        SEQAN_ASSERT(_isLeaf(leftMostChild));
        //if a virtual leaf is encountered this can only be a "suffix-leaf"
        //otherwise the iterator would already point to a partition
        unsigned const prefixLength = Index<TText, IndexSttd64<TIndexSpec> >::PREFIX_LENGTH;
        std::vector<unsigned> sufficesMatching;
        sufficesMatching = index.virtualPrefixTree
                .prefLenSufficesStartingWith(value(it).prefix,
                                             sufficesMatching,
                                             0);
        SEQAN_ASSERT_NOT(sufficesMatching.empty());

        return length(indexText(index)) - prefixLength + 1 + sufficesMatching[0];
    }
    return leftMostChild.node.lp - leftMostChild.parentRepLen;
    //return leftMostChild.node.lp - leftMostChild.node.cleanPointerOrDepth();
}

template <typename TText, typename TIndexSpec>
inline void
_getLeftMostChild(typename VertexDescriptor<Index<TText, IndexSttd64<TIndexSpec> > >::Type & leftMostChild,
                  typename VertexDescriptor<Index<TText, IndexSttd64<TIndexSpec> > >::Type const & nodeVD,
                  Index<TText, IndexSttd64<TIndexSpec> > const & index)
{
    //typedef Index<TText, IndexSttd64<TIndexSpec> > TIndex;
    //typename ValueSize<typename Value<TText>::Type>::Type const alphabetSize = ValueSize<typename Value<typename TIndex::TText>::Type>::VALUE;
    SEQAN_ASSERT_NOT(_isLeaf(nodeVD));
    leftMostChild = nodeVD; //set partitionNum etc. correctly (same as father)
    SEQAN_ASSERT_EQ(leftMostChild.prefix, nodeVD.prefix);    //check if assignment works
    if (nodeVD.inVirtual)
    {
        index.virtualPrefixTree.leftMostChild(leftMostChild.prefix, leftMostChild.isVirtualNode, leftMostChild.partitionNum, nodeVD.prefix);
        if (leftMostChild.isVirtualNode == false
                && leftMostChild.partitionNum != std::numeric_limits<unsigned>::max())
        {
            //unsigned partitionNum;
            //_getPartitionNumber(partitionNum, leftMostChild.prefix, alphabetSize, index.PREFIX_LENGTH);
            _getVertexOfPartitionFather(leftMostChild, leftMostChild.partitionNum, index);
        }
        return;
    }   //node in virtual

    //as soon as we are down below the partition root, there is no more need for the vertex-
    //-stored prefix
    if (nodeVD.nodePos == 0)
    {
        clear(leftMostChild.prefix);
    }
    _getLeftMostChildNode(leftMostChild.node, leftMostChild.nodePos, nodeVD, index);

    //it remains to set the edgeLen and parentRepLen etc. right
    if (leftMostChild.node.isLeaf())
    {
        leftMostChild.edgeLen = length(indexText(index)) - leftMostChild.node.lp;//TODO(aumann): SHOULD BE WRONG:*/leftMostChild.node.lp - nodeVD.node.lp;
    }
    else
    {
        Sttd64Node_ leftMostChildOfLeftMostChild;
        typename SAValue<Index<TText, IndexSttd64<TIndexSpec> > >::Type dummyChildChildPos;
        _getLeftMostChildNode(leftMostChildOfLeftMostChild, dummyChildChildPos, leftMostChild, index);
        leftMostChild.edgeLen = leftMostChildOfLeftMostChild.lp - leftMostChild.node.lp;
    }
    leftMostChild.parentRepLen = nodeVD.parentRepLen + nodeVD.edgeLen;
    SEQAN_ASSERT_EQ(leftMostChild.partitionNum, nodeVD.partitionNum);
}

//sets the edge len of a vertex
//if it already points to the right (loaded) partition
//and contains the right Sttd64node + position
template <typename TText, typename TIndexSpec>
inline void
_setVertexEdgeLen(typename VertexDescriptor<Index<TText, IndexSttd64<TIndexSpec> > >::Type & nodeVD,
                  Index<TText, IndexSttd64<TIndexSpec> > const & index)
{
    //it remains to set the edgeLen and parentRepLen etc. right
    if (nodeVD.node.isLeaf())
    {
        nodeVD.edgeLen = length(indexText(index)) - nodeVD.node.lp;
    }
    else
    {
        Sttd64Node_ leftMostChildOfLeftMostChild;
        typename SAValue<Index<TText, IndexSttd64<TIndexSpec> > >::Type dummyChildChildPos;
        _getLeftMostChildNode(leftMostChildOfLeftMostChild, dummyChildChildPos, nodeVD, index);
        nodeVD.edgeLen = leftMostChildOfLeftMostChild.lp - nodeVD.node.lp;
    }
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline String<typename SAValue<Index<TText, IndexSttd64<TIndexSpec> > >::Type const >
getOccurrences(Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TSpec> > const &it)
{
    String<typename SAValue<Index<TText, IndexSttd64<TIndexSpec> > >::Type > occurrences;
    typedef Index<TText, IndexSttd64<TIndexSpec> > TIndex;

    if (value(it).inVirtual)
    {
        return _getOccurrencesFromVirtual(it);
    }
    if (value(it).node.isLeaf())
    {
        appendValue(occurrences, value(it).node.lp - value(it).node.cleanPointerOrDepth());
        return occurrences;
    }
    //simply move sequentially through the suffix tree and append all leaf values
    //first get firstchild

    typename Iterator<typename TIndex::TPartTree const>::Type nodeIt;
    unsigned internalTreeFilePosition;
    if (container(it).inMem)
        internalTreeFilePosition = value(it).partitionNum;
    else
        internalTreeFilePosition = 0;
    nodeIt = iter(value(container(it).trees, internalTreeFilePosition),
                     value(it).node.cleanPointerOrDepth());
    unsigned nodesEncountered = 0;
    unsigned rightMostsPassed = 0;
    Sttd64Node_ currentNode;
    do
    {
        currentNode = *nodeIt;
        if (currentNode.isLeaf())
        {
            appendValue(occurrences, currentNode.lp - currentNode.cleanPointerOrDepth());
        }
        else
        {
            ++nodesEncountered;
        }
        if (currentNode.isRightMost())
        {
            ++rightMostsPassed;
        }
        ++nodeIt;
    } while (rightMostsPassed <= nodesEncountered);
    return occurrences;
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline typename Size<Index<TText, IndexSttd64<TIndexSpec> > >::Type
countOccurrences(Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TSpec> > const &it)
{
    return length(getOccurrences(it));
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline String<typename SAValue<Index<TText, IndexSttd64<TIndexSpec> > >::Type const >
_getOccurrencesFromVirtual(Iter<Index<TText, IndexSttd64<TIndexSpec> >, VSTree<TSpec> > const &it)
{
    SEQAN_ASSERT(value(it).inVirtual);
    String<typename SAValue<Index<TText, IndexSttd64<TIndexSpec> > >::Type > occurrences;
    typedef Index<TText, IndexSttd64<TIndexSpec> > TIndex;
    //go through all possible child partitions
    unsigned leftMostPrefPartition, rightMostPrefPartition;
    TIndex const & index = container(it);
    unsigned partChildrenCount
            = index.virtualPrefixTree
            .partitionChildren(leftMostPrefPartition, rightMostPrefPartition, value(it).prefix,
                               std::numeric_limits<unsigned>::max());
    // Supress warning regarding unused variable
    ignoreUnusedVariableWarning(partChildrenCount);
    
    if (leftMostPrefPartition != std::numeric_limits<unsigned>::max())
    {
        SEQAN_ASSERT_NEQ(rightMostPrefPartition, std::numeric_limits<unsigned>::max());
        std::vector<bool>::const_iterator potPartsIt = container(it).virtualPrefixTree.prefixesContained.begin();
        potPartsIt += leftMostPrefPartition;
        SEQAN_ASSERT(*potPartsIt == true);
        for (unsigned i = leftMostPrefPartition;
                potPartsIt <= container(it).virtualPrefixTree.prefixesContained.begin() + rightMostPrefPartition;
                ++potPartsIt, ++i)
        {
            if (!*potPartsIt) continue;
            Iter<TIndex, VSTree<TSpec> > subTreeIt(container(it), MinimalCtor());
            typename VertexDescriptor<TIndex>::Type subTreeVD;
            _getVertexOfPartitionFather(subTreeVD, i, container(it));
            value(subTreeIt) = subTreeVD;
            append(occurrences, getOccurrences(subTreeIt));
            _unloadTreeFile(index, i);
        }
    }

    //remaining are the occurrences in the final
    unsigned const prefixLength = Index<TText, IndexSttd64<TIndexSpec> >::PREFIX_LENGTH;
    std::vector<unsigned> sufficesMatching;
    sufficesMatching = index.virtualPrefixTree
            .prefLenSufficesStartingWith(value(it).prefix,
                                         sufficesMatching,
                                         0);
    //SEQAN_ASSERT_NOT(sufficesMatching.empty()); not holding as is, see below
    SEQAN_ASSERT(!sufficesMatching.empty() || partChildrenCount >= 2);

    unsigned const textLen = length(indexText(index));

    for (unsigned i = 0; i < sufficesMatching.size(); i++)
    {
        appendValue(occurrences, textLen - prefixLength + 1 + sufficesMatching[i]);
    }
    return occurrences;
}

template <typename TText, typename TIndexSpec>
inline void
_getVertexOfPartitionFather(typename VertexDescriptor<Index<TText, IndexSttd64<TIndexSpec> > >::Type & vd,
                            unsigned partitionNumber,
                            Index<TText, IndexSttd64<TIndexSpec> > const & index)
{
    vd.partitionNum = partitionNumber;
    vd.inVirtual = false;
//        vd.leftMostPrefPartition = std::numeric_limits<unsigned>::max(); //invalid
//        vd.prefLenSuffixMatches.clear();
    vd.parentRepLen =  0;
    vd.edgeLen = 0;
    vd.nodePos = 0;
    if (!_reopenTreeFile(index, partitionNumber))
    {
        std::cerr << "Could not reopen tree file (Partition: " << partitionNumber
                << ")." << std::endl;
        SEQAN_ASSERT(false);
        return;
    }
    _getLeftMostChildNode(vd.node, vd.nodePos, vd, index);
    if (vd.node.isLeaf()) {
        vd.edgeLen = length(indexText(index)) - vd.node.lp;
    }
    else
    {
        Sttd64Node_ firstChildOfRoot;
        unsigned internalFilePos;
        if (index.inMem)
            internalFilePos = vd.partitionNum;
        else
            internalFilePos = 0;
        firstChildOfRoot = getValue(getValue(index.trees, internalFilePos), 1);
        vd.edgeLen = firstChildOfRoot.lp - vd.node.lp;
    }

}

template<typename TText, typename TIndexSpec>
inline void _getLeftMostChildNode(
        Sttd64Node_ & leftMostChild,
        typename VertexDescriptor<Index<TText, IndexSttd64<TIndexSpec> > >::Type::TPosition & leftMostChildPos,
        typename VertexDescriptor<Index<TText, IndexSttd64<TIndexSpec> > >::Type const & nodeVD,
        Index<TText, IndexSttd64<TIndexSpec> > const & index)
{
    typedef Index<TText, IndexSttd64<TIndexSpec> > TIndex;
    Sttd64Node_ iteratorsNode;
    SEQAN_ASSERT_NOT(nodeVD.inVirtual);
    SEQAN_ASSERT_NEQ(nodeVD.partitionNum, std::numeric_limits<unsigned>::max());
    unsigned internalTreeFilePos;
    if (index.inMem)
        internalTreeFilePos = nodeVD.partitionNum;
    else
        internalTreeFilePos = 0;
    if (nodeVD.parentRepLen == 0 && nodeVD.edgeLen == 0 && TIndex::PREFIX_LENGTH > 0)
    {
        //the "leftmost child" is the root of the current partition (partition number must already be set correctly!!)
        SEQAN_ASSERT_EQ(nodeVD.nodePos, 0u);
        leftMostChildPos = 0;
        leftMostChild = getValue(getValue(index.trees, internalTreeFilePos), leftMostChildPos);
    }
    else
    {
        leftMostChildPos = nodeVD.node.cleanPointerOrDepth();
        leftMostChild = getValue(getValue(index.trees, internalTreeFilePos), leftMostChildPos);
    }
}

} //namespace seqan

#endif // SEQAN_INDEX_INDEX_STTD64_STREE_H_
