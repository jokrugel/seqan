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
// Author: Johannes Krugel <krugel@in.tum.de>
// ==========================================================================
// Finder for IndexDigest.
// ==========================================================================

#ifndef SEQAN_INDEX_FIND_INDEX_DIGEST_H_
#define SEQAN_INDEX_FIND_INDEX_DIGEST_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template<typename TVerification = True, typename TSpec = void>
struct FinderDigest {};

// ----------------------------------------------------------------------------
// Finder
// ----------------------------------------------------------------------------

// TVerification = True | False
template <typename TText, typename TSpec, typename TVerification, typename TSpecFinder>
class Finder<Index<TText, IndexDigest<TSpec> >, FinderDigest<TVerification, TSpecFinder> > {
    typedef Index<TText, IndexDigest<TSpec> >               TIndex;
    typedef typename TIndex::TSize                          TSize;
    typedef typename TIndex::TDividers                      TDividers;
    typedef typename Size<TDividers>::Type                  TDividersSize;
    typedef typename TIndex::TBufSize                       TBufSize;
    
    typedef typename Position<TIndex>::Type                 TPosition;
    typedef Map<Pair<TText, String<TPosition> > >           TPreparedPositions;
    typedef Holder<TPreparedPositions>                      TPreparedPositionsHolder;
    typedef typename Iterator<TPreparedPositions, Rooted>::Type TPreparedPositionsIter;

public:
    TSize data_length; // The length of the match = the length of the needle.
    Holder<TIndex> index;
    bool first;
    bool wasPrefix;
    TDividersSize currentTree;
    TBufSize leftmostLeaf;
    TBufSize rightmostLeaf; // actually one behind the rightmostLeaf
    TBufSize currentLeaf;

    TPreparedPositionsHolder preparedPositions;
    TPreparedPositionsIter preparedPositionsIter;

    Finder() {
        clear(*this);
    }

    Finder(TIndex &_index): index(_index) {
        clear(*this);
        setHost(*this, _index);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TText, typename TIndexSpec>
struct DefaultFinder<Index<TText, IndexDigest<TIndexSpec> > > {
    typedef FinderDigest<True>                              Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TVerification, typename TSpecFinder>
inline void
clear(Finder<Index<TText, IndexDigest<TSpec> >, FinderDigest<TVerification, TSpecFinder> > & finder) {
    finder.data_length = 0u;
    finder.first = true;
    finder.wasPrefix = true;
    finder.currentTree = 0u;
    finder.leftmostLeaf = 0u;
    finder.rightmostLeaf = 0u;
    finder.currentLeaf = 0u;
//    clear(finder.preparedPositions);
}

// ----------------------------------------------------------------------------
// empty
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TVerification, typename TSpecFinder>
inline bool
empty(Finder<Index<TText, IndexDigest<TSpec> >, FinderDigest<TVerification, TSpecFinder> > const & finder) {
    return finder.leftmostLeaf == finder.rightmostLeaf;
}

// ----------------------------------------------------------------------------
// setHost
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TVerification, typename TSpecFinder>
inline void
setHost(Finder<Index<TText, IndexDigest<TSpec> >, FinderDigest<TVerification, TSpecFinder> > & finder, typename Parameter_<Index<TText, IndexDigest<TSpec> > >::Type & index) {
    //typedef Index<TText, IndexDigest<TSpec> > TIndex;

    clear(finder);
    // finder.index = Holder<TIndex>(index);
    setValue(finder.index, index);
}


// ----------------------------------------------------------------------------
// beginPosition
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TVerification, typename TSpecFinder>
inline typename Position<Finder<Index<TText, IndexDigest<TSpec> >, FinderDigest<TVerification, TSpecFinder> > >::Type
beginPosition(Finder<Index<TText, IndexDigest<TSpec> >, FinderDigest<TVerification, TSpecFinder> > & finder) {
    typedef Index<TText, IndexDigest<TSpec> >               TIndex;
    typedef typename TIndex::TSuffixTree                    TSuffixTree;
    typedef typename TIndex::TLeaves                        TLeaves;

    if (empty(value(finder.preparedPositions))) {
        TIndex const & index = haystack(finder);
        TSuffixTree const & tree = index.suffixTrees[finder.currentTree];
        TLeaves const & leaves = tree.i2;
        return leaves[finder.currentLeaf];
    } else {
        return value(finder.preparedPositionsIter).i2[finder.currentLeaf];
    }
}

// Uses finder.currentTree and sets finder.leftmostLeaf/finder.rightmostLeaf
template <typename TText, typename TSpec, typename TVerification, typename TSpecFinder, typename TNeedle, typename TIndex, typename TSpecPattern>
inline bool
findInCurrentTree(Finder<Index<TText, IndexDigest<TSpec> >, FinderDigest<TVerification, TSpecFinder> > & finder, 
        Pattern<TNeedle, PatternDigest<TIndex, TSpecPattern> > & pattern) {
    typedef typename Size<TNeedle>::Type                    TNeedleSize;
    typedef typename Value<TNeedle>::Type                   TNeedleValue;
    typedef typename TIndex::TSuffixTree                    TSuffixTree;
    typedef typename TIndex::TNodes                         TNodes;
    typedef typename TIndex::TNode                          TNode;
    typedef typename TIndex::TBufSize                       TBufSize;
    
    // Shortcuts
    TIndex & index = haystack(finder);
    TSuffixTree const & tree = index.suffixTrees[finder.currentTree];
    TNodes const & nodes = tree.i1;

    //std::cerr << "|nodes| = " << length(nodes) << std::endl;
    
    TBufSize nodeId = 0u; // index of root node
    finder.leftmostLeaf = 0u;
    finder.rightmostLeaf = length(index.suffixTrees[finder.currentTree].i2);
    TNeedleSize currentDepth = 0u;
    clear(pattern.itBit);
    
    static const unsigned BITS_PER_VALUE = BitsPerValue<TNeedleValue>::VALUE;
    TNeedleSize needleBitLength = BITS_PER_VALUE * length(needle(pattern));

    //for(unsigned int i = 0u; i < needleBitLength; ++i) {
    //    std::cerr << *pattern.itBit << ".";
    //    pattern.itBit += 1;
    //}
    //std::cerr << "#";
    //return false;
    
    // Descend in the tree using the blind skip trick (and verify afterwards)
    while (currentDepth < needleBitLength) {
        TNode const & node = nodes[nodeId];
        //finder.currentLeaf = finder.leftmostLeaf;
        //std::cerr << "a " << nodeId << ", Bit " << *pattern.itBit << " at depth " << currentDepth << " < " << needleBitLength << "      " << infixWithLength(indexText(index), beginPosition(finder), length(needle(pattern))) << " at position " << beginPosition(finder) << " with currentLeaf = " << finder.currentLeaf << ", repLength = " << nodes[nodeId].repLength << ", leftChild = " << node.leftChild << ", rightChild = " << node.rightChild << std::endl;
        if ((!((bool) *pattern.itBit)) && node.leftChild != TNode::NONE) {
            nodeId = node.leftChild;
            finder.leftmostLeaf = nodes[nodeId].leftmostLeaf; // This normally has no effect, but is needed for epsilon edges.
            if (node.rightChild != TNode::NONE) {
                finder.rightmostLeaf = nodes[node.rightChild].leftmostLeaf;
            }
        } else if (*pattern.itBit && node.rightChild != TNode::NONE) {
            nodeId = node.rightChild;
            finder.leftmostLeaf = nodes[node.rightChild].leftmostLeaf;
        } else {
            return false;
        }
        //std::cerr << "b " << nodeId << ", Bit " << *pattern.itBit << " at depth " << currentDepth << " < " << needleBitLength << "      " << infixWithLength(indexText(index), beginPosition(finder), length(needle(pattern))) << " at position " << beginPosition(finder) << " with currentLeaf = " << finder.currentLeaf << ", repLength = " << nodes[nodeId].repLength << ", leftChild = " << node.leftChild << ", rightChild = " << node.rightChild << std::endl;
        pattern.itBit += nodes[nodeId].repLength - currentDepth;
        currentDepth = nodes[nodeId].repLength;
    }
    //std::cerr << "depth: " << currentDepth << " < " << needleBitLength << std::endl;
    //std::cerr << "leftmostLeaf: " << finder.leftmostLeaf << std::endl;
    //std::cerr << "currentLeaf: " << finder.currentLeaf << std::endl;
    
    return true;
}

// ----------------------------------------------------------------------------
// preparePiecePositions
// ----------------------------------------------------------------------------

// Precondition: pieces are sorted lexicographically
template <typename TText, typename TSpec, typename TVerification, typename TSpecFinder, typename TNeedle, typename TIndex, typename TSpecPattern, typename TPreparedPositions>
inline void
preparePiecePositions(Finder<Index<TText, IndexDigest<TSpec> >, FinderDigest<TVerification, TSpecFinder> > & finder, Pattern<TNeedle, PatternDigest<TIndex, TSpecPattern> > & pattern, TPreparedPositions & preparedPositions) {
    
    //typedef Pattern<TNeedle, PatternDigest<TIndex, TSpecPattern> > TPattern;
    typedef typename TIndex::TSize                          TSize;
    typedef typename TIndex::TDividers                      TDividers;
    typedef typename Iterator<TDividers>::Type              TDividersIterator;
    typedef Lexical<LexicalSuffixDigest<True, TIndex::PREFIX_LENGTH, TSize, TNeedle, TText> > TLexical;
    typedef typename Iterator<TPreparedPositions, Rooted>::Type TPreparedPositionsIter;
    //typedef typename TIndex::TSuffixTree                    TSuffixTree;
    typedef typename TIndex::TLeaves                        TLeaves;
    typedef typename Infix<TText>::Type                     TTextInfix;

    TIndex & index = haystack(finder);
    indexRequire(index, DigestSuffixTrees());

    //for (TPreparedPositionsIter pieceIter = begin(preparedPositions); ! atEnd(pieceIter); ++pieceIter) {
    //    std::cerr << "piece = " << (*pieceIter).i1 << std::endl;
    //}
    //TPreparedPositionsIter pieceIter = begin(preparedPositions);
    //std::cerr << "piece = " << (*pieceIter).i1 << std::endl;
    //TText piece1 = (*pieceIter).i1;
    //for (unsigned i = 0; i < 374; i++) pieceIter++;
    //std::cerr << "piece = " << (*pieceIter).i1 << std::endl;
    //TText piece2 = (*pieceIter).i1;
    //std::cerr << "piece1 < piece2 = " << piece1[0] << " < " << piece2[0] << (piece1[0] < piece2[0]) << std::endl;

    TPreparedPositionsIter pieceIterLeft = begin(preparedPositions);

    // Simulatenously traverse over the sorted dividers = trees and the sorted pattern pieces.
    // We therefore store in pieceIterLeft the smallest possible piece of the last tree and start the iteration there.
    // The only tricky part is that some pieces can be contained in several trees.
    finder.currentTree = 0u;
    for (TDividersIterator it = iter(haystack(finder).dividers, finder.currentTree); !atEnd(it, haystack(finder).dividers); ++it, ++finder.currentTree) {
        //std::cerr << "Divider = " << prefix(suffix(indexText(haystack(finder)), (*it).i1), 20u) << std::endl;
        for (TPreparedPositionsIter pieceIter = pieceIterLeft; ! atEnd(pieceIter); ++pieceIter) {
            //std::cerr << "piece = " << (*pieceIter).i1 << std::endl;
            //- if pieceIter <= pieceIterRight: DO      Certainly contained in current tree
            // if lexSmaller <= divider: pieceIterLeft++ DO  
            // if isPrefix divider: DO    Might also be contained in next tree
            // if isGreater /*pieceIterRight = pieceIter*/: continue with next tree
            
            clear(pattern);
            TTextInfix sfx = suffix((*pieceIter).i1, 0u);
            setNeedle(pattern, sfx);
            TLexical lex(needle(pattern), pattern.digestSuffix, indexText(haystack(finder)), *it);
            //std::cerr << "Needle = " << needle(pattern) << ", " << pattern.digestSuffix.i1 << ", " << pattern.digestSuffix.i2[0u] << std::endl;

            if (isPrefix(lex) || isLess(lex)) {
                //std::cerr << "isPrefix = " << isPrefix(lex) << ", isLess = " << isLess(lex) << ", isGreater = " << isGreater(lex) << ", pos = " << pieceIterLeft << std::endl;
                // For the next tree we do not have to consider this piece again
                if (! isPrefix(lex)) pieceIterLeft++;

                // Load the tree from disk
                if(findInCurrentTree(finder, pattern)) {
                    TLeaves const & leaves = index.suffixTrees[finder.currentTree].i2;
                    // Verify for an actual match using the text
                    if ((*pieceIter).i1 != infixWithLength(indexText(index), leaves[finder.leftmostLeaf], length((*pieceIter).i1))) continue;
                    //std::cerr << (*pieceIter).i1 << " occurrs at positions (" << infix(leaves, finder.leftmostLeaf, finder.rightmostLeaf) << ")" << std::endl;

                    // Append to the position list of the current piece
                    append((*pieceIter).i2, infix(leaves, finder.leftmostLeaf, finder.rightmostLeaf));
                //} else {
                //    std::cerr << (*pieceIter).i1 << " does not occur." << std::endl;
                }
            } else {   // If current piece is greater than divider
                //std::cerr << "next tree" << std::endl;
                break; // Continue with next tree
            }
        }
    }
    clear(finder); // especially set rightmostLeaf = 0u
}

// ----------------------------------------------------------------------------
// setPreparedPiecePositions
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TVerification, typename TSpecFinder, typename TPreparedPositions>
inline void
setPreparedPiecePositions(Finder<Index<TText, IndexDigest<TSpec> >, FinderDigest<TVerification, TSpecFinder> > & finder, TPreparedPositions & preparedPositions) {
    setValue(finder.preparedPositions, preparedPositions);
}

// ----------------------------------------------------------------------------
// _goFirstMatchingTree
// ----------------------------------------------------------------------------

// Result: finder.currentTree is set to the first matching tree
template <typename TText, typename TSpec, typename TVerification, typename TSpecFinder, typename TNeedle, typename TIndex, typename TSpecPattern>
inline bool
_goFirstMatchingTree(Finder<Index<TText, IndexDigest<TSpec> >, FinderDigest<TVerification, TSpecFinder> > & finder,
        Pattern<TNeedle, PatternDigest<TIndex, TSpecPattern> > & pattern) {
    //typedef Pattern<TNeedle, PatternDigest<TIndex, TSpecPattern> > TPattern;
    typedef typename TIndex::TSize                          TSize;
    typedef typename TIndex::TDividers                      TDividers;
    typedef typename Iterator<TDividers, Rooted>::Type      TDividersIterator;
    typedef Lexical<LexicalSuffixDigest<True, TIndex::PREFIX_LENGTH, TSize, TNeedle, TText> > TLexical;

    // Binary search in dividers
    
    TDividersIterator leftDivider = begin(haystack(finder).dividers);
    TDividersIterator rightDivider = end(haystack(finder).dividers);
    TDividersIterator midDivider = begin(haystack(finder).dividers);
    //std::cerr << "init: " << position(leftDivider) << ", " << position(rightDivider) << std::endl;
    
    while (position(leftDivider) < position(rightDivider)) {
        setPosition(midDivider, (position(leftDivider) + position(rightDivider)) / 2);
        TLexical lex(needle(pattern), pattern.digestSuffix, indexText(haystack(finder)), *midDivider);

        //std::cerr << position(leftDivider) << ", " << position(rightDivider) << std::endl;
        //std::cerr << "Needle = " << needle(pattern) << ", " << pattern.digestSuffix.i1 << ", " << length(pattern.digestSuffix.i2) << std::endl;
        //std::cerr << "Divider = " << suffix(indexText(haystack(finder)), (*midDivider).i1) << std::endl;
        if (isPrefix(lex) || isLess(lex)) {
            rightDivider = midDivider;
        } else {
            leftDivider = midDivider + 1u;
        }
    }

    //std::cerr << "currentTree = " << position(leftDivider) << std::endl;
    finder.currentTree = position(leftDivider);
    return ! atEnd(leftDivider, haystack(finder).dividers);
}

// ----------------------------------------------------------------------------
// _goNextMatchingTree
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TVerification, typename TSpecFinder, typename TNeedle, typename TIndex, typename TSpecPattern>
inline bool
_goNextMatchingTree(Finder<Index<TText, IndexDigest<TSpec> >, FinderDigest<TVerification, TSpecFinder> > & finder,
        Pattern<TNeedle, PatternDigest<TIndex, TSpecPattern> > & pattern) {
    typedef typename TIndex::TSize                          TSize;
    typedef Lexical<LexicalSuffixDigest<True, TIndex::PREFIX_LENGTH, TSize, TNeedle, TText> > TLexical;

    ++finder.currentTree;
    if (finder.currentTree >= length(haystack(finder).dividers)) return false;

    TLexical lex(needle(pattern), pattern.digestSuffix, indexText(haystack(finder)), haystack(finder).dividers[finder.currentTree]);
    //std::cerr << "Needle = " << needle(pattern) << ", " << pattern.digestSuffix.i1 << ", " << length(pattern.digestSuffix.i2) << std::endl;
    //std::cerr << "Divider = " << suffix(indexText(haystack(finder)), (*it).i1) << std::endl;
    if (isPrefix(lex) || isLess(lex)) return true;
    else return false;
}

// ----------------------------------------------------------------------------
// find
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TVerification, typename TSpecFinder, typename TNeedle, typename TIndex, typename TSpecPattern>
inline bool
find(Finder<Index<TText, IndexDigest<TSpec> >, FinderDigest<TVerification, TSpecFinder> > & finder, 
        Pattern<TNeedle, PatternDigest<TIndex, TSpecPattern> > & pattern) {

    if (empty(value(finder.preparedPositions))) {
        return findRegular(finder, pattern);
    } else {
        // If prepared, simply lookup
        return findUsingPreparePositions(finder, pattern);
    }
}

template <typename TText, typename TSpec, typename TVerification, typename TSpecFinder, typename TNeedle, typename TIndex, typename TSpecPattern>
inline bool
findUsingPreparePositions(Finder<Index<TText, IndexDigest<TSpec> >, FinderDigest<TVerification, TSpecFinder> > & finder, 
        Pattern<TNeedle, PatternDigest<TIndex, TSpecPattern> > & pattern) {

    // We abuse the currentLeaf pointers here for the prepared positions
    
    // To instead use the regular find function when testing
    //return findRegular(finder, pattern);

    if (finder.rightmostLeaf == 0u) {
        finder.currentLeaf = 0u;
        finder.preparedPositionsIter = find(value(finder.preparedPositions), needle(pattern));
        //std::cerr << length(value(finder.preparedPositions)) << "; \"" << needle(pattern) << "\"" << hasKey(value(finder.preparedPositions), needle(pattern)) << ".";
        if (atEnd(finder.preparedPositionsIter)) { //, value(finder.preparedPositions))) {
        //if (hasKey(value(finder.preparedPositions), needle(pattern))) {        
            //std::cerr << "l = " << length((*finder.preparedPositionsIter).i2);
            finder.rightmostLeaf = 0u;
            return false;
        }

        finder.rightmostLeaf = length((*finder.preparedPositionsIter).i2);  // length of the positions array of this pattern
    } else {
        ++finder.currentLeaf;
    }
    // Do we still have occurrences in the current sequence of prepared positions?
    return finder.currentLeaf < finder.rightmostLeaf;
}

template <typename TText, typename TSpec, typename TVerification, typename TSpecFinder, typename TNeedle, typename TIndex, typename TSpecPattern>
inline bool
findRegular(Finder<Index<TText, IndexDigest<TSpec> >, FinderDigest<TVerification, TSpecFinder> > & finder, 
        Pattern<TNeedle, PatternDigest<TIndex, TSpecPattern> > & pattern) {

    // Shortcut
    TIndex & index = haystack(finder);
    
    // We still have occurrences in the current tree
    if (! finder.first && finder.currentLeaf + 1u < finder.rightmostLeaf) {
        ++finder.currentLeaf;
        return true;
    }

    if (finder.first) {
        finder.data_length = length(needle(pattern));
        indexRequire(index, DigestSuffixTrees());
        if (! _goFirstMatchingTree(finder, pattern)) return false;
        finder.first = false;
    } else {
        if (! _goNextMatchingTree(finder, pattern)) return false;
    }

    //std::cerr << "currentTree = " << finder.currentTree << " / " << length(haystack(finder).suffixTrees) << std::endl;
    //std::cerr << "pattern = " << needle(pattern) << std::endl;
    
    finder.currentLeaf = 0u;
    if (! findInCurrentTree(finder, pattern)) return false;

    // Verify for an actual match (has to be done only once for all leafs)
    finder.currentLeaf = finder.leftmostLeaf;
    
    //std::cerr << "Verifying positions from suffix array (" << finder.leftmostLeaf << ", " << finder.rightmostLeaf << ")... " << std::endl; // << suffixArray
    //std::cerr << "infixWithLength(" << "indexText(index)" << ", " << beginPosition(finder) << ", " << length(needle(pattern)) << ") = " << infixWithLength(indexText(index), beginPosition(finder), length(needle(pattern))) << std::endl;
    //std::cerr << "leaves[" << front(index.suffixTrees[finder.currentTree].i2) << " .. " << back(index.suffixTrees[finder.currentTree].i2) << "]" << std::endl;
    //std::cerr << "leaves[" << infixWithLength(indexText(index), front(index.suffixTrees[finder.currentTree].i2), 20) << " .. " << infixWithLength(indexText(index), back(index.suffixTrees[finder.currentTree].i2), 20) << "]" << std::endl;

    // If the caller told us not to verify we simply return the candidate match (which might be a false positive).
    if (! TVerification::VALUE) return true;
    
    // Verify
    if (needle(pattern) == infixWithLength(indexText(index), beginPosition(finder), length(needle(pattern)))) return true;

    //std::cerr << "VERIFY failed." << std::endl;

    return false;
}

}  // namespace seqan

#endif  // SEQAN_INDEX_FIND_INDEX_DIGEST_H_
