// ==========================================================================
//                             find_index_approx
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
/*
cd ~/patternMatching/SeqAn/
touch demos/find_index_approx.cpp
make --directory=build/Release/demos/
./build/Release/bin/demo_find_index_approx
*/
// ==========================================================================
// Demo for module find_index_approx
// (Approximate Pattern Matching with Index Structures)
// ==========================================================================
// Author: Johannes Krugel <krugel@in.tum.de>
// ==========================================================================

#include <iostream>

// SeqAn core modules
#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/sequence.h>

// New modules
#include <seqan/find_index_approx.h>
#include <seqan/index2.h>

using namespace seqan;

int main(int /*argc*/, char const ** /*argv*/) {

    // Define the type of string, similarity measure and index structure
    typedef String<char>                                        TString;
    typedef EditDistanceScore                                   TScore;
    typedef Index<TString, IndexEsa<> >                         TIndex;
//  typedef Index<TString, IndexWotd<> >                        TIndex;
//  typedef Index<TString, IndexSttd64<> >                      TIndex;
//  typedef Index<TString, IndexDigest<> >                      TIndex;
//  typedef Index<TString, FMIndex<> >                          TIndex;
//  typedef Index<TString, IndexSadakane<> >                    TIndex;
//  typedef Index<TString, IndexLZ<> >                          TIndex;
//  typedef Index<TString, IndexQGram<UngappedShape<3> > >      TIndex;
//  typedef Index<TString, IndexQGram2L<UngappedShape<3> > >    TIndex;
    typedef DefaultIndexPattern<TString, TIndex>::Type TIndexPattern;
    
    TString txt("1234567890abcdefghijklXYZ1234567890abcdefghijklXYZ1234567890aPc"
        "defghijklXYZ1234567890");      // The haystack text in which we search.
    TIndex idx(txt);                    // The index structure.
    TString ndl("aPcQeRgSijkl");        // The needle we search for.
    int limit = -4;                     // The search tolerance.

    // Search with an online algorithm (without the index) and output the matches.
    {
        std::cout << std::endl << "Myers (online):" << std::endl;
        Finder<TString> fdr(txt);
        Pattern<TString, Myers<> > ptn(ndl, limit);
        while (find(fdr, ptn)) {
            while (findBegin(fdr, ptn)) {
                std::cout << "Found '" << infix(fdr) << "' with a score of ";
                std::cout << getBeginScore(ptn) << " (pos " << beginPosition(fdr);
                std::cout << " to " << endPosition(fdr) << ")." << std::endl;
            }
        }
    }
    // Search with all four index-based algorithms and output the matches.
    {
        std::cout << std::endl << "Backtracking:" << std::endl;
        Finder<TIndex, DPBacktracking<TScore> > fdr(idx);
        Pattern<TString, DPBacktracking<TScore> > ptn(ndl, limit);
        while (find(fdr, ptn)) {
            while (findBegin(fdr, ptn)) {
                std::cout << "Found '" << infix(fdr) << "' with a score of ";
                std::cout << getBeginScore(ptn) << " (pos " << beginPosition(fdr);
                std::cout << " to " << endPosition(fdr) << ")." << std::endl;
            }
        }
    }
    {
        std::cout << std::endl << "Partitioning into exact search:" << std::endl;
        Finder<TIndex, FinderPartitioning<> > fdr(idx);
        Pattern<TString, Partitioning<IntoExactSearch, TScore, TIndexPattern> > 
            ptn(ndl, limit);
        while (find(fdr, ptn)) {
            while (findBegin(fdr, ptn)) {
                std::cout << "Found '" << infix(fdr) << "' with a score of ";
                std::cout << getBeginScore(ptn) << " (pos " << beginPosition(fdr);
                std::cout << " to " << endPosition(fdr) << ")." << std::endl;
            }
        }
    }
    {
        std::cout << std::endl << "Intermediate partitioning:" << std::endl;
        Finder<TIndex, FinderPartitioning<DPBacktracking<> > > fdr(idx);
        Pattern<TString, Partitioning<Intermediate, TScore, DPBacktracking<> > >
            ptn(ndl, limit);
        // We explicitely set the number of pieces
        setNumberOfPieces(ptn, 2u);
        while (find(fdr, ptn)) {
            while (findBegin(fdr, ptn)) {
                std::cout << "Found '" << infix(fdr) << "' with a score of ";
                std::cout << getBeginScore(ptn) << " (pos " << beginPosition(fdr);
                std::cout << " to " << endPosition(fdr) << ")." << std::endl;
            }
        }
    }
    {
        std::cout << std::endl << "Partitioning with hierarchical verification:"
            << std::endl;
        Finder<TIndex, FinderPartitioning<FinderPartitioning<> > > fdr(idx);
        Pattern<TString, Partitioning<PartitioningHierarchical, TScore,
          Partitioning<IntoExactSearch, TScore, TIndexPattern> > > ptn(ndl, limit);
        // We explicitely set the number of pieces
        setNumberOfPieces(ptn, 2u);
        while (find(fdr, ptn)) {
            while (findBegin(fdr, ptn)) {
                std::cout << "Found '" << infix(fdr) << "' with a score of ";
                std::cout << getBeginScore(ptn) << " (pos " << beginPosition(fdr);
                std::cout << " to " << endPosition(fdr) << ")." << std::endl;
            }
        }
    }

    return 0;
}
