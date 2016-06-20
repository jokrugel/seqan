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
    typedef String<char>                               TString;
    typedef String<char, MMap<> >                      TFileString;
    typedef EditDistanceScore                          TScore;
    typedef Index<TString, IndexEsa<> >                TIndex;
    typedef DefaultIndexPattern<TString, TIndex>::Type TIndexPattern;
    
    TFileString fileString;
    open(fileString, "D:/Arbeit/Presentation/2015-12 Verteidigung/Dissertation.txt", OPEN_RDONLY);
    TString txt = fileString;
    TIndex idx(txt);               // The index structure.
    TString ndl("Mutations");      // The needle we search for.
    int limit = -1;                // The search tolerance.

    // Search with Backtracking and output the matches.
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

    // Search with Partitioning and output the matches.
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

    return 0;
}
