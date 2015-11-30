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
// Author: Alexander Aumann <aumann@in.tum.de>
// ==========================================================================
// Tests for the SeqAn module index2.
// ==========================================================================

#include "test_index_suffix_trees.h"
#include "test_index_digest.h"
#include "test_index_sttd64.h"

using namespace seqan;

SEQAN_BEGIN_TESTSUITE(test_index_suffix_trees)
{
    // Call tests.
    SEQAN_CALL_TEST(test_index_digest_bits);
    SEQAN_CALL_TEST(test_index_digest_bits_dna);
    SEQAN_CALL_TEST(test_index_digest_lexical);
    SEQAN_CALL_TEST(test_index_digest_construct_find);
    SEQAN_CALL_TEST(test_index_digest_construct_find_2);
    //SEQAN_CALL_TEST(test_index_digest_construct_find_file);
    //SEQAN_CALL_TEST(test_index_digest_open_save);

    SEQAN_CALL_TEST(test_index_sttd64_basicConstruction1);
    SEQAN_CALL_TEST(test_index_sttd64_basicConstruction2);
    SEQAN_CALL_TEST(test_index_sttd64_basicConstruction3);
    SEQAN_CALL_TEST(test_index_sttd64_basicConstruction4);
    SEQAN_CALL_TEST(test_index_sttd64_printTree);
    //SEQAN_CALL_TEST(test_index_sttd64_readFile);
    SEQAN_CALL_TEST(test_index_sttd64_iterateTopDown1);
    SEQAN_CALL_TEST(test_index_sttd64_iterateTopDown2);
    SEQAN_CALL_TEST(test_index_sttd64_iterateTopDown3);
    SEQAN_CALL_TEST(test_index_sttd64_iterateTopDown4);
    SEQAN_CALL_TEST(test_index_sttd64_iterateTopDown5);
    SEQAN_CALL_TEST(test_index_sttd64_iterateTopDownPrefLen2);
    SEQAN_CALL_TEST(test_index_sttd64_iterateTopDownPrefLen3);
    //SEQAN_CALL_TEST(test_index_sttd64_iterateTopDownPrefLen2Files);
    //SEQAN_CALL_TEST(test_index_sttd64_iterateTopDownPrefLen0Files);
    SEQAN_CALL_TEST(test_index_sttd64_iterateTopDownLongInput);
}
SEQAN_END_TESTSUITE
