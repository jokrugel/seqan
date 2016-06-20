// ==========================================================================
//                             find_index_approx
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
// Tests for the SeqAn module find_index_approx.
// ==========================================================================

#include "test_find_index_approx.h"

using namespace seqan;

SEQAN_BEGIN_TESTSUITE(test_find_index_approx)
{
    // Call tests.
    SEQAN_CALL_TEST(test_find_index_approx_dpbacktracking_esa);
    SEQAN_CALL_TEST(test_find_index_approx_dpbacktracking_wotd);

    SEQAN_CALL_TEST(test_find_index_approx_partitioning_intoexactsearch_esa);
    SEQAN_CALL_TEST(test_find_index_approx_partitioning_intoexactsearch_wotd);
    SEQAN_CALL_TEST(test_find_index_approx_partitioning_intoexactsearch_sttd64);
    SEQAN_CALL_TEST(test_find_index_approx_partitioning_intoexactsearch_digest);
    SEQAN_CALL_TEST(test_find_index_approx_partitioning_intoexactsearch_fm);
    SEQAN_CALL_TEST(test_find_index_approx_partitioning_intoexactsearch_sadakane);
    SEQAN_CALL_TEST(test_find_index_approx_partitioning_intoexactsearch_lz);
    SEQAN_CALL_TEST(test_find_index_approx_partitioning_intoexactsearch_qgram);
    // q-sample cannot work correctly because needles are smaller than Q + stepSize - 1
    // SEQAN_CALL_TEST(test_find_index_approx_partitioning_intoexactsearch_qsample);
    SEQAN_CALL_TEST(test_find_index_approx_partitioning_intoexactsearch_qgram2l);
    SEQAN_CALL_TEST(test_find_index_approx_partitioning_intoexactsearch_esa_myers);
    SEQAN_CALL_TEST(test_find_index_approx_partitioning_intoexactsearch_esa_dpsearch);

    SEQAN_CALL_TEST(test_find_index_approx_partitioning_intermediate_esa);
    SEQAN_CALL_TEST(test_find_index_approx_partitioning_intermediate_wotd);
    SEQAN_CALL_TEST(test_find_index_approx_partitioning_intermediate_esa_numberofpieces);

    // Hierarchical partitioning with PieceFinder = Partitioning<IntoExactSearch>
    SEQAN_CALL_TEST(test_find_index_approx_partitioning_hierarchical_esa);

    // TODO(krugel) Implement Iterator<TopDown<ParentLinks> >
    //SEQAN_CALL_TEST(test_find_index_approx_dpbacktracking_sttd64);
    //SEQAN_CALL_TEST(test_find_index_approx_partitioning_intermediate_sttd64);

    SEQAN_CALL_TEST(test_find_index_approx_partitioning_intoexactsearch_digest_prepare);

    SEQAN_CALL_TEST(test_find_index_on_segments);

    // TODO(krugel) Test copycon etc.
    SEQAN_CALL_TEST(test_pattern_copycon);
    SEQAN_CALL_TEST(test_pattern_assign);

#ifdef SEQAN_CXX11_STANDARD
    SEQAN_CALL_TEST(test_pattern_movecon);
    SEQAN_CALL_TEST(test_pattern_moveassign);
    SEQAN_CALL_TEST(test_pattern_set_host);
#endif  // SEQAN_CXX11_STANDARD
}
SEQAN_END_TESTSUITE
