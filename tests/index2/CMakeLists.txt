# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/index2/CMakeLists.txt
#
# CMakeLists.txt file for the index2 module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_tests_index2 CXX)
message (STATUS "Configuring tests/index2")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your test.

add_executable (test_find_index
                test_find_index.cpp
                test_find_index.h)

add_executable (test_index_qgram_ext
                test_index_qgram_ext.cpp
                test_index_qgram_ext.h)

add_executable (test_index_suffix_trees
                test_index_suffix_trees.cpp
                test_index_suffix_trees.h
                test_index_digest.h
                test_index_sttd64.h)

add_executable (test_index_stree_iterators_2
                test_index_stree_iterators_2.cpp
                test_index_stree_iterators_2.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_find_index ${SEQAN_LIBRARIES})
target_link_libraries (test_index_qgram_ext ${SEQAN_LIBRARIES})
target_link_libraries (test_index_suffix_trees ${SEQAN_LIBRARIES})
target_link_libraries (test_index_stree_iterators_2 ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_find_index COMMAND $<TARGET_FILE:test_find_index_approx>)
add_test (NAME test_test_index_qgram_ext COMMAND $<TARGET_FILE:test_index_qgram_ext>)
add_test (NAME test_test_index_suffix_trees COMMAND $<TARGET_FILE:test_index_suffix_trees>)
add_test (NAME test_test_index_stree_iterators_2 COMMAND $<TARGET_FILE:test_index_stree_iterators_2>)
