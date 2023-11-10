cmake_minimum_required(VERSION 3.14)
#set(CMAKE_C_COMPILER "gcc")
#set(CMAKE_CXX_COMPILER "g++")

project(povu)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
enable_testing()


if (MSVC)
    # warning level 4
    add_compile_options(/W4)
else()
    # additional warnings
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Include libhandlegraph library
add_subdirectory(
# deps/gfakluge
  deps/libhandlegraph
)

ADD_LIBRARY(LibsModule
  # io
  src/io/gfa_to_digraph.cpp

  # graph
  src/graph/tree.cpp
  src/graph/u_graph.cpp
  src/graph/digraph.cpp
  src/graph/bidirected.cpp
  src/graph/spanning_tree.cpp
  # 
  #src/vst/vst.cpp
  #src/vst/pst.cpp
  src/pvst/pvst.cpp
  src/pst/pst.cpp
  # genomics
  src/genomics/variant_calling.cpp
  src/genomics/vcf.cpp

  # algorithms
  src/algorithms/cycle_equiv.cpp

  # cli
  src/cli/cli.cpp
)

#include_directories(
#    ./deps/gfakluge/src
#    ./deps/gfakluge/src/tinyFA
#    ./deps/libhandlegraph/src/include
#)

include_directories(
  # deps
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/args/
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/gfakluge/src
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/gfakluge/src/tinyFA
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/libhandlegraph/src/include
)

set(SOURCE_FILES
  src/main.cpp
)



# source
add_executable(povu
  ${SOURCE_FILES}
)

# recommended but does not work
# Specify the include directories for the target
#target_include_directories(povu
#    PRIVATE
#    ${CMAKE_CURRENT_SOURCE_DIR}/deps/gfakluge/src
#    ${CMAKE_CURRENT_SOURCE_DIR}/deps/gfakluge/src/tinyFA
#    ${CMAKE_CURRENT_SOURCE_DIR}/deps/libhandlegraph/src/include
#)

#target_include_directories(povu PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/build/libs/libhandlegraph)


#
# Link the library (LibsModule) and libhandlegraph to your executable
#target_link_libraries(
#  povu
#  LibsModule
#  handlegraph_static
#)

# Link the library (LibsModule) and libhandlegraph to your executable
target_link_libraries(povu PRIVATE LibsModule handlegraph_shared)

set(BINARY_DIR ./bin)

file(MAKE_DIRECTORY ${BINARY_DIR}/)


add_custom_command(TARGET povu
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:povu> ../${BINARY_DIR}/)

add_executable( povu_test
  tests/main_tests.cc
)


#target_link_libraries(
#  povu_test
#  LibsModule
#  handlegraph
#  GTest::gtest_main
#)

# Link the library (LibsModule), libhandlegraph, and Google Test to your test executable
target_link_libraries(povu_test PRIVATE LibsModule handlegraph_static GTest::gtest_main)


include(GoogleTest)
gtest_discover_tests(povu_test)
