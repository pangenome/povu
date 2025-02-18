cmake_minimum_required(VERSION 3.14)

# Project name and language
project(povu LANGUAGES C CXX)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check if CMAKE_BUILD_TYPE is set; if not, default to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build (Debug, Release)" FORCE)
endif()


# Set Linker flags
#set(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++")

# common flags for all build types
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic")
# Debug-specific flags
set(DEBUG_FLAGS " -Wpadded -Wunused -Wno-unknown-pragmas -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-private-field -Wno-unused-local-typedefs")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG ${DEBUG_FLAGS}") # let cmake.sanitizers set debug level to -O1
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DRELEASE")

# Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake 3.24
# https://cmake.org/cmake/help/latest/policy/CMP0135.html
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

# Sanitizers are only enabled in Debug mode to avoid performance overhead in Release builds.
include(sanitizers.cmake)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Directories
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)

# add_subdirectory(deps/libhandlegraph)
# add_subdirectory(deps/WFA2-lib)
add_subdirectory(deps/liteseq)

# deps
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/liteseq/include
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/args/
  #${CMAKE_CURRENT_SOURCE_DIR}/deps/libhandlegraph/src/include
  #${CMAKE_CURRENT_SOURCE_DIR}/deps/WFA2-lib
  #${CMAKE_CURRENT_SOURCE_DIR}/deps/WFA2-lib/system
  #${CMAKE_CURRENT_SOURCE_DIR}/deps/WFA2-lib/wavefront
  #${CMAKE_CURRENT_SOURCE_DIR}/deps/WFA2-lib/bindings/cpp
)

ADD_LIBRARY(LibsModule

  # src/io/vcf.cpp

  # genomics
  ${INCLUDE_DIR}/genomics/untangle.cpp
  ${INCLUDE_DIR}/genomics/variants.cpp
  ${INCLUDE_DIR}/genomics/genomics.cpp
  #src/genomics/vcf.cpp

  # graph
  ${INCLUDE_DIR}/graph/bidirected.cpp
  #${INCLUDE_DIR}/graph/biedged.cpp
  ${INCLUDE_DIR}/graph/flubble_tree.cpp
  #src/graph/graph.cpp
  ${INCLUDE_DIR}/graph/spanning_tree.cpp
  ${INCLUDE_DIR}/graph/bracket_list.cpp

  # common
  ${INCLUDE_DIR}/common/types.cpp
  ${INCLUDE_DIR}/common/utils.cpp

  # algorithms
  ${INCLUDE_DIR}/algorithms/algorithms.cpp

  # align
  # ${INCLUDE_DIR}/align/align.cpp
)


set(SOURCE_FILES
  # cli
  ${SRC_DIR}/cli/cli.cpp
  ${SRC_DIR}/cli/app.cpp

  # io
  ${SRC_DIR}/io/bub.cpp
  ${SRC_DIR}/io/from_gfa.cpp
  ${SRC_DIR}/io/common.cpp

  # subcommand
  ${SRC_DIR}/subcommand/common.cpp
  ${SRC_DIR}/subcommand/deconstruct.cpp
  ${SRC_DIR}/subcommand/info.cpp
  ${SRC_DIR}/subcommand/call.cpp

  ${SRC_DIR}/main.cpp
)

# source
add_executable(povu
  ${SOURCE_FILES}
)

# Link the library (LibsModule) and libhandlegraph to your executable
target_link_libraries(povu
  PRIVATE
  LibsModule
  #handlegraph_shared
  #wfa2cpp
  liteseqlib
)

# create a bin directory (bin) and copy the binary to the bin directory
set(BINARY_DIR ./bin)
file(MAKE_DIRECTORY ${BINARY_DIR}/)
add_custom_command(TARGET povu
		   POST_BUILD
		   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:povu> ../${BINARY_DIR}/)
