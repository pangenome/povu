cmake_minimum_required(VERSION 3.14)

# Project name and language
project(povu LANGUAGES C CXX)

option(BUILD_STATIC
  "Build a fully static executable (including libstdc++, libgcc, glibc, etc.)"
  OFF
)

# 2) If static was requested, force static libs & flags
if(BUILD_STATIC)
  # Make all add_library() calls STATIC by default
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
  # Only search for .a
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a)

  # Pass -static to the linker
  if(UNIX AND NOT APPLE)
    set(CMAKE_EXE_LINKER_FLAGS
      "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
  endif()
endif()

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


add_subdirectory(deps/liteseq)

# deps
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/liteseq/include
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/args/
)

ADD_LIBRARY(povulib

  # algorithms
  ${INCLUDE_DIR}/algorithms/slubbles.cpp
  ${INCLUDE_DIR}/algorithms/mubbles.cpp
  ${INCLUDE_DIR}/algorithms/flubbles.cpp

  # align
  ${INCLUDE_DIR}/align/align.cpp

  # common
  ${INCLUDE_DIR}/common/types.cpp
  ${INCLUDE_DIR}/common/utils.cpp
  ${INCLUDE_DIR}/common/tree_utils.cpp

  # genomics
  ${INCLUDE_DIR}/genomics/untangle.cpp
  ${INCLUDE_DIR}/genomics/variants.cpp
  ${INCLUDE_DIR}/genomics/genomics.cpp
  ${INCLUDE_DIR}/genomics/vcf.cpp

  # graph
  ${INCLUDE_DIR}/graph/bidirected.cpp
  ${INCLUDE_DIR}/graph/bracket_list.cpp
  ${INCLUDE_DIR}/graph/spanning_tree.cpp
)


set(SOURCE_FILES
  # cli
  ${SRC_DIR}/cli/cli.cpp
  ${SRC_DIR}/cli/app.cpp

  # io
  ${SRC_DIR}/io/pvst.cpp
  ${SRC_DIR}/io/from_gfa.cpp
  ${SRC_DIR}/io/common.cpp
  ${SRC_DIR}/io/to_vcf.cpp

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

target_link_libraries(povu PRIVATE povulib liteseqlib)

# create a bin directory (bin) and copy the binary to the bin directory
set(BINARY_DIR ./bin)
file(MAKE_DIRECTORY ${BINARY_DIR}/)
add_custom_command(TARGET povu
		   POST_BUILD
		   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:povu> ../${BINARY_DIR}/)
