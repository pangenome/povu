# === CMake Setup ===
cmake_minimum_required(VERSION 3.14) # Minimum version required for FetchContent and modern features

# --- CMake Policy Handling ---
# Suppress warnings in CMake 3.24+ about DOWNLOAD_EXTRACT_TIMESTAMP
# See: https://cmake.org/cmake/help/latest/policy/CMP0135.html
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

# === Project Metadata ===
project(povu LANGUAGES C CXX)

# === C++ Standard Configuration ===
# Require C++17 for all targets
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# === Build Type Defaults ===
# Default to Release if user hasn't set a build type
# Can be overridden via -DCMAKE_BUILD_TYPE=Debug
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build (Debug, Release)" FORCE)
endif()


# === Compile Commands for Tooling ===
# Generate compile_commands.json for use with clangd, clang-tidy, etc.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# === Platform & Toolchain Configuration ===
#
# Optional static build toggle.
# Enable with: -DBUILD_STATIC=ON
# This forces fully static linking for portable binaries (primarily Linux).
option(BUILD_STATIC
  "Build a fully static executable (including libstdc++, libgcc, glibc, etc.)"
  OFF
)

if(BUILD_STATIC)
  # Use static libraries where possible
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

  # Linker should prefer .a files
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a)

  # Pass static linking flags to the linker (Linux only)
  if(UNIX AND NOT APPLE)
    set(CMAKE_EXE_LINKER_FLAGS
      "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
  endif()
endif()

# === Build Configuration ===
#
# This section configures compile-time warnings, optimization levels,
# and debug/release-specific behavior using generator expressions.
# Flags are applied conditionally based on the active build configuration.
# Users can add extra flags via CMake options if needed:
#   cmake -DCMAKE_C_FLAGS="-Wpadded" -DCMAKE_CXX_FLAGS="-Wpadded" ..

# --- Build Flags Configuration ---

# Warnings that apply to all build types
set(COMMON_WARNINGS
  -Wall
)

# Additional warnings enabled only for Debug builds
set(DEBUG_WARNINGS
  -Wextra
  -Wpedantic
  -Wshadow
  -Wunused
  -Wunknown-pragmas
  -Wunused-parameter
  -Wunused-variable
  -Wunused-function
  -Wunused-but-set-variable
  -Wmissing-field-initializers
  -Wunused-local-typedefs
  -Wuninitialized
  -Wmaybe-uninitialized
  # -Wpadded # Uncomment or gate by compiler if needed
)

# Apply compiler flags per build type
# Define macros for use in #ifdef DEBUG / #ifdef RELEASE
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(-g)
  add_compile_options(${COMMON_WARNINGS})
  add_compile_options(${DEBUG_WARNINGS})
  add_compile_definitions(DEBUG)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_compile_options(${COMMON_WARNINGS})
  add_compile_options(-O3)
  add_compile_definitions(RELEASE)
endif()


# --- CI Mode Configuration ---
# CI is a build environment, not a build type.
# Enable stricter checks (e.g., -Werror) using -DCI=ON or if ENV{CI} is set
option(CI "Enable stricter build rules for CI (e.g., -Werror)" OFF)

if(CI OR DEFINED ENV{CI})
  message(STATUS "CI mode enabled")
  add_compile_options(-Werror)
endif()


# === Enable Sanitizers for Debug Builds ===
#
# This file typically enables AddressSanitizer, UndefinedBehaviorSanitizer, etc.
# It's included unconditionally, but should internally check for Debug mode and compiler support.
# You can make it optional with 'include(... OPTIONAL)' if needed.
include(sanitizers.cmake)

# === Project Directory Structure ===
#
# These variables define paths for source code, headers, and tests.
# Used later in the project to organize sources and targets.
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)

# === povulib Library ===
#
# Build the core 'povulib' library, composed of implementation files from various modules.
# This library contains the core algorithms, data structures, and utilities used by povu.
# The library is built from the source files located in the 'include' directory.

ADD_LIBRARY(povulib

  # algorithms
  ${INCLUDE_DIR}/algorithms/concealed.cpp
  ${INCLUDE_DIR}/algorithms/midi.cpp
  ${INCLUDE_DIR}/algorithms/smothered.cpp
  ${INCLUDE_DIR}/algorithms/parallel.cpp
  ${INCLUDE_DIR}/algorithms/tiny.cpp
  ${INCLUDE_DIR}/algorithms/flubbles.cpp

  # align
  ${INCLUDE_DIR}/align/align.cpp

  # common
  ${INCLUDE_DIR}/common/types/graph.cpp
  ${INCLUDE_DIR}/common/types/pvst.cpp
  ${INCLUDE_DIR}/common/utils.cpp
  ${INCLUDE_DIR}/common/tree_utils.cpp


  # genomics
  ${INCLUDE_DIR}/genomics/untangle.cpp
  ${INCLUDE_DIR}/genomics/variants.cpp
  ${INCLUDE_DIR}/genomics/vcf.cpp
  ${INCLUDE_DIR}/genomics/utils.cpp

  # graph
  ${INCLUDE_DIR}/graph/bidirected.cpp
  ${INCLUDE_DIR}/graph/bracket_list.cpp
  ${INCLUDE_DIR}/graph/spanning_tree.cpp
)

# === Command Line Executable Sources ===
# Main executable source files are listed here.

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
  ${SRC_DIR}/subcommand/gfa2vcf.cpp

  ${SRC_DIR}/main.cpp
)

# === Build Executable ===
# Define the main executable target
add_executable(povu ${SOURCE_FILES})


# === Dependencies ===
# Use FetchContent to fetch and build dependencies from GitHub
include(FetchContent)


# --- args (command-line parser) ---
# Disable building example and test binaries that are unnecessary for povu
set(ARGS_BUILD_UNITTESTS OFF CACHE BOOL "" FORCE)
set(ARGS_BUILD_EXAMPLE OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
  args
  GIT_REPOSITORY https://github.com/Taywee/args.git
  GIT_TAG        114200a9ad5fe06c8dea76e15d92325695cf3e34 # v6.4.7
)
FetchContent_MakeAvailable(args)

# --- liteseqlib (lightweight genomics data library) ---
FetchContent_Declare(
  liteseq
  GIT_REPOSITORY https://github.com/pangenome/povu.git
  GIT_TAG        b2abf40c1ad334019604dc085cb38cde772f76d2
)
FetchContent_MakeAvailable(liteseq)

# --- fmt (fast formatting library) ---
FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG        e69e5f977d458f2650bb346dadf2ad30c5320281 # v10.2.1
)
FetchContent_MakeAvailable(fmt)


# --- Silence compiler warnings for dependencies ---
function(silence_warnings target)
  if(TARGET ${target})
    get_target_property(type ${target} TYPE)
    if(type STREQUAL "INTERFACE_LIBRARY")
      target_compile_options(${target} INTERFACE -w)
    else()
      target_compile_options(${target} PRIVATE -w)
    endif()
  endif()
endfunction()

silence_warnings(fmt)
silence_warnings(liteseq)
silence_warnings(args)


# === Link Dependencies ===
# Link fmt to povulib (used in the library implementation)
target_link_libraries(povulib PRIVATE fmt::fmt)

# Link everything needed for the executable
target_link_libraries(povu
  PRIVATE
  povulib        # internal library
  liteseq     # gfa handling
  fmt::fmt       # formatting
  taywee::args   # argument parsing
)

# === Post-Build: Copy Binary to bin/ ===
# Create a bin/ directory and copy the compiled binary to it
set(BINARY_DIR ./bin)
file(MAKE_DIRECTORY ${BINARY_DIR}/)
add_custom_command(TARGET povu
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:povu> ../${BINARY_DIR}/
)
