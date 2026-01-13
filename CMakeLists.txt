# === CMake Setup ===
cmake_minimum_required(VERSION 3.14) # Minimum version required for FetchContent and modern features

# --- CMake Policy Handling ---
# Suppress warnings in CMake 3.24+ about DOWNLOAD_EXTRACT_TIMESTAMP
# See: https://cmake.org/cmake/help/latest/policy/CMP0135.html
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

# === Project Metadata ===
project(povu VERSION 0.0.1 LANGUAGES C CXX)

# === C++ Standard Configuration ===
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# === Compile Commands for Tooling ===
# Generate compile_commands.json for use with clangd, clang-tidy, etc.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# === Build Configuration Options ===

# Optional testing toggle.
# Enable with -DPOVU_ENABLE_TESTING=ON
# This will build and run tests using GoogleTest.
option(POVU_ENABLE_TESTING "Build and run tests" OFF)

# --- Build Type Defaults ---
# Default to Release if user hasn't set a build type
# Can be overridden via -DCMAKE_BUILD_TYPE=Debug
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build (Debug, Release)" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

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



# This section configures compile-time warnings, optimization levels,
# and debug/release-specific behavior using generator expressions.
# Flags are applied conditionally based on the active build configuration.
# Users can add extra flags via CMake options if needed:
#   cmake -DCMAKE_C_FLAGS="-Wpadded" -DCMAKE_CXX_FLAGS="-Wpadded" ..

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
  # -Wmaybe-uninitialized
  # -Wpadded # Uncomment or gate by compiler if needed
)

# Add GCC-specific warnings
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  list(APPEND DEBUG_WARNINGS
    -Wmaybe-uninitialized)
endif()

# Add Clang-specific suppressions (especially for macOS)
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  list(APPEND DEBUG_WARNINGS
    -Wno-uninitialized
    -Wno-gnu-zero-variadic-macro-arguments
  )
endif()

# --- CI Mode Configuration ---
# CI is a build environment, not a build type.
# Enable stricter checks (e.g., -Werror) using -DCI=ON or if ENV{CI} is set
option(CI "Enable stricter build rules for CI (e.g., -Werror)" OFF)

if(CI OR DEFINED ENV{CI})
  message(STATUS "CI mode enabled")
  add_compile_options(-Werror)
endif()


# ========================== Dependencies ===============================

# --- Enable Sanitizers for Debug Builds ---
#
# This file typically enables AddressSanitizer, UndefinedBehaviorSanitizer, etc.
# It's included unconditionally, but should internally check for Debug mode and compiler support.
# sanitizers.cmake checks if build mode is Debug and if not it returns early.
include(cmake/sanitizers.cmake)

# --- Depds from CPM.cmake --
include(cmake/deps.cmake)

# ============================ povu config =============================

# These variables define paths for source code, headers, and tests.
# Used later in the project to organize sources and targets.
set(APP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/app)
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(TESTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)
set(BIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(POVULIB_SOURCES_DIR ${SRC_DIR}/povu)
set(ZIENLIB_SOURCES_DIR ${SRC_DIR}/zien)

# When set uses uses build dir for output binaries
# set(BIN_DIR ${CMAKE_BINARY_DIR}/bin)

# -------- povu Library & executable ----------------
#
# Build the core 'povulib' library, composed of implementation files from various modules.
# This library contains the core algorithms, data structures, and utilities used by povu.
# The library is built from the source files located in the 'include' directory.

#  ........ Define the library target ............

ADD_LIBRARY(povulib
  # algorithms
  ${POVULIB_SOURCES_DIR}/algorithms/concealed.cpp
  ${POVULIB_SOURCES_DIR}/algorithms/midi.cpp
  ${POVULIB_SOURCES_DIR}/algorithms/smothered.cpp
  ${POVULIB_SOURCES_DIR}/algorithms/parallel.cpp
  ${POVULIB_SOURCES_DIR}/algorithms/tiny.cpp
  ${POVULIB_SOURCES_DIR}/algorithms/flubbles.cpp

  # align
  ${POVULIB_SOURCES_DIR}/align/align.cpp

  # common
  ${POVULIB_SOURCES_DIR}/common/utils.cpp

  #overlays
  ${POVULIB_SOURCES_DIR}/overlay/overlay.cpp
  ${POVULIB_SOURCES_DIR}/overlay/sne.cpp

  # variation
  ${POVULIB_SOURCES_DIR}/variation/rov.cpp
  ${POVULIB_SOURCES_DIR}/variation/color.cpp

  # genomics
  #${POVULIB_SOURCES_DIR}/genomics/allele.cpp
  ${POVULIB_SOURCES_DIR}/genomics/genomics.cpp
  ${POVULIB_SOURCES_DIR}/genomics/graph.cpp
  ${POVULIB_SOURCES_DIR}/genomics/untangle.cpp
  ${POVULIB_SOURCES_DIR}/genomics/vcf.cpp


  # trees
  ${POVULIB_SOURCES_DIR}/tree/slice_tree.cpp

  # graph
  ${POVULIB_SOURCES_DIR}/graph/types.cpp
  ${POVULIB_SOURCES_DIR}/graph/bidirected.cpp
  ${POVULIB_SOURCES_DIR}/graph/bracket_list.cpp
  ${POVULIB_SOURCES_DIR}/graph/spanning_tree.cpp
  ${POVULIB_SOURCES_DIR}/graph/tree_utils.cpp

  # io
  ${POVULIB_SOURCES_DIR}/io/from_pvst.cpp
  ${POVULIB_SOURCES_DIR}/io/to_pvst.cpp
  ${POVULIB_SOURCES_DIR}/io/from_gfa.cpp
  ${POVULIB_SOURCES_DIR}/io/to_gfa.cpp
  ${POVULIB_SOURCES_DIR}/io/common.cpp
  ${POVULIB_SOURCES_DIR}/io/to_vcf.cpp
  ${POVULIB_SOURCES_DIR}/io/from_vcf.cpp
)

ADD_LIBRARY(zien
  ${ZIENLIB_SOURCES_DIR}/tui/tui.cpp

  #validate
  ${ZIENLIB_SOURCES_DIR}/validate/validate.cpp
)

# ..... executable .......
#
# Main executable source files are listed here.
# Define the main executable target
add_executable(povu
  # cli
  ${APP_DIR}/cli/cli.cpp

  # subcommand
  ${APP_DIR}/subcommand/decompose.cpp
  ${APP_DIR}/subcommand/info.cpp
  ${APP_DIR}/subcommand/call.cpp
  ${APP_DIR}/subcommand/prune.cpp
  ${APP_DIR}/subcommand/gfa2vcf.cpp
  ${APP_DIR}/subcommand/vcf.cpp
  ${APP_DIR}/main.cpp
)

function(add_compile_options_for_build_type target)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(${target} PRIVATE
      ${COMMON_WARNINGS}
      ${DEBUG_WARNINGS}
      -g -Wno-system-headers)
    target_compile_definitions(${target} PRIVATE DEBUG)
  elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(${target} PRIVATE
      ${COMMON_WARNINGS}
      -O3 -Wno-system-headers)
    target_compile_definitions(${target} PRIVATE RELEASE)
  endif()
endfunction()

# Apply compile options only to your own targets
add_compile_options_for_build_type(povulib)
add_compile_options_for_build_type(povu)


target_include_directories(povulib
  PUBLIC
  ${INCLUDE_DIR}
  ${APP_DIR}
)

# ---- linking -----

target_link_libraries(povulib
  PRIVATE
  fmt::fmt # formatting
  liteseq  # gfa handling
)

target_link_libraries(zien
  PRIVATE
  povulib
  liteseq      # gfa handling
  fmt::fmt     # formatting
  ncurses
)

target_link_libraries(povu
  PRIVATE
  povulib      # internal library
  liteseq      # gfa handling
  fmt::fmt     # formatting
  taywee::args # argument parsing
  zien        # zien library
  ncurses
)

# =================== Set Build-Time Output Directory ========================

# Set the runtime output directory for the binary
set_target_properties(povu PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
)
message(STATUS "Build-time binary output directory: ${BIN_DIR}")

# ========================== Installation ====================================
# Installation configuration for install() actions (system-wide install)
# triggered by 'cmake --install .'

# Install configuration
include(GNUInstallDirs)

# Install the povulib library
install(TARGETS povulib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}  # Static Libraries
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}  # Shared Libraries
)

# Install the povu binary into the system bin directory
install(TARGETS povu
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}  # Executables
)

# Install header files (if applicable)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.h"
)

# Log where items are being installed
message(STATUS "Install include directory: ${CMAKE_INSTALL_INCLUDEDIR}")
message(STATUS "Install library directory: ${CMAKE_INSTALL_LIBDIR}")
message(STATUS "Install binary directory: ${CMAKE_INSTALL_BINDIR}")

# ================================ Tests =====================================

if (POVU_ENABLE_TESTING)
  # Add the global TESTING definition
  add_compile_definitions(TESTING)

  enable_testing()
  add_subdirectory(${TESTS_DIR})
endif()
