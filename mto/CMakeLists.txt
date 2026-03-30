# === CMake Setup ===
cmake_minimum_required(VERSION 3.14) # Minimum version required for FetchContent and modern features

# --- CMake Policy Handling ---
# Suppress warnings in CMake 3.24+ about DOWNLOAD_EXTRACT_TIMESTAMP
# See: https://cmake.org/cmake/help/latest/policy/CMP0135.html
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

# === Project Metadata ===
project(mto VERSION 0.0.1 LANGUAGES C CXX)

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
include(../cmake/sanitizers.cmake)

# --- Depds from CPM.cmake --
include(../cmake/deps.cmake)

# ============================ povu config =============================

# These variables define paths for source code, headers, and tests.
# Used later in the project to organize sources and targets.
set(APP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/app)
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(MTO_SOURCES_DIR ${SRC_DIR}/mto)

# When set uses uses build dir for output binaries
# set(BIN_DIR ${CMAKE_BINARY_DIR}/bin)

# -------- povu Library & executable ----------------
#
# Build the core 'povulib' library, composed of implementation files from various modules.
# This library contains the core algorithms, data structures, and utilities used by povu.
# The library is built from the source files located in the 'include' directory.

#  ........ Define the library target ............


ADD_LIBRARY(mto # io library
  ${MTO_SOURCES_DIR}/from_pvst.cpp
  ${MTO_SOURCES_DIR}/to_pvst.cpp
  ${MTO_SOURCES_DIR}/from_gfa.cpp
  ${MTO_SOURCES_DIR}/to_gfa.cpp
  ${MTO_SOURCES_DIR}/common.cpp
  ${MTO_SOURCES_DIR}/to_vcf.cpp
  ${MTO_SOURCES_DIR}/from_vcf.cpp
)

target_include_directories(ita
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
