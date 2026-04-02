# === CMake Setup ===
# Minimum version required for FetchContent and modern features
cmake_minimum_required(VERSION 3.14)

# --- CMake Policy Handling ---
# Suppress warnings in CMake 3.24+ about DOWNLOAD_EXTRACT_TIMESTAMP
# See: https://cmake.org/cmake/help/latest/policy/CMP0135.html
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

option(ENABLE_CUDA "Enable CUDA acceleration (build CUDA targets)" ON)
if(ENABLE_CUDA)
  message(STATUS "CUDA support enabled")
  set(MEZA_ENABLE_CUDA ON)
else()
  message(STATUS "CUDA support disabled")
endif()
#option(MEZA_BUILD_EXAMPLE "Build examples for meza" OFF)

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
#option(POVU_ENABLE_TESTING "Build and run tests" OFF)
option(ENABLE_ALL_TESTS "Enable all tests (including dependencies)" OFF)
#option(ENABLE_TESTS "Enable CTest testing framework" ON)
if(ENABLE_ALL_TESTS)
  message(STATUS "Testing internal deps")

  set(TEST_MEZA ON)
  set(TEST_OZA ON)
endif()

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


set(BIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

add_subdirectory(quilt)
add_subdirectory(meza)

target_link_libraries(meza
  PUBLIC
  quilt)

add_subdirectory(log)
add_subdirectory(liteseq)

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

silence_warnings(log)
target_link_libraries(liteseq PRIVATE log)

add_subdirectory(dynamo)
add_subdirectory(oza)
add_subdirectory(ita)
add_subdirectory(mto)
add_subdirectory(zien)
add_subdirectory(povu)

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

add_compile_options_for_build_type(povu)
add_compile_options_for_build_type(oza)

target_link_libraries(quilt
  INTERFACE
  fmt::fmt
)

target_link_libraries(oza
  PRIVATE
  fmt::fmt # formatting
  liteseq  # gfa handling
  quilt
  meza
)

target_link_libraries(ita
  PRIVATE
  oza
  fmt::fmt     # formatting
  liteseq      # gfa handling
  meza
  dynamo
)

target_link_libraries(mto
  PRIVATE
  oza
  fmt::fmt     # formatting
  liteseq      # gfa handling
  ita
  meza
)

target_link_libraries(zien
  PRIVATE
  oza
  ita
  liteseq      # gfa handling
  fmt::fmt     # formatting
  ncurses
  meza
)

target_link_libraries(povu
  PRIVATE

  # internal libraries
  oza
  zien
  ita
  mto

  # external libraries
  liteseq      # gfa handling
  fmt::fmt     # formatting
  taywee::args # argument parsing
  ncurses
  meza
)


# # =================== Set Build-Time Output Directory ========================

# Set the runtime output directory for the binary
set_target_properties(povu PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
)
message(STATUS "Build-time binary output directory: ${BIN_DIR}")

# ============================ Testing Configuration ===========================
include(CTest)
