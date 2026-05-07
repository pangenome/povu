# ================================= CMake Setup ==============================
# Minimum version required for FetchContent and modern features
cmake_minimum_required(VERSION 3.14)


# --- CMake Policy Handling ---
# Suppress warnings in CMake 3.24+ about DOWNLOAD_EXTRACT_TIMESTAMP
# See: https://cmake.org/cmake/help/latest/policy/CMP0135.html
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()


# ========================= Build Configuration Options =====================

option(MEZA_ENABLE_CUDA "Enable CUDA acceleration (build CUDA targets)" OFF)
if(MEZA_ENABLE_CUDA)
  message(STATUS "CUDA support enabled")
endif()

option(MEZA_BUILD_EXAMPLE "Build examples" OFF)
if (MEZA_BUILD_EXAMPLE)
  message(STATUS "Build examples.")
else()
  message(STATUS "Not building examples.")
endif()

# enable or disable testing. OFF by default
option(MEZA_ENABLE_TESTS "Build tests" OFF)
if (MEZA_ENABLE_TESTS)
  message(STATUS "Testing is enabled.")
else()
  message(STATUS "Testing is disabled.")
endif()


# ============================ Project Definition ===========================
project(meza VERSION 0.0.1 LANGUAGES CXX)


# === C++ Standard Configuration ===
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# === Compile Commands for Tooling ===
# Generate compile_commands.json for use with clangd, clang-tidy, etc.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


if(MEZA_ENABLE_CUDA)
  enable_language(CUDA)

  # Prefer the modern package when available
  find_package(CUDAToolkit REQUIRED)

  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

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
  # add_compile_options(-Werror)
endif()


# ========================== Dependencies ===============================

# --- Enable Sanitizers for Debug Builds ---
#
# This file typically enables AddressSanitizer, UndefinedBehaviorSanitizer, etc.
# It's included unconditionally, but should internally check for Debug mode and compiler support.
# sanitizers.cmake checks if build mode is Debug and if not it returns early.
#include(cmake/sanitizers.cmake)

# --- Depds from CPM.cmake --
#include(cmake/deps.cmake)

# ============================ config =============================

# These variables define paths for source code, headers, and tests.
# Used later in the project to organize sources and targets.

set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(TESTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)
set(BIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(MEZA_SOURCES_DIR ${SRC_DIR}/meza)
set(EXAMPLES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/examples)

#  ........ Define the library target ............

ADD_LIBRARY(meza
  ${MEZA_SOURCES_DIR}/ops/ops.cpp
  ${MEZA_SOURCES_DIR}/pool/hap_comp.cpp
  ${MEZA_SOURCES_DIR}/pool/pool_ops.cpp
  ${MEZA_SOURCES_DIR}/owned/matrix.cpp
)

set_target_properties(meza PROPERTIES
  POSITION_INDEPENDENT_CODE ON
)

target_compile_definitions(meza PUBLIC
  $<$<CONFIG:Debug>:MEZA_DEBUG=1>
  $<$<CONFIG:Release>:MEZA_RELEASE=1>
  $<$<CONFIG:RelWithDebInfo>:MEZA_RELWITHDEBINFO=1>
  $<$<CONFIG:MinSizeRel>:MEZA_MINSIZEREL=1>
)

target_compile_options(meza PRIVATE
  $<$<COMPILE_LANGUAGE:C>:${COMMON_WARNINGS}>
  $<$<COMPILE_LANGUAGE:CXX>:${COMMON_WARNINGS}>

  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:${DEBUG_WARNINGS}>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:${DEBUG_WARNINGS}>

  $<$<CONFIG:Debug>:-g>
  $<$<CONFIG:Release>:-O3>
)

target_include_directories(meza
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# ---- linking -----

if(MEZA_ENABLE_CUDA)
  target_sources(meza PRIVATE
    ${MEZA_SOURCES_DIR}/info.cu
    ${MEZA_SOURCES_DIR}/ops/ops.cu
  )

  target_link_libraries(meza
    PRIVATE
    CUDA::cudart
  )

  # TODO: make this PRIVATE
  target_include_directories(meza INTERFACE
    ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
  )

  target_compile_options(meza PRIVATE
    # 1. Wrap warnings for the host compiler
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-Wno-pedantic>

    # 2. Fix the __builtin_dynamic_object_size error (Fortify Source)
    $<$<COMPILE_LANGUAGE:CUDA>:-U_FORTIFY_SOURCE>
    $<$<COMPILE_LANGUAGE:CUDA>:-D_FORTIFY_SOURCE=0>

    # 3. Fix the nvcc fatal: Value '-g' is not defined for option 'Werror'
    # We ensure -g is only passed to CXX, and use -g -G for CUDA debug
    $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:Debug>>:-g -G>
  )

  set_source_files_properties(
    ${MEZA_SOURCES_DIR}/mat.cu
    ${MEZA_SOURCES_DIR}/info.cu
    ${MEZA_SOURCES_DIR}/ops/ops.cu
    PROPERTIES LANGUAGE CUDA
  )
endif()

target_compile_definitions(meza
  PUBLIC
  MEZA_USE_CUDA=$<BOOL:${MEZA_ENABLE_CUDA}>
)

# ============================= Tests =====================================

if(TEST_MEZA)
  message(STATUS "Tests enabled for meza")
  # Add the global TESTING definition
  add_compile_definitions(TESTING)

  enable_testing()
  add_subdirectory(${TESTS_DIR})
endif()

# === Example Binary ===
if(MEZA_BUILD_EXAMPLE)
  add_subdirectory(${EXAMPLES_DIR})
endif()
