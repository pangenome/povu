# ==============================================================================
# CMAKE SETUP & PROJECT METADATA
# ==============================================================================
cmake_minimum_required(VERSION 3.14)

# Handle policy for FetchContent (CMake 3.24+)
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

project(povu VERSION 0.0.1 LANGUAGES C CXX)

# Standard C++17 configuration
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Generate compile_commands.json for LSP (clangd, etc.)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ==============================================================================
# BUILD FEATURES & OPTIONS
# ==============================================================================

# CUDA Support
option(ENABLE_CUDA "Enable CUDA acceleration" ON)
if(ENABLE_CUDA)
  message(STATUS "CUDA support enabled")
  set(MEZA_ENABLE_CUDA ON)
else()
  message(STATUS "CUDA support disabled")
endif()

# Testing Toggles
option(ENABLE_ALL_TESTS "Enable all tests (including dependencies)" OFF)
if(ENABLE_ALL_TESTS)
  message(STATUS "Testing internal deps")
  set(TEST_MEZA ON)
  set(TEST_OZA ON)
endif()

# Default Build Type (Release if not specified)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
endif()
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

# ==============================================================================
# LIBRARY DEFAULTS & LINKING STRATEGY
# ==============================================================================

# 1. THE "STATIC BY DEFAULT" RULE
# By default, 'add_library(name sources...)' creates STATIC libraries (.a).
# This is because the global variable BUILD_SHARED_LIBS defaults to OFF.
# Static linking copies the library code into the target, making the final
# binary more portable but larger. To change this default project-wide,
# one would set(BUILD_SHARED_LIBS ON).

# 2. POSITION INDEPENDENT CODE (PIC)
# Required for monorepos when linking static libs into shared libs (.so).
# In a monorepo, you often link Static libraries into Shared libraries
# (e.g., linking static 'quilt' into shared 'oza').
# On Linux, Shared libraries (.so) REQUIRE their dependencies to be
# Position Independent. If this is OFF, you get "relocation" linker errors.
# We enable this globally to ensure compatibility regardless of library type.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# 2. FULLY STATIC BINARY OPTION
# Attempts to bundle system runtimes (glibc, libstdc++) for maximum portability.
option(BUILD_STATIC "Build a fully static executable" OFF)

if(BUILD_STATIC)
    message(STATUS "Build Mode: FULLY STATIC (Experimental/Portable)")

    # Defaults all add_library calls to STATIC
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

    # Force search for .a instead of .so
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)

    # Linux-specific static flags
    if(UNIX AND NOT APPLE)
        set(CMAKE_EXE_LINKER_FLAGS
            "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
    endif()
else()
    message(STATUS "Build Mode: DYNAMIC (Standard/Shared)")
endif()

if(CI OR DEFINED ENV{CI})
  message(STATUS "CI mode enabled")
  add_compile_options(-Werror)
endif()


# ==============================================================================
# COMPILER WARNINGS & HARDENING
# ==============================================================================

set(COMMON_WARNINGS -Wall)

set(DEBUG_WARNINGS
  -Wextra -Wpedantic -Wshadow -Wunused -Wunknown-pragmas
  -Wunused-parameter -Wunused-variable -Wunused-function
  -Wunused-but-set-variable -Wmissing-field-initializers
  -Wunused-local-typedefs -Wuninitialized
)

# GCC-specific warnings
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  list(APPEND DEBUG_WARNINGS -Wmaybe-uninitialized)
endif()

# Clang/macOS specific suppressions
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  list(APPEND DEBUG_WARNINGS -Wno-uninitialized -Wno-gnu-zero-variadic-macro-arguments)
endif()

# CI Mode: Stricter rules
option(CI "Enable stricter build rules for CI" OFF)
if(CI OR DEFINED ENV{CI})
  message(STATUS "CI mode enabled")
endif()

# ==============================================================================
# EXTERNAL DEPENDENCIES
# ==============================================================================

# Sanitizers for Debugging
include(cmake/sanitizers.cmake)

# External package management (CPM.cmake)
include(cmake/deps.cmake)

# ==============================================================================
# PROJECT SUBDIRECTORIES (Internal Libraries)
# ==============================================================================

add_subdirectory(quilt)
add_subdirectory(log)
add_subdirectory(liteseq)

# Helper function to silence noise from dependencies
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

# Internal sub-modules
add_subdirectory(dynamo)
add_subdirectory(oza)
add_subdirectory(ita)
add_subdirectory(mto)
add_subdirectory(zien)
add_subdirectory(povu)
add_subdirectory(meza)

# ==============================================================================
# COMPILER FLAG APPLICATION
# ==============================================================================

# Helper function to apply warnings/optimizations based on build type
function(add_compile_options_for_build_type target)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(${target} PRIVATE ${COMMON_WARNINGS} ${DEBUG_WARNINGS} -g -Wno-system-headers)
    target_compile_definitions(${target} PRIVATE DEBUG)
  elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(${target} PRIVATE ${COMMON_WARNINGS} -O3 -Wno-system-headers)
    target_compile_definitions(${target} PRIVATE RELEASE)
  endif()
endfunction()

add_compile_options_for_build_type(povu)
add_compile_options_for_build_type(oza)

# ==============================================================================
# TARGET LINKING & DEPENDENCY GRAPH
# ==============================================================================

target_link_libraries(liteseq PRIVATE log)

target_link_libraries(quilt PUBLIC log fmt::fmt)

target_link_libraries(meza PUBLIC log quilt)

target_link_libraries(oza PRIVATE fmt::fmt liteseq quilt meza log)

target_link_libraries(ita PRIVATE log oza fmt::fmt liteseq meza dynamo quilt)

target_link_libraries(mto PRIVATE log oza fmt::fmt liteseq ita meza)

target_link_libraries(zien PRIVATE log oza ita liteseq fmt::fmt ncurses meza)

target_link_libraries(povu PRIVATE log oza zien ita mto liteseq fmt::fmt taywee::args ncurses meza)

# ==============================================================================
# OUTPUT & TESTING CONFIGURATION
# ==============================================================================

set(BIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

set_target_properties(povu PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})
message(STATUS "Build-time binary output directory: ${BIN_DIR}")

include(CTest)
