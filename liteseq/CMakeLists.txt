# === CMake Setup ===
cmake_minimum_required(VERSION 3.14)

# Project name, version and languages
project(liteseq VERSION 0.0.1 LANGUAGES C CXX)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/liteseq/version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/generated/liteseq/version.h
  @ONLY)

# === BUILD OPTIONS ===
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
# Check if CMAKE_BUILD_TYPE is set; if not, default to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build (Debug, Release)" FORCE)
endif()

# enable or disable testing. OFF by default
option(LITESEQ_ENABLE_TESTING "Build tests" OFF)
if (LITESEQ_ENABLE_TESTING)
  message(STATUS "Testing is enabled.")
else()
  message(STATUS "Testing is disabled.")
endif()

option(LITESEQ_BUILD_EXAMPLE "Build examples" OFF)
if (LITESEQ_BUILD_EXAMPLE)
  message(STATUS "Build examples.")
else()
  message(STATUS "LQ not building examples disabled.")
endif()

# === Platform & Toolchain Configuration ===

# Directories
set(EXAMPLES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/examples)
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SRC_INTERNAL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/internal)
set(TESTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)
#set(BIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# --- liteseq Library ---

# Add the library
add_library(liteseq
  ${SRC_INTERNAL_DIR}/lq_utils.c
  ${SRC_INTERNAL_DIR}/lq_io.c
  ${SRC_DIR}/gfa.c
  ${SRC_DIR}/gfa_l.c
  ${SRC_DIR}/gfa_s.c
  ${SRC_DIR}/refs/ref_impl.c
  ${SRC_DIR}/refs/ref_walk.c
  ${SRC_DIR}/refs/ref_name.c
)

# C standard & basic props
set_target_properties(liteseq PROPERTIES
  C_STANDARD 17
  C_STANDARD_REQUIRED YES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED YES
  C_EXTENSIONS OFF
)

# Compiler options
# Debug-specific flags

set(COMMON_WARNINGS
  -Wall
)

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
  # -Wpadded # Uncomment or gate by compiler if needed
)

# compiler-specific warnings (check the CMAKE_C_COMPILER_ID)
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  message(STATUS "Configuring warnings for GCC...")
  list(APPEND DEBUG_WARNINGS
    -Wunused-but-set-variable
    -Wmaybe-uninitialized
  )
elseif (CMAKE_C_COMPILER_ID STREQUAL "Clang")
  message(STATUS "Configuring warnings for Clang...")
  # Add any Clang-specific warnings here (if necessary)
elseif (CMAKE_C_COMPILER_ID STREQUAL "MSVC")
  message(STATUS "Configuring warnings for MSVC... (none configured yet)")
  # Add MSVC-specific flags here if needed
endif()


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


# Linker options for sanitizers
# Linker options for sanitizers
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  # Default to ON in Debug mode, but allow manual override
  option(USE_SANITIZER "Enable sanitizers (ON by default in Debug mode)" ON)
else()
  # Default to OFF in non-Debug modes, and allow manual override
  option(USE_SANITIZER "Enable sanitizers (ON by default in Debug mode)" OFF)
endif()

# Include sanitizers only if USE_SANITIZER is ON
if (USE_SANITIZER)
  include(cmake/sanitizers.cmake)
endif()

# ============== Dependencies ===================
include(cmake/deps.cmake)

# # --- Silence compiler warnings for dependencies ---
# function(silence_warnings target)
#   if(TARGET ${target})
#     get_target_property(type ${target} TYPE)
#     if(type STREQUAL "INTERFACE_LIBRARY")
#       target_compile_options(${target} INTERFACE -w)
#     else()
#       target_compile_options(${target} PRIVATE -w)
#     endif()
#   endif()
# endfunction()

# silence_warnings(log)

# target_link_libraries(liteseq PRIVATE log)

# --- Math Library (libm) ---

message("Looking for math library (libm)...")
find_library(MATH_LIB m)
if(MATH_LIB)
  message(STATUS "Found math library: ${MATH_LIB}")
  # link math library libm
  target_link_libraries(liteseq PRIVATE m)
else()
  message(FATAL_ERROR "Math library (libm) not found.")
endif()

# --- Threading Configuration ---

message(STATUS "Looking for an appropriate threading library...")
if (APPLE)
  # Log the start of macOS SDK retrieval
  message(STATUS "Retrieving macOS SDK path using xcrun...")

  # Get the macOS SDK path dynamically
  execute_process(
    COMMAND xcrun --show-sdk-path
    OUTPUT_VARIABLE MACOSX_SDK_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Check if the macOS SDK path was successfully retrieved
  if (NOT MACOSX_SDK_PATH)
    message(FATAL_ERROR "Failed to retrieve macOS SDK path. Ensure that Xcode tools are installed and configured.")
  else()
    # Log the retrieved macOS SDK path
    message(STATUS "Using macOS SDK path: ${MACOSX_SDK_PATH}")

    # Add the SDK's /usr/include directory to the CMake include paths
    message(STATUS "Adding ${MACOSX_SDK_PATH}/usr/include to the include directories.")
    include_directories(${MACOSX_SDK_PATH}/usr/include)
  endif()
endif()

# Find and link the Threads library
find_package(Threads REQUIRED)
target_link_libraries(liteseq PRIVATE Threads::Threads)

# Log success if Threads library was linked
message(STATUS "Threads library linked successfully.")

# Include directories for the library
# expose the given directories to anything linking against liteseqlib
target_include_directories(liteseq
  PUBLIC
  $<BUILD_INTERFACE:${INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/generated>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
  ${SRC_INTERNAL_DIR}
)

# include(GNUInstallDirs)
# install(TARGETS liteseq
#   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#   LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
# )
# install(DIRECTORY include/liteseq/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/liteseq)


# message(STATUS "Install include dir: ${CMAKE_INSTALL_INCLUDEDIR}")
# message(STATUS "Install lib dir: ${CMAKE_INSTALL_LIBDIR}")
# message(STATUS "Install bin dir: ${CMAKE_INSTALL_BINDIR}")

# === Tests ===

if (LITESEQ_ENABLE_TESTING)
  # Add the global TESTING definition
  add_compile_definitions(TESTING)

  enable_testing()
  add_subdirectory(${TESTS_DIR})
endif()

# === Example Binary ===

if(LITESEQ_BUILD_EXAMPLE)
  add_subdirectory(${EXAMPLES_DIR})
endif()
