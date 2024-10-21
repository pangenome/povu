cmake_minimum_required(VERSION 3.14)
#set(CMAKE_C_COMPILER "gcc")
#set(CMAKE_CXX_COMPILER "g++")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  # set(CMAKE_BUILD_TYPE Debug)
endif()

# Set Linker flags
#set(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++")

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DRELEASE")

project(povu)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake 3.24
# https://cmake.org/cmake/help/latest/policy/CMP0135.html
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()


include(sanitizers.cmake)

if (MSVC)
    # warning level 4
    add_compile_options(/W4)
else()
    # additional warnings
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()


set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# add_subdirectory( deps/gfakluge )
add_subdirectory( deps/libhandlegraph )
add_subdirectory( deps/WFA2-lib )

ADD_LIBRARY(LibsModule
  # io
  src/io/bub.cpp
  src/io/from_gfa.cpp
  src/io/txt.cpp
  src/io/vcf.cpp

  # genomics
  src/genomics/untangle.cpp
  src/genomics/variant_calling.cpp
  src/genomics/vcf.cpp

  # graph
  src/graph/bidirected.cpp
  src/graph/biedged.cpp
  src/graph/flubble_tree.cpp
  src/graph/graph.cpp
  src/graph/spanning_tree.cpp
  src/graph/bracket_list.cpp

  # common
  src/common/types.cpp
  src/common/utils.cpp

  # algorithms
  src/algorithms/algorithms.cpp
    src/algorithms/align.cpp

  # subcommand
  src/subcommand/deconstruct.cpp

  # cli
  src/cli/cli.cpp
  src/cli/app.cpp
)


include_directories(
  # deps
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/args/
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/gfakluge/src
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/gfakluge/src/tinyFA
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/libhandlegraph/src/include
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/WFA2-lib
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/WFA2-lib/system
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/WFA2-lib/wavefront
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/WFA2-lib/bindings/cpp
)


set(SOURCE_FILES
  src/main.cpp
)


# source
add_executable(povu
  ${SOURCE_FILES}
)

# Link the library (LibsModule) and libhandlegraph to your executable
target_link_libraries(povu PRIVATE LibsModule handlegraph_shared wfa2cpp)

set(BINARY_DIR ./bin)

file(MAKE_DIRECTORY ${BINARY_DIR}/)


add_custom_command(TARGET povu
		   POST_BUILD
		   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:povu> ../${BINARY_DIR}/)
