cmake_minimum_required(VERSION 3.14)

project(quilt VERSION 0.0.1 LANGUAGES CXX)

add_library(quilt INTERFACE)

target_include_directories(quilt
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_compile_features(quilt INTERFACE cxx_std_17)

# # download CPM.cmake
# file(
#   DOWNLOAD
#   https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.40.8/CPM.cmake
#   ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
#   EXPECTED_HASH SHA256=78ba32abdf798bc616bab7c73aac32a17bbd7b06ad9e26a6add69de8f3ae4791
# )

# include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

# CPMAddPackage(
#   NAME fmt
#   GITHUB_REPOSITORY fmtlib/fmt
#   GIT_TAG           12.1.0
#   VERSION           12.1.0
#   OPTIONS
#   "CMAKE_POSITION_INDEPENDENT_CODE ON"
# )
