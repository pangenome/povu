cmake_minimum_required(VERSION 3.14)

project(dynamo VERSION 0.0.1 LANGUAGES CXX)

add_library(dynamo INTERFACE)

target_include_directories(dynamo
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_compile_features(dynamo INTERFACE cxx_std_17)
