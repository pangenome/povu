set(CPM_SOURCE_CACHE
  "${CMAKE_SOURCE_DIR}/.cpm-cache"
  CACHE PATH "Directory for caching CPM sources"
)

# download CPM.cmake
file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.40.8/CPM.cmake
  ${CPM_SOURCE_CACHE}/CPM.cmake
  EXPECTED_HASH SHA256=78ba32abdf798bc616bab7c73aac32a17bbd7b06ad9e26a6add69de8f3ae4791
)
include(${CPM_SOURCE_CACHE}/CPM.cmake)

set(CPM_USE_LOCAL_PACKAGES ON)

CPMAddPackage(
  NAME              args
  GITHUB_REPOSITORY Taywee/args
  GIT_TAG           6.4.7
  VERSION           6.4.7
  OPTIONS
  "ARGS_BUILD_UNITTESTS OFF"
)

CPMAddPackage(
  NAME              fmt
  GITHUB_REPOSITORY fmtlib/fmt
  GIT_TAG           12.1.0
  VERSION           12.1.0
  OPTIONS
  "CMAKE_POSITION_INDEPENDENT_CODE ON" # Ensure -fPIC is enabled
)

if (ENABLE_ALL_TESTS)
CPMAddPackage(
  NAME              googletest
  GITHUB_REPOSITORY google/googletest
  GIT_TAG           v1.17.0
  VERSION           1.17.0
  OPTIONS
  "BUILD_GMOCK OFF"
  "INSTALL_GTEST OFF"
  "gtest_force_shared_crt ON"
)
endif()
