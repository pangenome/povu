# download CPM.cmake
file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.40.8/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
  EXPECTED_HASH SHA256=78ba32abdf798bc616bab7c73aac32a17bbd7b06ad9e26a6add69de8f3ae4791
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

CPMAddPackage(
  NAME              args
  GITHUB_REPOSITORY Taywee/args
  GIT_TAG           114200a9ad5fe06c8dea76e15d92325695cf3e34 # v6.4.7
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

# CPMAddPackage(
#   NAME              liteseq
#   GITHUB_REPOSITORY urbanslug/liteseq
#   GIT_TAG           3678c8b7d1c5e1bfaa9e3f07133e6e3f00567015
#   OPTIONS
#   "LITESEQ_BUILD_EXAMPLE OFF"
#   "LITESEQ_ENABLE_TESTING OFF"
# )

# CPMAddPackage(
#   NAME              quilt
#   GITHUB_REPOSITORY urbanslug/quilt
#   GIT_TAG           0d0c7e238887518a2e602d11b0f0d16cddfce9f0 # v0.0.1
# )

# CPMAddPackage(
#   NAME              meza
#   SOURCE_DIR        /home/sluggie/src/pangenomics/meza
#   OPTIONS
#   "MEZA_ENABLE_TESTS OFF"
#   "MEZA_BUILD_EXAMPLE OFF"
#   "MEZA_ENABLE_CUDA ON"
# )

# CPMAddPackage(
#   NAME              meza
#   GITHUB_REPOSITORY urbanslug/meza
#   GIT_TAG           c95166f3fad8bbef697dfb91d44a9301aad485a8
#   OPTIONS
#   "MEZA_ENABLE_TESTS OFF"
#   "MEZA_BUILD_EXAMPLE OFF"
#   "MEZA_ENABLE_CUDA ON"
# )

if (POVU_ENABLE_TESTING)
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
