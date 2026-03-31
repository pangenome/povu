# download CPM.cmake
file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.40.8/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
  EXPECTED_HASH SHA256=78ba32abdf798bc616bab7c73aac32a17bbd7b06ad9e26a6add69de8f3ae4791
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

# CPMAddPackage(
#   log
#   GITHUB_REPOSITORY urbanslug/log
#   GIT_TAG           7f1b62d70bb1b30ecfa86df55d2c82523e92fbfa
#   # Override with -DLOG_FULL_FILE_NAME=ON
#   OPTIONS
#   "LOG_FULL_FILE_NAME OFF"
# )

if (LITESEQ_ENABLE_TESTING)
CPMAddPackage(
  googletest
  GITHUB_REPOSITORY google/googletest
  GIT_TAG        v1.17.0
  VERSION        1.17.0
  OPTIONS
  "BUILD_GMOCK OFF"
  "INSTALL_GTEST OFF"
  "gtest_force_shared_crt ON"
)
endif()
