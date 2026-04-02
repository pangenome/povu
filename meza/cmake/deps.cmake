# download CPM.cmake
file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.40.8/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
  EXPECTED_HASH SHA256=78ba32abdf798bc616bab7c73aac32a17bbd7b06ad9e26a6add69de8f3ae4791
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

if (MEZA_ENABLE_TESTS)
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
