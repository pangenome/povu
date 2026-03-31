# Create the test executable
add_executable(test_meza
  main_tests.cc
)

# Add include directories required for tests
target_include_directories(test_meza
  PRIVATE
  ${SRC_DIR}          # Grants access to the src directory
)

target_link_libraries(test_meza
  PRIVATE
  quilt
  meza
  gtest_main
)

# Discover GoogleTest tests
include(GoogleTest)
gtest_discover_tests(test_meza)
