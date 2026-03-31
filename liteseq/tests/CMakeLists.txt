# Create the test executable
add_executable(test_liteseq
  main_tests.cc  # Replace this with your actual test source files
)

# Add include directories required for tests
target_include_directories(test_liteseq
  PRIVATE
  ${SRC_DIR}          # Grants access to the src directory
  ${SRC_INTERNAL_DIR} # Grants access to the src/internal directory
)


target_link_libraries(test_liteseq
  PRIVATE
  gtest_main
  liteseq  # Replace with your actual target
  log
)

# Discover GoogleTest tests
include(GoogleTest)
gtest_discover_tests(test_liteseq)
