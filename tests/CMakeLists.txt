add_executable(test_povu main_tests.cc)
target_link_libraries(test_povu
  PRIVATE
  gtest_main
  povulib  # Replace with your actual target
  liteseq      # gfa handling
  fmt::fmt     # formatting
  taywee::args # argument parsing
  indicators   # progress bars
)

include(GoogleTest)
gtest_discover_tests(test_povu)
