add_executable(test_povu main_tests.cc)
add_dependencies(test_povu povu)
target_compile_definitions(test_povu
  PRIVATE
  POVU_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)
target_link_libraries(test_povu
  PRIVATE
  gtest_main
  povulib  # Replace with your actual target
  ita
  mto
  liteseq      # gfa handling
  fmt::fmt     # formatting
  taywee::args # argument parsing
)

include(GoogleTest)
gtest_discover_tests(test_povu)
