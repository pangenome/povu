add_executable(test_povu main_tests.cc)
target_link_libraries(test_povu
  PRIVATE
  gtest_main
  povulib  # Replace with your actual target
)

include(GoogleTest)
gtest_discover_tests(test_povu)
