add_executable(test_oza
  main_tests.cc
)

target_link_libraries(test_oza
  PRIVATE

  gtest
  gtest_main

  fmt::fmt
  taywee::args

  oza
  liteseq
  meza
)

include(GoogleTest)
gtest_discover_tests(test_oza
  TEST_PREFIX "oza."
)
