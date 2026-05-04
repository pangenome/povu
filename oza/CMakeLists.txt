set(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(OZA_SOURCES_DIR ${SRC_DIR}/povu)
set(TESTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)

ADD_LIBRARY(oza
  # algorithms
  ${OZA_SOURCES_DIR}/algorithms/concealed.cpp
  ${OZA_SOURCES_DIR}/algorithms/midi.cpp
  ${OZA_SOURCES_DIR}/algorithms/smothered.cpp
  ${OZA_SOURCES_DIR}/algorithms/parallel.cpp
  ${OZA_SOURCES_DIR}/algorithms/tiny.cpp
  ${OZA_SOURCES_DIR}/algorithms/flubbles.cpp

  ${OZA_SOURCES_DIR}/refs/refs.cpp

  # graph
  ${OZA_SOURCES_DIR}/graph/bidirected.cpp
  ${OZA_SOURCES_DIR}/graph/bracket_list.cpp
  ${OZA_SOURCES_DIR}/graph/spanning_tree.cpp
  ${OZA_SOURCES_DIR}/graph/tree_utils.cpp
)

target_include_directories(oza
  PUBLIC
  $<BUILD_INTERFACE:${INC_DIR}>
  $<INSTALL_INTERFACE:include>
)


target_link_libraries(oza
  PRIVATE
  log
)

if (TEST_OZA)
  message(STATUS "Tests enabled for oza")
  # Add the global TESTING definition
  add_compile_definitions(TESTING)

  enable_testing()
  add_subdirectory(${TESTS_DIR})
endif()
