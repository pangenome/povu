cmake_minimum_required(VERSION 3.14)
project(domibubble)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
enable_testing()

FetchContent_Declare(
  libhandlegraph
  URL https://github.com/vgteam/libhandlegraph/archive/refs/heads/master.zip
)
FetchContent_MakeAvailable(libhandlegraph)

if (MSVC)
    # warning level 4
    add_compile_options(/W4)
else()
    # additional warnings
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

ADD_LIBRARY(LibsModule
  src/graph/tree.cpp
  src/graph/graph.cpp
  src/graph/digraph.cpp
  src/graph/spanning_tree.cpp
  src/vst/vst.cpp
  src/vst/pst.cpp
)

set(SOURCE_FILES
  src/main.cpp
)

# source
add_executable( domibubble
${SOURCE_FILES}
)

target_link_libraries(
  domibubble
  LibsModule)


set(BINARY_DIR ./bin)

file(MAKE_DIRECTORY ${BINARY_DIR}/)


add_custom_command(TARGET domibubble
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:domibubble> ../${BINARY_DIR}/)

add_executable( domibubble_test
  tests/main_tests.cc
)


target_link_libraries(
  domibubble_test
  LibsModule
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(domibubble_test)
