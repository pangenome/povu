set(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(QUILT_SOURCES_DIR ${SRC_DIR}/quilt)

ADD_LIBRARY(quilt
  ${QUILT_SOURCES_DIR}/graph_types.cpp
  ${QUILT_SOURCES_DIR}/utils.cpp
)

target_include_directories(quilt
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_compile_features(quilt PUBLIC cxx_std_17)
