set(ITA_INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(ITA_SOURCES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/ita)

ADD_LIBRARY(ita
  # align
  ${ITA_SOURCES_DIR}/align/align.cpp

  # genomics
  ${ITA_SOURCES_DIR}/genomics/genomics.cpp
  ${ITA_SOURCES_DIR}/genomics/vcf.cpp

  # trees
  ${ITA_SOURCES_DIR}/graph/slice_tree.cpp
  ${ITA_SOURCES_DIR}/graph/graph.cpp

  # convolutions
  ${ITA_SOURCES_DIR}/convolutions/convolutions.cpp
  ${ITA_SOURCES_DIR}/convolutions/trip.cpp

  # traversals
  ${ITA_SOURCES_DIR}/traversals/at_matrix.cpp
  ${ITA_SOURCES_DIR}/traversals/at_matrix_tangled.cpp
  ${ITA_SOURCES_DIR}/traversals/at_matrix_no_tangle.cpp
  ${ITA_SOURCES_DIR}/traversals/depth_matrix.cpp
  ${ITA_SOURCES_DIR}/traversals/traversals.cpp
  ${ITA_SOURCES_DIR}/traversals/untangle.cpp

  # variation
  ${ITA_SOURCES_DIR}/variation/sne.cpp
  ${ITA_SOURCES_DIR}/variation/rov.cpp
  ${ITA_SOURCES_DIR}/variation/color.cpp
)

target_include_directories(ita
  PUBLIC
  $<BUILD_INTERFACE:${ITA_INC_DIR}>
  $<INSTALL_INTERFACE:include>
)
