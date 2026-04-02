set(MTO_SOURCES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/mto)

ADD_LIBRARY(mto # io library
  ${MTO_SOURCES_DIR}/from_pvst.cpp
  ${MTO_SOURCES_DIR}/to_pvst.cpp
  ${MTO_SOURCES_DIR}/from_gfa.cpp
  ${MTO_SOURCES_DIR}/to_gfa.cpp
  ${MTO_SOURCES_DIR}/common.cpp
  ${MTO_SOURCES_DIR}/to_vcf.cpp
  ${MTO_SOURCES_DIR}/from_vcf.cpp
)

target_include_directories(ita
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
