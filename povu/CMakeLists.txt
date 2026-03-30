set(POVU_INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(POVU_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/povu)

add_executable(povu
  # cli
  ${POVU_SRC_DIR}/cli.cpp

  # subcommands
  ${POVU_SRC_DIR}/decompose.cpp
  ${POVU_SRC_DIR}/info.cpp
  ${POVU_SRC_DIR}/call.cpp
  ${POVU_SRC_DIR}/prune.cpp
  ${POVU_SRC_DIR}/gfa2vcf.cpp
  ${POVU_SRC_DIR}/vcf.cpp
  ${POVU_SRC_DIR}/view.cpp

  # main
  ${POVU_SRC_DIR}/main.cpp
)

target_include_directories(povu
  PRIVATE
  ${POVU_INC_DIR}
)
