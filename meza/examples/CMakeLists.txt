set(BIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# Source files for binary
set(SOURCES
${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
${CMAKE_CURRENT_SOURCE_DIR}/hap_xors.cpp
${CMAKE_CURRENT_SOURCE_DIR}/shared.cpp
)

# create the BIN_DIR if it doesn't exist
file(MAKE_DIRECTORY ${BIN_DIR})

# Add the executable
add_executable(meza-example ${SOURCES})

# Link the library to the executable
target_link_libraries(meza-example
  PRIVATE
  meza
)

# Set the runtime output directory for the example binary
set_target_properties(meza-example PROPERTIES
RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
)

# Install the example binary
install(TARGETS meza-example
RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
