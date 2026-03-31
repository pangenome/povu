set(BIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# Source files for binary
set(SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example.c
)

# create the BIN_DIR if it doesn't exist
file(MAKE_DIRECTORY ${BIN_DIR})

# Add the executable
add_executable(liteseq-example ${SOURCES})

# Link the library to the executable
target_link_libraries(liteseq-example PRIVATE liteseq log)

# Set the runtime output directory for the example binary
set_target_properties(liteseq-example PROPERTIES
RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
)

# Install the example binary
install(TARGETS liteseq-example
RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
