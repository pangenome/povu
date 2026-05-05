# === CMake Setup ===
cmake_minimum_required(VERSION 3.14)

# Project name, version and languages
project(log LANGUAGES C CXX)

# === BUILD OPTIONS ===
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
# Check if CMAKE_BUILD_TYPE is set; if not, default to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build (Debug, Release)" FORCE)
endif()

# Automatically enable the generation of compile_commands.json in Debug mode.
# The compile_commands.json file is useful for tools like clangd or IDEs
# that rely on compile commands for features like autocompletion and diagnostics.
# In Release mode (or other build types), this feature is turned off by default
# to avoid unnecessary output, but it can still be manually enabled by users.
# For other build types (e.g., Release), it defaults to `OFF` unless explicitly enabled by the user.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Enable compile_commands.json generation" FORCE)
else()
  set(CMAKE_EXPORT_COMPILE_COMMANDS OFF CACHE BOOL "Disable compile_commands.json generation" FORCE)
endif()

# === Platform & Toolchain Configuration ===
include(GNUInstallDirs)
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SAMPLES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/samples)
set(BIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)


# --- log Library ---
add_library(log
  ${SRC_DIR}/log.c
)

# C standard & basic props
set_target_properties(log PROPERTIES
  C_STANDARD 99
  CXX_STANDARD 11
)

# Compiler options
# Debug-specific flags

set(COMMON_WARNINGS
  -Wall
)

set(DEBUG_WARNINGS
  -Wextra
  -Wpedantic
  -Wshadow
  -Wunused
  -Wunknown-pragmas
  -Wunused-parameter
  -Wunused-variable
  -Wunused-function
  -Wunused-but-set-variable
  -Wmissing-field-initializers
  -Wunused-local-typedefs
  -Wuninitialized
  -Wmaybe-uninitialized
  -Wreturn-type
  # -Wpadded # Uncomment or gate by compiler if needed
)


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(-g)
  add_compile_options(${COMMON_WARNINGS})
  add_compile_options(${DEBUG_WARNINGS})
  add_compile_definitions(DEBUG)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_compile_options(${COMMON_WARNINGS})
  add_compile_options(-O3)
  add_compile_definitions(RELEASE)
endif()

# ===================== Threading Configuration ======================

message("Looking for an appropriate threading library...")
if (APPLE)
  # Log the start of macOS SDK retrieval
  message(STATUS "Retrieving macOS SDK path using xcrun...")

  # Get the macOS SDK path dynamically
  execute_process(
    COMMAND xcrun --show-sdk-path
    OUTPUT_VARIABLE MACOSX_SDK_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Check if the macOS SDK path was successfully retrieved
  if (NOT MACOSX_SDK_PATH)
    message(FATAL_ERROR "Failed to retrieve macOS SDK path. Ensure that Xcode tools are installed and configured.")
  else()
    # Log the retrieved macOS SDK path
    message(STATUS "Using macOS SDK path: ${MACOSX_SDK_PATH}")

    # Add the SDK's /usr/include directory to the CMake include paths
    message(STATUS "Adding ${MACOSX_SDK_PATH}/usr/include to the include directories.")
    include_directories(${MACOSX_SDK_PATH}/usr/include)
  endif()
endif()

# Find and link the Threads library
find_package(Threads REQUIRED)
target_link_libraries(log
  PRIVATE
  Threads::Threads)

# Log success if Threads library was linked
message(STATUS "Threads library linked successfully.")

# =========================== Logging Configuration ============================

# Allow users to provide their own settings by defining LOGC__USER_SETTINGS
add_compile_definitions(LOGC__USER_SETTINGS)

# Enable color output by default; can be disabled by setting LOG_USE_COLOR to OFF
option(LOG_USE_COLOR "Enable color output" ON)
if(LOG_USE_COLOR)
  target_compile_definitions(log PRIVATE LOG_USE_COLOR)
endif()

option(LOG_FULL_FILE_NAME "Log entire file path" OFF)
if(LOG_FULL_FILE_NAME)
  target_compile_definitions(log PRIVATE LOG_FULL_FILE_NAME)
endif()

option(LOG_FILE_INFO "Log file and line number" OFF)
if(LOG_FILE_INFO)
  target_compile_definitions(log PRIVATE LOG_FILE_INFO)
endif()

# --- Set log output stream ---
# Choice between stdout and stderr
set(LOG_DEFAULT_STREAM "stderr" CACHE STRING "Default stream for log output (stdout/stderr)")
# Validate value
set_property(CACHE LOG_DEFAULT_STREAM PROPERTY STRINGS stdout stderr)
# Expose as a macro
target_compile_definitions(log PRIVATE LOG_DEFAULT_STREAM=${LOG_DEFAULT_STREAM})

# ==================== Include Directories Configuration =======================

# Include directories for the library
# expose the given directories to anything linking against loglib
target_include_directories(log
  PUBLIC
  ${INCLUDE_DIR}
)


if(LOG_BUILD_SAMPLES)
  add_executable(log_callbacks
    ${SAMPLES_DIR}/callbacks/main.c
  )

  add_executable(log_multithreading
    ${SAMPLES_DIR}/threading/main.c
  )

   add_executable(log_basic
     ${SAMPLES_DIR}/basic/main.c
   )

  # create the BIN_DIR if it doesn't exist
  file(MAKE_DIRECTORY ${BIN_DIR})

  target_link_libraries(log_multithreading PRIVATE log)
  target_link_libraries(log_callbacks PRIVATE log)
  target_link_libraries(log_basic PRIVATE log)

  set_target_properties(log_multithreading PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
  )
  set_target_properties(log_callbacks PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
  )
  set_target_properties(log_basic PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
  )
endif()
