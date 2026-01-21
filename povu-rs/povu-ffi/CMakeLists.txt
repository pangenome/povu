# FFI library for Rust bindings
add_library(povu_ffi STATIC
    povu_ffi.cpp
)

target_include_directories(povu_ffi PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(povu_ffi PUBLIC
    povulib
    fmt::fmt
    liteseq
)

# Ensure C++17 standard
target_compile_features(povu_ffi PUBLIC cxx_std_17)
