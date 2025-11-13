use std::env;
use std::path::PathBuf;

fn main() {
    // Build the C++ library using CMake
    let dst = cmake::Config::new("..")
        .define("POVU_ENABLE_TESTING", "OFF")
        .define("BUILD_SHARED_LIBS", "OFF")
        .define("POVU_BUILD_FFI", "ON")
        .build();

    // Tell cargo to tell rustc to link the povu libraries
    println!("cargo:rustc-link-search=native={}/build/povu-rs/povu-ffi", dst.display());
    println!("cargo:rustc-link-search=native={}/build", dst.display());
    println!("cargo:rustc-link-lib=static=povu_ffi");
    println!("cargo:rustc-link-lib=static=povulib");

    // Link dependencies
    println!("cargo:rustc-link-search=native={}/build/_deps/liteseq-build", dst.display());
    println!("cargo:rustc-link-search=native={}/build/_deps/fmt-build", dst.display());
    println!("cargo:rustc-link-search=native={}/build/_deps/log-build", dst.display());
    println!("cargo:rustc-link-lib=static=liteseq");
    println!("cargo:rustc-link-lib=static=fmtd");
    println!("cargo:rustc-link-lib=static=log");

    // Link C++ standard library
    let target = env::var("TARGET").unwrap();
    if target.contains("apple") {
        println!("cargo:rustc-link-lib=dylib=c++");
    } else if target.contains("linux") {
        println!("cargo:rustc-link-lib=dylib=stdc++");
    }

    // Tell cargo to invalidate the built crate whenever the wrapper changes
    println!("cargo:rerun-if-changed=povu-ffi/povu_ffi.h");
    println!("cargo:rerun-if-changed=povu-ffi/povu_ffi.cpp");

    // Generate bindings
    let bindings = bindgen::Builder::default()
        .header("povu-ffi/povu_ffi.h")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .allowlist_function("povu_.*")
        .allowlist_type("Povu.*")
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
