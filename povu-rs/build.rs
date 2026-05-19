use std::env;
use std::fs;
use std::path::{Path, PathBuf};

fn emit_rerun_if_changed(path: impl AsRef<Path>) {
    let path = path.as_ref();
    if path.is_file() {
        println!("cargo:rerun-if-changed={}", path.display());
        return;
    }

    if !path.is_dir() {
        return;
    }

    let Ok(entries) = fs::read_dir(path) else {
        return;
    };

    for entry in entries.flatten() {
        emit_rerun_if_changed(entry.path());
    }
}

fn main() {
    println!("cargo:rerun-if-env-changed=CARGO_FEATURE_FFI");
    if env::var_os("CARGO_FEATURE_FFI").is_none() {
        return;
    }

    // Build the C++ library using CMake
    let dst = cmake::Config::new("..")
        .define("POVU_ENABLE_TESTING", "OFF")
        .define("BUILD_SHARED_LIBS", "OFF")
        .define("POVU_BUILD_FFI", "ON")
        .build();

    // Tell cargo to tell rustc to link the povu libraries
    println!(
        "cargo:rustc-link-search=native={}/build/povu-rs/povu-ffi",
        dst.display()
    );
    println!("cargo:rustc-link-search=native={}/build", dst.display());
    println!("cargo:rustc-link-lib=static=povu_ffi");
    // The FFI GFA loader calls the migrated mto::from_gfa API, so Cargo must
    // link libmto explicitly in addition to the core povulib archive.
    println!("cargo:rustc-link-lib=static=mto");
    println!("cargo:rustc-link-lib=static=povulib");

    // Link dependencies
    println!(
        "cargo:rustc-link-search=native={}/build/_deps/liteseq-build",
        dst.display()
    );
    println!(
        "cargo:rustc-link-search=native={}/build/_deps/fmt-build",
        dst.display()
    );
    println!(
        "cargo:rustc-link-search=native={}/build/_deps/log-build",
        dst.display()
    );
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

    // Tell Cargo to rerun the native build when the C++ library or wrapper changes.
    for path in [
        "../CMakeLists.txt",
        "../cmake",
        "../include",
        "../src",
        "povu-ffi",
    ] {
        emit_rerun_if_changed(path);
    }

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
