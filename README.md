# povu &middot; [![Test Status](https://github.com/pangenome/povu/actions/workflows/test.yml/badge.svg?label=tests)](https://github.com/pangenome/povu/actions/workflows/test.yml) [![C++17 Build (CMake, macOS)](https://github.com/pangenome/povu/actions/workflows/macos.yml/badge.svg)](https://github.com/pangenome/povu/actions/workflows/macos.yml) [![C++17 Build (CMake, Ubuntu)](https://github.com/pangenome/povu/actions/workflows/build_cmake.yml/badge.svg)](https://github.com/pangenome/povu/actions/workflows/build_cmake.yml) ![c++17](https://img.shields.io/badge/C++-17-informational.svg?style=flat&logo=c%2B%2B&logoColor=white) [![License: MIT](https://img.shields.io/badge/License-MIT-informational.svg)](https://opensource.org/licenses/MIT)

A toolkit for exploring regions of genomic variation

## Table of Contents
- [Usage and Examples](#usage-and-examples)
- [Rust Bindings](#rust-bindings)
- [Building povu](#building-povu)
  - [Installing with Guix](#installing-with-guix)
  - [Building specific target](#building-specific-target)
  - [Development](#development)
- [Repository Structure](#repository-structure)


## Usage and Examples

For general help, run:

```bash
./bin/povu -h
# or simply
./bin/povu
```

The table below summarizes the subcommands currently available:

| Subcommand | Description                            |
|------------|--------------------------------------- |
| gfa2vcf    | Convert GFA to VCF (decompose + call)  |
| decompose  | Identify regions of variation          |
| call       | Generate VCF from regions of variation |
| info       | Print graph information                |

For detailed documentation on each subcommand, refer to the [docs/](./docs) directory.

### Quick Start: GFA to VCF

The simplest way to call variants is using `gfa2vcf`:

```bash
# Convert GFA to VCF using path prefix
./bin/povu gfa2vcf -i input.gfa -P HG > output.vcf

# Using a reference list file
./bin/povu gfa2vcf -i input.gfa -r ref_list.txt > output.vcf
```

The `gfa2vcf` command internally handles all intermediate steps and outputs a combined VCF to stdout.

### Two-Step Workflow

For more control, use the separate `decompose` and `call` commands:

```bash
# Step 1: Identify regions of variation
./bin/povu decompose -i input.gfa -o regions/

# Step 2a: Generate separate VCF files (one per reference)
./bin/povu call -i input.gfa -f regions/ -r ref_list.txt -o vcf_output/

# Step 2b: Generate single combined VCF to stdout
./bin/povu call -i input.gfa -f regions/ -r ref_list.txt --stdout > output.vcf
```

## Rust Bindings

Povu provides Rust bindings for embedding pangenome variation analysis in Rust applications. The bindings are located in the `povu-rs/` directory and offer both high-level convenience methods and detailed topology access.

### Quick Start (Rust)

Add to your `Cargo.toml`:
```toml
[dependencies]
povu = { git = "https://github.com/pangenome/povu", branch = "rust" }
```

### Example: Build a graph in memory

```rust
use povu::{PovuGraph, Orientation};

// Create an empty graph
let mut graph = PovuGraph::new(10, 15, 0);

// Add vertices
graph.add_vertex(1, "AAAA")?;
graph.add_vertex(2, "GGGG")?;
graph.add_vertex(3, "TTTT")?;
graph.add_vertex(4, "CCCC")?;

// Add edges to create a diamond (bubble)
graph.add_edge(1, Orientation::Forward, 2, Orientation::Forward)?;
graph.add_edge(1, Orientation::Forward, 3, Orientation::Forward)?;
graph.add_edge(2, Orientation::Forward, 4, Orientation::Forward)?;
graph.add_edge(3, Orientation::Forward, 4, Orientation::Forward)?;

// Finalize and analyze
graph.finalize();
println!("Graph has {} vertices and {} edges",
         graph.vertex_count(), graph.edge_count());
```

### Example: Load and analyze a GFA file

```rust
use povu::PovuGraph;

// Load from GFA
let graph = PovuGraph::load("graph.gfa")?;

// Query topology
let vertices = graph.vertices()?;
let edges = graph.edges()?;
let paths = graph.paths()?;

// Analyze for variation
let analysis = graph.analyze()?;
println!("Found {} variation regions", analysis.flubble_count());
```

For more examples, see the `povu-rs/examples/` directory.

## Building povu

Prerequisites:
- CMake (3.0+ recommended)
- C compiler (e.g., GCC or Clang)

**Clone the repository:**
```bash
git clone https://github.com/pangenome/povu.git
cd povu
```

**Standard build:**
```bash
cmake -H. -Bbuild && cmake --build build -- -j 3
```

The binary will be in `./bin/povu`

### Installing with Guix

**Enter Guix development shell:**
```bash
guix shell -C -N -D -f guix.scm
```

**Build inside the shell:**
```bash
cmake -H. -Bbuild -D CMAKE_BUILD_TYPE=Release && cmake --build build -- -j 3
```

### Building specific target

Building only the povu library

```
cmake -H. -DCMAKE_BUILD_TYPE=Debug -Bbuild && cmake --build build --target povulib -- -j 8
```

Building only the povu binary

```
cmake -H. -DCMAKE_BUILD_TYPE=Debug -Bbuild && cmake --build build --target povu -- -j 8
```

### Development

To compile povu with debug symbols and address sanitizer:

```bash
cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER=address -H. -Bbuild && cmake --build build -- -j 3
```

Running tests

1. Configure the build with testing enabled:

```
cmake -Bbuild -DPOVU_ENABLE_TESTING=ON -DCMAKE_BUILD_TYPE=Debug
```

2. Build the project:
```
cmake --build build
```

3. Run tests with CTest

```
ctest --test-dir build
```

## Repository Structure

```
povu/
├── bin/              # Compiled binaries (povu CLI)
├── docs/             # Documentation
├── include/          # C++ headers
│   └── povu/        # Public API headers
├── src/             # C++ implementation
│   └── povu/        # Core library sources
├── tests/           # C++ test data
├── povu-rs/         # Rust bindings
│   ├── src/         # Rust high-level API
│   ├── povu-ffi/    # C FFI bridge layer
│   │   ├── povu_ffi.h     # C API header
│   │   └── povu_ffi.cpp   # C++ implementation
│   ├── examples/    # Rust usage examples
│   └── tests/       # Rust integration tests
└── CMakeLists.txt   # Build configuration
```

### Key Components

- **C++ Core** (`src/`, `include/`): The main Povu library implementing the variation detection algorithms
- **CLI Tool** (`bin/`): Command-line interface for gfa2vcf, decompose, call, and info subcommands
- **Rust Bindings** (`povu-rs/`): Safe Rust wrappers around the C++ library
  - High-level API for graph construction, querying, and analysis
  - C FFI bridge providing a stable interface between C++ and Rust
  - Comprehensive test suite ensuring correctness

The Rust bindings live alongside the C++ code following the pattern used by many projects that provide multiple language interfaces (e.g., ripgrep, numpy, node.js).

## Name

The etymology of the name is rooted in profound philosophy 🤔. "Povu," is [Kiswahili](https://en.wikipedia.org/wiki/Swahili_language) for "foam." Foam, by nature, comprises countless flubbles.
