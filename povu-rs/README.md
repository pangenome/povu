# povu-rs

Rust bindings for [Povu](https://github.com/pangenome/povu), a toolkit for exploring regions of genomic variation in pangenomes.

## Overview

`povu-rs` provides both low-level FFI bindings and a high-level idiomatic Rust API for working with pangenome variation graphs. It's designed for embedding Povu functionality into Rust-based pangenome analysis tools.

## Features

- Load GFA files into bidirected pangenome graphs
- Query graph topology (vertices, edges, reference paths)
- Detect flubbles (regions of variation/bubbles)
- Access hierarchical PVST (Pangenome Variation Structure Tree)
- Generate VCF variant calls
- Thread-safe API

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
povu-rs = { path = "povu-rs" }
```

### Build Requirements

- Rust 1.70+ (2021 edition)
- CMake 3.14+
- C++17 compatible compiler
- clang/libclang (for bindgen)

## Quick Start

### Simple Usage

```rust
use povu::PovuGraph;

// Load a graph
let graph = PovuGraph::load("input.gfa")?;

// Query topology
println!("Graph has {} vertices and {} edges",
         graph.vertex_count(), graph.edge_count());

// Get vertices, edges, and paths
let vertices = graph.vertices()?;
let edges = graph.edges()?;
let paths = graph.paths()?;

// Analyze for variation
let analysis = graph.analyze()?;
println!("Found {} flubbles", analysis.flubble_count());
```

### Working with Graph Topology

```rust
use povu::PovuGraph;

let graph = PovuGraph::load("input.gfa")?;

// Iterate through vertices
for vertex in graph.vertices()? {
    println!("Vertex {}: {} bp", vertex.id, vertex.len());
}

// Examine edges
for edge in graph.edges()? {
    println!("Edge: {}{} -> {}{}",
             edge.from_id, edge.from_orientation,
             edge.to_id, edge.to_orientation);
}

// Inspect reference paths
for path in graph.paths()? {
    if let Some(sample) = path.sample_name() {
        println!("Sample: {}, Haplotype: {:?}",
                 sample, path.haplotype());
    }
}
```

### Variation Analysis

```rust
use povu::PovuGraph;

let mut graph = PovuGraph::load("input.gfa")?;

// Select references by prefix
graph.set_references_from_prefixes(&["HG002", "HG003"])?;

// Analyze variation
let analysis = graph.analyze()?;

// Access PVST tree
let pvst = analysis.pvst_tree();
println!("PVST has {} vertices", pvst.vertex_count());

// Generate VCF (when implemented)
// analysis.write_vcf("output.vcf")?;
```

## Examples

Run the included examples:

```bash
# Basic graph statistics
cargo run --example simple tests/data/simple.gfa

# Detailed topology analysis
cargo run --example topology tests/data/simple.gfa
```

## Testing

Run the test suite (compares Rust outputs to native Povu):

```bash
cargo test
```

Tests require GFA test files in `../tests/data/`.

## Architecture

The library is structured in three layers:

1. **FFI Layer** (`povu-ffi/`): C-compatible bridge to C++ Povu library
2. **Low-level Bindings** (`src/ffi.rs`): Auto-generated bindings via bindgen
3. **High-level API** (`src/*.rs`): Idiomatic Rust interface

### Build Process

1. CMake builds the `povulib` static library
2. CMake builds the `povu_ffi` C bridge library
3. Cargo's build script (`build.rs`) invokes CMake
4. bindgen generates Rust FFI bindings from `povu_ffi.h`
5. High-level Rust API wraps the FFI bindings

## API Documentation

Generate and view the full API documentation:

```bash
cargo doc --open
```

## Limitations & TODOs

- VCF generation is not yet fully implemented in the FFI layer
- PVST tree traversal APIs are minimal
- Detailed flubble iteration not yet exposed

## Contributing

This is part of the [Povu](https://github.com/pangenome/povu) project. Contributions welcome!

## License

Same as Povu (check parent repository)

## Citation

If you use povu-rs in your research, please cite the original Povu paper (when published).
