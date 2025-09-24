# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

povu is a lightweight C++ tool for exploring regions of genomic variation. The name "povu" means "foam" in Kiswahili, referencing the foam's composition of countless "flubbles" (a novel concept for representing genomic variation).

## Build Commands

### Standard build
```bash
cmake -H. -Bbuild && cmake --build build -- -j 3
```

### Debug build with sanitizers
```bash
cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER=address -H. -Bbuild && cmake --build build -- -j 3
```

The compiled binary is automatically copied to `./bin/povu` after building.

### Clean build
```bash
rm -rf build bin/povu
```

## Architecture

### Core Concepts
- **Flubbles**: The fundamental unit of genomic variation in povu
- **PVST (Povu Variation Structure Tree)**: The main data structure for representing genomic variation
- **Graph-based representation**: Uses GFA format to represent genomic sequences as graphs

### Module Structure
- `include/algorithms/` - Core algorithms for flubble detection and processing
- `include/align/` - Alignment-related functionality
- `include/genomics/` - Genomic-specific operations (extract_canonical_flubbles)
- `include/graph/` - Graph operations and structures
- `include/common/` - Shared utilities and data structures
- `src/subcommand/` - Implementation of CLI subcommands (deconstruct, call, info)
- `src/io/` - Input/output operations for various file formats

### Key Dependencies
- `args`: Command-line argument parsing (git submodule in deps/args)
- `liteseq`: Lightweight sequence handling library (git submodule in deps/liteseq)

### File Formats
- **Input**: GFA (Graphical Fragment Assembly) files with unique numeric segment names
- **Output**: 
  - `.flb` files: 3-column flubble tree format (parent_id, child_id, descriptor)
  - VCF files: Standard variant call format

## Testing

Currently, the project has test files written with Google Test framework but lacks build configuration for running them. Test files exist in:
- `tests/compute_pvst.cc` - PVST functionality tests
- `tests/genomics.cc` - Genomics module tests
- `test_data/` - Contains synthetic and real test data

To run tests, the CMakeLists.txt would need to be updated to include Google Test integration.

## Development Notes

- C++20 standard is required
- Extensive compiler warnings are enabled (`-Wall -Wextra -Wpedantic`)
- Debug builds include additional warnings and can use various sanitizers (address, memory, undefined, thread, leak)
- No automated code formatting or linting tools are currently configured