# povu

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Test Status (CMake)](https://github.com/pangenome/povu/actions/workflows/build_cmake.yml/badge.svg?label=build)](https://github.com/pangenome/povu/actions/workflows/build_cmake.yml)
[![Test Status](https://github.com/pangenome/povu/actions/workflows/test.yml/badge.svg?label=tests)](https://github.com/pangenome/povu/actions/workflows/test.yml)


A toolkit for exploring regions of genomic variation

## Table of Contents
- [Usage and Examples](#usage-and-examples)
- [Building povu](#building-povu)
  - [Installing with Guix](#installing-with-guix)
  - [Building specific target](#building-specific-target)
  - [Development](#development)


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

## Name

The etymology of the name is rooted in profound philosophy 🤔. "Povu," is [Kiswahili](https://en.wikipedia.org/wiki/Swahili_language) for "foam." Foam, by nature, comprises countless flubbles.
