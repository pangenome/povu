# meza

`meza` is a matrix library with support for both owning matrices and non-owning 
matrix views, with optional NVIDIA CUDA acceleration.

## Features

- Owning matrix types
- Non-owning matrix views
- Pool-based allocation utilities
- Optional CUDA support for GPU-enabled builds

## Requirements

To build `meza`, you will need:

- CMake
- A C++ compiler
- Optionally, an NVIDIA GPU and CUDA toolkit

## Building

### CUDA build

Configure:

with an NVIDIA GPU
```bash
cmake -B build -S . -DMEZA_ENABLE_CUDA=ON -DCMAKE_BUILD_TYPE=Debug
```

CPU only build
```bash
cmake -B build -S . -DMEZA_ENABLE_CUDA=OFF -DCMAKE_BUILD_TYPE=Debug
```


Build
```
cmake --build build
```
