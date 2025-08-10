# povu

A lightweight tool for exploring regions of genomic variation

## Table of Contents
 - [Usage and Examples](#usage-and-examples)
 - [Building povu](#building-povu)
   * [Building specific target](#building-specific-target)
   * [Development](#development)


## Usage and Examples

For general help, run:

```bash
./bin/povu -h
# or simply
./bin/povu
```

The table below summarizes the subcommands currently available:


| Subcommand | Description                                 |
|------------|---------------------------------------------|
| decompose | Identifies regions of variation in the graph |
| call      | Call variants                                |
| info      | Provides a summary of the input GFA          |


For detailed documentation on each subcommand, refer to the [docs/](./docs) directory.



## Building povu
Prerequisites:

  - CMake (3.0+ recommended)
  - C compiler (e.g., GCC or Clang)


1. **Fetch the source code**
```
git clone --recursive https://github.com/urbanslug/povu.git
```

2. **Compile**
```
cmake -H. -Bbuild && cmake --build build -- -j 3
```

The binary should be in `./bin/povu`

### Installing with Guix

1. **Clone the repository**
```
git clone --recursive https://github.com/urbanslug/povu.git
cd povu
```

2. **Enter Guix development shell**
```
guix shell -C -N -D -f guix.scm
```

3. **Build inside the shell**
```
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

To compile povu with debug symbols and with address sanitizer

```
cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER=address  -H. -Bbuild && cmake --build build -- -j 3
```


## Name

The etymology of the name is rooted in profound philosophy 🤔. "Povu," is [Kiswahili](https://en.wikipedia.org/wiki/Swahili_language) for "foam." Foam, by nature, comprises countless flubbles.
