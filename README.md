# povu
A variant caller based on cycle equivalence

## Install

### Compile from source

1. Fetch the source code
```
git clone --recursive git@github.com:pangenome/domibubble.git
```

2. Compile
```
cmake -H. -Bbuild && cmake --build build -- -j 3
```

3. The binary should be in `./bin/povu`


## Run

Using the dataset in the test_data folder we can call variants related to the `HG02572__LPA__tig00000001` reference
like so:

```
./bin/povu -v 2 call -i test_data/LPA.max120.gfa  --  HG02572__LPA__tig00000001
```

This creates a `HG02572__LPA__tig00000001.vcf` which contains variants relative to that reference path in the graph.

For the help text run `./bin/povu -h` or just `./bin/povu`

## Input
Input GFA

Expect the segments in the input GFA to have unique numeric [segment names](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#s-segment-line)

## Development

Compile with debug symbols and with address sanitizer

```
cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER=address  -H. -Bbuild && cmake --build build -- -j 3
```
