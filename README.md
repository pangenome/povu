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

For the help text run `./bin/povu -h` or just `./bin/povu`

### Examples

#### 1

Using the dataset in the `test_data` folder we can call variants related to the
`HG02572__LPA__tig00000001` reference like so:

We can pass an individual reference

```
./bin/povu -v 2 call -i test_data/LPA.max120.gfa -- HG02572__LPA__tig00000001
```

This creates a `HG02572__LPA__tig00000001.vcf` which contains variants relative to that reference path in the graph.


#### 2

Extract all references from a GFA file and save in a file
```
grep '^P' test_data/HPRC/chrY.hprc-v1.0-pggb.gfa | cut -f 2 > haps.txt
```

Pass the references with `-p`
```
povu -v 4 call -i test_data/HPRC/chrY.hprc-v1.0-pggb.gfa -p haps.txt
```
This creates a vcf for each reference in the file `haps.txt`

## Input
Input GFA

Expect the segments in the input GFA to have unique numeric [segment names](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#s-segment-line)

## Development

Compile with debug symbols and with address sanitizer

```
cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER=address  -H. -Bbuild && cmake --build build -- -j 3
```
