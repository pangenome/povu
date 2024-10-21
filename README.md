# povu
Find regions of variation in a graph by finding Eulerian cycles.


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


## Usage

For general help text run `./bin/povu -h` or just `./bin/povu`

### Call

**help:** `./bin/povu call`

The `call` sub-command takes a graph and a set of reference/haplotype paths and returns the variants
relative to each haplotype in a VCF format.


### Examples

#### Example 1

*Passing a haplotype directly.*

Using the dataset in the `./test_data` folder we can call variants related to the
`HG02572__LPA__tig00000001` haplotype like so:


```
./bin/povu -v 2 call -i test_data/LPA.max120.gfa -- HG02572__LPA__tig00000001
```

This creates a `HG02572__LPA__tig00000001.vcf` which contains variants relative to that reference path in the graph.


#### Example 2

*Passing a set of line separated haplotypes in a file.*

Extract all references from a GFA file and save them in a file
```
grep '^P' test_data/HPRC/chrY.hprc-v1.0-pggb.gfa | cut -f 2 > haps.txt
```

Pass the references with `-p`
```
povu -v 4 call -i test_data/HPRC/chrY.hprc-v1.0-pggb.gfa -p haps.txt
```

This creates a VCF file for each haplotype in the `haps.txt`


### Deconstruct

**help:** `./bin/povu deconstruct`

For any graph, the `deconstruct` sub-command generates the bubble forest: a set of bubble trees. Currently, each bubble tree is stored in its own bubble file.

#### Bub Format

The bub format is a plain-text file representing bubble tree(s) and ends with the `.bub` extension. It is tab-separated and consists of 3 columns where a dot `.` represents a null value.
Below is a table with the specifics of each column.

| column    | type             | description                                                                                                                                                                                   |
|-----------|------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| vertex id | unsigned numeric | A unique identifier (can be zero & cannot be null).                                                                                                                                           |
| range     | string           | The start & end vertices, as well as their strand. <br> Null in dummy vertices. <br> For example, `>4946,>4948` refers to a bubble starting at 4946 and ending at 4948 in the forward strand. |
| children  | string           | A comma seperated string of unsigned integers which are the child vertices <br> Null if the vertex is a leaf.                                                                                 |



## Input
Input GFA

Expect the segments in the input GFA to have unique numeric [segment names](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#s-segment-line).


## Development

To compile povu with debug symbols and with address sanitizer

```
cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER=address  -H. -Bbuild && cmake --build build -- -j 3
```

## Name

The etymology of the name is rooted in profound philosophy 🤔. "Povu," is [Kiswahili](https://en.wikipedia.org/wiki/Swahili_language) for "foam." Foam, by nature, comprises countless bubbles.
