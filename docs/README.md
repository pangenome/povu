# Subcommands

## Table of Contents
 - [Deconstruct](#Deconstruct)
 - [Flubble Tree](#flubble-tree)
   * [flb Format](#flb-format)
 - [Call](#Call)
 - [Input](#Input)

## Overview

povu provides a set of tools that we can use to analyse genomic variation. Currently these are under two main subcommands.

 - `call`:
 - `deconstruct`:

## Deconstruct

**help:** `./bin/povu deconstruct -h`

The `deconstruct` sub-command finds flubbles, reports hairpin inversion boundaries, and generates the flubble forest: a set of flubble trees.
Currently, each flubble tree is stored in its own flubble file in a directory specified by the user using the `-o` CLI argument.

Example 1.
Lines starting with `#` are comments
```
# create an output directory
mkdir results/

# run povu on the yeast dataset
./bin/povu deconstruct -v 0 -i ./test_data/cerevisiae.fa.gz.d1a145e.417fcdf.7493449.smooth.final.gfa -o results
```

Currently hairpin boundaries are printed by `povu deconstruct` at runtime, if none is printed then none was found.


## Flubble Tree

A tree representation of the hierarchy and nesting relationship between flubbles

### flb Format

The flb format is a plain-text file representing flubble tree(s) and ends with the `.flb` extension. It is tab-separated and consists of 3 columns where a dot `.` represents a null value.
Below is a table with the specifics of each column.

| column    | type             | description                                                                                                                                                                                   |
|-----------|------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| vertex id | unsigned numeric | A unique identifier (can be zero & cannot be null).                                                                                                                                           |
| range     | string           | The start & end vertices, as well as their strand. <br> Null in dummy vertices. <br> For example, `>4946,>4948` refers to a flubble starting at 4946 and ending at 4948 in the forward strand. |
| children  | string           | A comma seperated string of unsigned integers which are the child vertices <br> Null if the vertex is a leaf.                                                                                 |


## Call

Using the LPA dataset as an example. Generate a list of references from the P lines.
```
grep '^P' test_data/real/LPA.gfa | cut -f 2 > ~/Data/povu/results/refs.txt
```

Generate the flubble tree for the same dataset

```
povu deconstruct -t 4 -v 2 -i test_data/real/LPA.gfa -o ~/Data/povu/results
```

Generate the VCF files related to the LPA dataset
```
povu call -i test_data/real/LPA.gfa -f ~/Data/povu/results/flb -r ~/Data/povu/results/refs.txt -o ~/Data/povu/results/vcf
```

## Input
Input GFA

Expect the segments in the input GFA to have unique numeric [segment names](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#s-segment-line).

