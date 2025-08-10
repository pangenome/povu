# Subcommands

## Table of Contents
 - [Decompose](#Decompose)
 - [Pangenome Variation Structure Tree](#pangenome-variation-structure-tree)
   * [PVST Format](#pvst-format)
 - [Call](#Call)
    * [The Ref List](#the-ref-list)
 - [Input](#Input)

## Overview

povu provides a set of tools that we can use to analyse genomic variation. Currently these are under two main subcommands.

 - `decompose`: find regions of variation in the graph
 - `call`: generate a VCF file

## Decompose

> **help:** `./bin/povu decompose -h`

The `decompose` sub-command finds flubbles, reports hairpin inversion boundaries, and generates the flubble forest: a set of flubble trees.
Currently, each flubble tree is stored in its own flubble file in a directory specified by the user using the `-o` CLI argument.

Example 1.
Lines starting with `#` are comments
```
# create an output directory
mkdir results/

# run povu on the yeast dataset
./bin/povu decompose -v 0 -i ./test_data/cerevisiae.fa.gz.d1a145e.417fcdf.7493449.smooth.final.gfa -o results
```

Currently hairpin boundaries are printed by `povu decompose` at runtime, if none is printed then none was found.


## Pangenome Variation Structure Tree

A tree representation of the hierarchy and nesting relationship between regions of variation

### PVST format

The pvst format is a plain-text file representing flubble tree(s) and ends with the `.pvst` extension. 


## Call

> **help:** `./bin/povu call -h`


Using the LPA dataset as an example. You have several options to specify reference paths:

1. Generate a list of references from the P lines:
```
grep '^P' test_data/LPA.gfa | cut -f 2 > ~/Data/povu/results/refs.txt
```

2. Or use path prefixes to automatically select paths:
```
# Single prefix
povu call -i test_data/LPA.gfa -f ~/Data/povu/results/flb -P HG -o ~/Data/povu/results/vcf

# Multiple prefixes
povu call -i test_data/LPA.gfa -f ~/Data/povu/results/flb -P HG -P NA -o ~/Data/povu/results/vcf
```

Generate the flubble tree for the same dataset

```
povu decompose -t 4 -v 2 -i test_data/LPA.gfa -o ~/Data/povu/results
```

Generate the VCF files related to the LPA dataset (using reference file):
```
povu call -i test_data/LPA.gfa -f ~/Data/povu/results/flb -r ~/Data/povu/results/refs.txt -o ~/Data/povu/results/vcf
```

Or using path prefixes (multiple prefixes can be specified):
```
povu call -i test_data/LPA.gfa -f ~/Data/povu/results/flb -P HG -P NA -o ~/Data/povu/results/vcf
```

### The Ref List


The ref list is the list of P lines against which we will call variants

## VCF

For the VCF files produced read the [detailed docs here](./vcf.md)


## Input
povu can currently supports GFA version 1.0
Input GFA

Expect the segments in the input GFA to have unique numeric [segment names](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#s-segment-line).

