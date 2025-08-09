# VCF Format

Output in [VCF 4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf) it represents variation as follows

## Deletions
When we finds a deletion of one or more bases, we include at least one “anchor” base so that neither REF nor ALT becomes empty

 - **REF** holds the anchor base plus the deleted bases.
 - **ALT** holds only the anchor base.
 - **POS** is the coordinate of that anchor base.

Suppose the reference sequence is `…A C G T…` at positions 99–102, and you want to delete C G (positions 100–101). In VCF we write:


```
#CHROM POS   ID  REF  ALT
chr1   100  .   ACG  A

```

## Insertions

Just like deletions, every insertion must include an anchor base so that neither REF nor ALT is empty.

 - **REF** holds the anchor base only.
 - **ALT** holds the anchor base followed by the inserted sequence.
 - **POS** is the coordinate of the anchor base.

For example, suppose the reference is

```
... A   T   G ...

   100 101 102
```

and you want to insert “CCC” between A (pos 100) and T (pos 101). You anchor on A (pos 100):

```
#CHROM  POS   ID  REF  ALT
chr1    100   .   A    ACCC
```

This means “at position 100, after the A, insert CCC.”



## Substitution

A simple substitution replaces one (or more) bases with another sequence. It does not require special anchoring beyond the usual rules:

 - **REF** is the reference sequence being replaced.
 - **ALT** is the new sequence.
 - **POS** is the coordinate of the first base of the REF string.

For a single-base substitution—say C → G at position 150—you’d write:

```
#CHROM  POS   ID  REF  ALT
chr1    150   .   C    G
```

For a multi-base substitution—e.g. replacing “CGT” at positions 200–202 with “TAA”—you do:



```
#CHROM  POS   ID  REF   ALT
chr1    200   .   CGT   TAA

```

Here POS = 200, REF = CGT, ALT = TAA, indicating “at 200–202, replace CGT with TAA.”

## Qual

Fixed at 60

## Filter

Fixed at PASS

## Info field

| field name | VCF description                                       |                                                                  |
|------------|-------------------------------------------------------|------------------------------------------------------------------|
| `AN`       | Total number of alleles in called genotypes           | Total number of called alleles (ref + alt), across all samples   |
| `AC`       | Total number of alternate alleles in called genotypes | Number of alternate alleles only, across all samples             |
| `AF`       |                                                       | AF=AC/AN                                                         |
| `AT`       |                                                       | not always (but tries to be) minimal walks for ref and alt calls |
| `NS`       |                                                       | how many samples at this variant site have genotype data         |
|            |                                                       |                                                                  |

## Genotyping

Genotyping is based on [PanSN](https://github.com/pangenome/PanSN-spec)


We use `P` lines to determine PanSN genotype information

If a ref does not go through a rov then it is marked as `.`


Either the label of a `P` line in the GFA will contain the following information

sample_name := string
haplotype_id := number
contig_or_scaffold_name := string


If a sample name does not follow the PanSN spec then the entire P line will be considered to be a sample name 



The [ref list](./README.md#the-ref-list) should contain the exact P lines to
call variants against
