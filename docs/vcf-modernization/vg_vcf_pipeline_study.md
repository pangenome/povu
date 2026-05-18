# vg / vcfbub / vcfwave VCF pipeline study

Task: `vcf-modern-vg-pipeline-study`
Date: 2026-05-18

## Scope

This note studies the commonly used pangenomics VCF cleanup path:

```text
graph with reference and haplotype paths
  -> vg deconstruct
  -> vcfbub
  -> bcftools norm and/or vcfwave
  -> downstream VCF consumers
```

The goal is not to require povu to clone that pipeline.  The useful part for
povu is the contract boundary: what is best emitted by the core graph-to-VCF
algorithm, what should remain a post-processing normalization option, and which
odd upstream choices are practical enough that downstream users may expect them.

## Sources and versions checked

Checked on 2026-05-18.

| Component | Version / source checked | Notes |
| --- | --- | --- |
| `vg` upstream source | `vgteam/vg` `master` HEAD `d3e3e661451be7477b6ec0b0915c06255fdf7416`, commit date 2026-05-15 | Source files checked: `src/subcommand/deconstruct_main.cpp`, `src/deconstructor.cpp`, and `test/nesting/*.gfa`. |
| Local `vg` binary | `vg version v1.70.0 "Zebedassi"` | Used for tiny command reproductions below. Local help includes experimental `-n/--nested` and `-f/--nested-fasta`; current `master` help source checked above does not show those options, so nested-mode CLI details require confirmation against the exact release povu wants to emulate. |
| `vg deconstruct` docs | `https://github.com/vgteam/vg/wiki/VCF-export-with-vg-deconstruct`, last wiki edit observed from GitHub UI as 2024-06-12 | The wiki documents nested examples and `-n`; it may lag or lead source help in particular releases. |
| `vcfbub` | `pangenome/vcfbub` `main` HEAD and tag `v0.1.2`, commit `77289654b246a4e3422902d04277e258d9fabe9a`, date 2025-06-12 | Source checked: `README.md` and `src/main.rs`. |
| `vcfwave` | `vcflib/vcflib` `master` HEAD `b118a9bfd99b07da9d40d0bd8b3c2bdc4523b568`, commit date 2026-03-20 | Docs checked: `doc/vcfwave.md`, doc last modified by commit `9e8c0192d677bddd68eb2743ff9ec194995e92b7`, date 2025-05-18. Source checked: `src/vcfwave.cpp`. |
| Minigraph-Cactus pipeline docs | `ComparativeGenomicsToolkit/cactus` `master` HEAD `632f051b7cabdc061a458fa165a53f4d8dd6ccdf`; `doc/pangenome.md` last changed at `d8c062bf5975b2a5936568d398f7031759aaac56`, date 2026-02-21 | Used for the pipeline-level contract: raw VCF, vcfbub VCF, optional vcfwave VCF, clipping/filtering, and `bcftools norm`. |
| HPRC pangenome resources | `human-pangenomics/hpp_pangenome_resources` `main` HEAD `74553c422af5521e3d297d6511da0a33fcf3a744` | Used only as context for public Minigraph-Cactus graph outputs. |

Local command reproduction was limited to `vg`; `vcfbub`, `vcfwave`, and
`bcftools` were not installed in this environment.  Their behavior below is
therefore based on upstream docs and source unless marked otherwise.

## Pipeline contracts

### 1. Graph and path contract before VCF export

`vg deconstruct` assumes the graph already contains a reference coordinate
system and non-reference path or GBWT/GBZ haplotype traversals.  The reference
can be chosen by exact path (`-p`) or prefix (`-P`), or inferred from reference
and generic paths.  If every graph path would be selected as reference and no
alternative haplotypes remain, `vg deconstruct` errors rather than emitting an
empty or all-reference callset.

Important contract points:

- The VCF reference coordinate system is selected before site decomposition.
  Allele spelling is not reference-neutral.
- Path metadata matters.  `vg deconstruct` parses sample and haplotype
  structure from path names or GBWT metadata, and it infers phasing from
  haplotypes when available.
- GBZ input supplies both graph and GBWT haplotype data.  When GBWT/GBZ
  haplotypes are used, the current source disables the context-Jaccard path
  remapping window because step-level path positions are not available for those
  traversals.
- Node ID translation can be applied to VCF snarl IDs and traversal fields with
  `-T` or, for GBZ, `-O`.

For povu, this argues for a clear VCF-capable graph contract distinct from a
topology-only flubble contract: labels, references, path/sample metadata, and
stable traversal IDs are part of VCF emission even though they are not needed
to detect flubbles.

### 2. `vg deconstruct`: decompose graph snarls into VCF records

`vg deconstruct` reports graph snarls as VCF sites relative to a reference
path.  In current source, it loads or computes snarls, then deconstructs them
with a `Deconstructor` configured for:

- top-level snarls only by default;
- all snarls with `-a/--all-snarls`;
- optional traversal untangling with `-u`;
- conflict handling with `-K` and `-S`;
- optional traversal clustering with `-L`;
- optional star allele handling with `-R`, only when all/nested snarls are
  requested in the current source.

The normal VCF header includes `GT`, `AC`, `AF`, `NS`, `AN`, and either:

- `AT`: graph traversal strings for REF and ALT alleles; or
- `UT`: traversal strings annotated with reference start/end positions when
  `-u` is used.

When nested/all-snarl mode is used, the source emits:

- `LV`: level in the snarl tree, where `0` is top-level;
- `PS`: ID of the parent snarl;
- `RC`, `RS`, `RD`: reference contig and reference interval of the top-level
  containing site.

Local `vg v1.70.0` also exposes experimental `-n/--nested`, whose output adds
context tags such as `PA`, `PL`, `PR`, and `RL` in the tested tiny examples.
Those fields were not present in the `master` source help checked on
2026-05-18, so they should not be specified for Rust output without a fixture
task against the intended `vg` release.

#### Local nested reproduction

Command, using upstream `test/nesting/nested_snp_in_del.gfa` and converting GFA
through `vg convert` because process substitution has no `.gfa` extension for
VPKG format inference:

```bash
vg deconstruct -p x -n \
  <(curl -fsSL https://raw.githubusercontent.com/vgteam/vg/master/test/nesting/nested_snp_in_del.gfa \
    | vg convert -g -)
```

Local `vg v1.70.0` emitted a parent top-level deletion/MNP site and an inner
SNP site:

```text
x  1  >1>6  CATG  CAAG,C  ...  LV=0                       GT  1|2
x  3  >2>5  T     A       ...  LV=1;PS=>1>6               GT  1|.
```

With `-R`, the inner site used a star allele for the haplotype that traverses
the parent but not the child:

```text
x  3  >2>5  T  A,*  ...  AT=>2>3>5,>2>4>5,.;LV=1;PS=>1>6  GT  1|2
```

The practical choice is clear: upstream supports both "missing because this
haplotype does not traverse this nested site" and "explicit spanning deletion
with `*`" semantics, depending on options.

#### `-a` versus `-n` in the local binary

For the same tiny deletion, local `vg deconstruct -a` emitted the nested child
with `LV/PS`, but not the richer experimental nested context tags produced by
`-n`.  This makes `-a` look like the stable flattenable raw mode, while `-n`
is a richer nested representation with extra context fields and possible
off-reference child coordinates.  Because the 2026-05-18 `master` source help
does not show `-n`, the exact status of this split is an open confirmation
item.

#### Off-reference nested coordinates

Local `vg v1.70.0` on upstream `nested_snp_in_ins2.gfa` with `-n` emitted the
inner SNP on an inserted haplotype contig rather than the top-level reference
contig:

```text
a#1#y0  3  >2>5  A  T  ...  LV=1;PS=>1>6  GT  0|1  0|.
x       1  >1>6  C  CAAG,CATG  ...  LV=0  GT  1|2  1|0
```

This is surprising but practical: a child bubble inside an insertion has no
ordinary coordinate on the outer reference.  `vg` can expose it on the path
that actually contains the nested sequence.  Downstream VCF consumers that
expect all records to share the primary reference contig may dislike this;
graph-aware consumers get a coordinate that preserves the nested allele.

### 3. `vcfbub`: pop large or unwanted parent bubbles

`vcfbub` is intentionally narrow.  It reads `vg deconstruct` VCFs containing
snarl nesting tags (`LV` and `PS`) and filters records by:

- maximum snarl level (`-l/--max-level`);
- maximum allele length (`-a/--max-allele-length`);
- maximum reference allele length (`-r/--max-ref-length`).

The non-obvious behavior is its "pop parent, keep child" rule.  In source,
`vcfbub` first marks each record as keep/drop.  A record is dropped if its
`LV` is greater than the maximum level, or if its max allele length or ref
allele length exceeds the configured maximum.  Records dropped by size are
stored as popped bubbles.  On the second pass, a child record is kept if its
parent `PS` is one of the popped bubbles and the child itself was not popped.

So `vcfbub -l 0 -r 10000 input.vcf` is not merely "keep `LV=0` and ref length
under 10 kb."  It keeps top-level sites unless large, and if a large parent is
popped, it can keep nested children even though those children have `LV > 0`.

The upstream README motivates this with large inversions, multi-megabase
structural variants, and segmental duplications: the giant parent may be poor
input for imputation panels or haplotype panels, while the smaller nested
events inside it are useful.

### 4. `bcftools norm`: reference-based left normalization

Minigraph-Cactus documents that, by default since Cactus v2.8.2, non-raw VCF
output is normalized with `bcftools norm -f`.  This is a post-deconstruct,
reference-FASTA-dependent step that left-aligns and normalizes indels.  It is
not a graph decomposition algorithm; it rewrites VCF alleles after graph sites
have already been selected.

For povu, this is an output-normalization option or recommended downstream
command, not part of the core flubble algorithm.

### 5. `vcfwave`: align complex alleles into simpler VCF records

`vcfwave` is in `vcflib`.  Its docs describe it as reducing complex alleles by
pairwise alignment against the reference allele with BiWFA/WFA.  Its source and
docs show this contract:

- Read existing VCF records.
- Skip records that are already simple enough or exceed `--max-length`.
- Align each ALT allele to REF.
- Split complex alleles into simpler SNP/MNP/INS/DEL records.
- Add `TYPE`, `LEN`, and an origin tag, default `ORIGIN`.
- Preserve and adjust genotypes; update `AC`, `AF`, `AN`, and `AT` where they
  are present.
- Detect long inversions and emit `INV=YES` rather than decomposing them into
  many primitive records.
- For deletions, overlapping SNP/MNP calls in samples carrying the deletion can
  become missing/haploid at the overlapping site.

The documented workflow is:

```bash
vcfwave -L 1000 raw-or-bubbed.vcf \
  | bcftools norm -m- ...
```

followed by optional `vcfcreatemulti` to reassemble compatible records into
multi-allelic VCF records.  The ordering matters: `vcfwave` is an allele-level
realignment/decomposition pass after graph-site selection and, in the Cactus
pipeline, after or alongside the `vcfbub` path.

### 6. Minigraph-Cactus pipeline-level contract

The Cactus pangenome docs describe two VCFs for selected graph types:

- a raw VCF containing nested variants, indicated by `LV` and `PS`;
- a vcfbub-processed VCF that removes nested sites and sites above a size
  threshold, defaulting in the docs to 100 kb.

Cactus also documents optional `--vcfwave` output and default `bcftools norm`
normalization for non-raw VCFs since v2.8.2.  In other words, the mature
pipeline keeps graph-aware raw output and consumer-friendly normalized output
separate.

The graph construction side also matters.  Cactus clipping/filtering can remove
unaligned sequence stretches, dangling non-reference tips, or low-haplotype
coverage nodes before VCF export.  Those graph-level choices are upstream of
`vg deconstruct` and change the set of snarls available to deconstruct.

## Behaviors povu should emulate, avoid, or expose

| Area | Recommendation for povu | Upstream basis |
| --- | --- | --- |
| Raw nested VCF tags | Emulate stable `LV` and `PS` for nested flubbles.  Add parent/top-level context fields only after choosing a concrete schema and release target. | `vg deconstruct -a` and `vcfbub` use `LV/PS` as the minimal nesting contract.  Local `-n` adds richer tags, but source/doc status is inconsistent. |
| Allele traversal provenance | Emit a graph traversal field by default.  Prefer an `AT`-compatible path-string field for graph-aware users, and consider an optional untangled coordinate field similar to `UT`. | `vg deconstruct` always emits `AT` unless `-u` requests untangled traversal coordinates.  `vcfwave` can preserve/update `AT`. |
| REF/ALT anchoring | Use VCF-compliant anchored alleles for indels, with deterministic reference coordinate rules.  Keep graph traversal IDs so the allele can be traced back even after sequence normalization. | `vg deconstruct` anchors to the reference traversal through a snarl.  `vcfwave` may later shift/split sequence alleles while preserving origin fields. |
| Nested child genotypes | Expose an option for nested spanning behavior: missing genotype for non-traversal versus `*` allele for spanning parent deletions.  Do not silently mix policies. | Local `vg v1.70.0` emits `1|.` without `-R` and `1|2` with ALT `*` under `-R`. |
| Pop-large-parent behavior | Implement or expose a vcfbub-like post-pass, not as part of core flubble detection.  It needs parent IDs, levels, and allele/ref lengths. | `vcfbub`'s useful behavior depends on filtering parent records while rescuing children of popped parents. |
| Complex allele decomposition | Treat vcfwave-like WFA decomposition as downstream normalization.  Do not make core povu VCF output depend on WFA unless a later task explicitly scopes it. | `vcfwave` works on VCF REF/ALT strings after graph-site selection and rewrites records/genotypes. |
| Repeats and redundancy | Keep repetitive-sequence untangling and large-bubble popping downstream of the core flubble algorithm.  Provide thresholds/options for export profiles. | Cactus uses graph clipping/filtering, vcfbub thresholds, `bcftools norm`, and optional `vcfwave` as post/core-adjacent pipeline controls, not as snarl detection semantics. |
| Source/version metadata | Include `##source=povu` plus optional `##povuCommand`, `##povuVersion`, and normalization/provenance tags. | The upstream pipeline relies on tool-specific fields (`AT`, `ORIGIN`, `TYPE`, `LEN`, `INV`) to keep rewritten VCFs auditable. |

## Edge cases and practical surprises

### Nested insertion child records may not be on the primary reference contig

When a child site exists inside an insertion, local `vg -n` can place the child
record on the inserted path, not on the outer reference contig.  This is graph
faithful but VCF-consumer surprising.

Rust implication: decide whether povu's first practical VCF spec allows
off-reference `CHROM` values for nested children.  If not, nested insertion
children probably need either no separate VCF record, an `LV/PS` record with
top-level coordinates and explicit graph traversal, or a downstream fixture
task to test what users actually need.

### `-a`, `-n`, and `-R` are not one semantic switch

Local `vg v1.70.0` distinguishes:

- top-level only: no `LV/PS`;
- `-a`: all snarls with `LV/PS`;
- `-n`: richer experimental nested context;
- `-R`: explicit `*` alleles for haplotypes spanning parent but not child.

Current `master` source help checked on 2026-05-18 does not show `-n`.  This
must be confirmed before copying the richer nested fields.

### `vcfbub` can output records above the requested level

Because children of popped parents are rescued, output from `vcfbub -l 0 -r N`
can include `LV > 0` records.  That is intentional.  It makes the output more
useful in complex SVs, but consumers must not interpret `--max-level 0` as an
absolute postcondition.

Rust implication: if povu implements a "flat" export profile, name the
postcondition precisely.  "Top-level only" and "vcfbub-style popped output" are
different modes.

### Missing genotypes can be an allele-overlap statement, not missing data

`vg` nested output can use `.` for a haplotype that does not traverse a nested
site.  `vcfwave` can also convert overlapping SNP/MNP genotypes to missing for
samples carrying a deletion.  These are not ordinary low-confidence no-calls;
they are structural overlap consequences.

Rust implication: the VCF model should distinguish:

- no sample data;
- no traversal through this nested site;
- allele absent because a deletion spans the site;
- explicit `*` spanning allele.

The text VCF may collapse some of these to `.`, but the internal Rust model
should not.

### Allele clustering is lossy but useful

`vg deconstruct -L F` clusters traversals by graph-handle Jaccard similarity
when similarity is at least `F`.  With clustering enabled, source can emit
`TS` and `TL` fields describing sample path similarity and length difference.
This is not normalization in the VCF-left-alignment sense; it is a lossy graph
allele grouping step to keep nearly equivalent traversals together.

Rust implication: keep any clustering option explicit and visibly lossy.
The default should preserve distinct flubble traversals unless the output spec
has a separate "collapsed alleles" profile.

### Inversions are not always decomposed

`vcfwave` detects inversions and can emit `INV=YES` with `TYPE=mnp` instead of
splitting an inversion into many primitive edits.  That is practical: a long
reverse-complement allele is more interpretable as one inversion than as a
field of substitutions after pairwise alignment.

Rust implication: povu's current `SUBR`/hairpin convention should be mapped
deliberately to either graph-native inversion records or downstream `vcfwave`
normalization.  Do not expect a generic primitive decomposition pass to preserve
inversion semantics automatically.

## Repetitive-sequence untangling and normalization

The upstream pipeline treats repetitive or redundant sequence as a layered
problem:

1. Graph construction may clip unaligned sequence and dangling non-reference
   tips, and can filter low-haplotype coverage subgraphs for mapping indexes.
2. `vg deconstruct` exports graph snarls and can optionally untangle traversal
   positions with `-u`, cluster similar graph traversals with `-L`, and carry
   path-level provenance in `AT` or `UT`.
3. `vcfbub` removes or pops very large parent bubbles while preserving smaller
   child variation inside them.
4. `bcftools norm` left-aligns and normalizes VCF alleles against a reference
   FASTA.
5. `vcfwave` uses WFA to split complex REF/ALT strings into simpler records and
   updates genotypes.

For povu, the core flubble algorithm should remain a faithful graph topology
and hierarchy detector.  Repetitive-sequence untangling should be downstream of
that core for three reasons:

- It is representation-policy heavy.  A large nested repeat may be a valid
  flubble even if a consumer wants it popped out of a VCF panel.
- It depends on thresholds and use case.  Cactus documents graph filtering,
  vcfbub size thresholds, and optional vcfwave because mapping, imputation, and
  graph visualization need different trade-offs.
- It rewrites records.  `bcftools norm` and `vcfwave` can change POS, REF, ALT,
  record count, and genotype columns after the core graph site has already been
  selected.

Recommended split:

- Core Rust/povu output: deterministic flubble records with stable hierarchy
  (`LV/PS`), traversal provenance (`AT` or equivalent), anchored REF/ALT, and
  clear genotype semantics.
- Export profile `raw-graph`: preserve nested hierarchy and complex alleles.
- Export profile `popped`: vcfbub-like parent popping by level/ref length/max
  allele length.
- Export profile `normalized`: run or recommend `bcftools norm` and optionally
  `vcfwave`, with explicit provenance that records were rewritten.

## Open questions for Rust implementation

1. Which `vg` release is the compatibility target: local `v1.70.0`, current
   `master`, Cactus-bundled `vg`, or an HPRC-pinned container version?  This
   matters for `-n`, `-f`, `PA/PL/PR/RL`, and star-allele behavior.

2. Should povu implement only the stable `LV/PS/AT` contract first, or also
   model richer nested-context tags from experimental `vg -n`?  The latter
   needs fixtures for insertion-inside-reference, deletion-with-child,
   inversion, cyclic reference, and off-reference nested child coordinates.

3. What is povu's exact genotype semantics for nested sites?  The Rust model
   should distinguish missing data, missing because no traversal reaches the
   child, and explicit VCF `*` spanning deletion alleles.

4. Should `AT` be byte-for-byte compatible with `vg` traversal strings, or
   should povu keep its existing traversal/provenance field names and provide a
   compatibility alias?  Compatibility helps `vcfwave` and graph-aware tooling,
   but field naming also needs to align with existing `ES` and `VARTYPE`
   conventions in povu.

5. How should a nested child inside an insertion be represented if povu does
   not want off-reference `CHROM` records?  Candidate policies are: suppress
   the child from VCF, emit it with top-level reference coordinates and graph
   provenance, or permit inserted-path contigs in raw graph mode only.

6. Should vcfbub-like popping be implemented in Rust as a library pass over
   povu VCF records, or left to external `vcfbub` compatibility?  A Rust pass
   would need parent snarl IDs, level, reference length, max allele length, and
   a two-pass child-rescue algorithm.

7. Should `vcfwave` be an optional external normalization step in documentation,
   or should a future Rust task implement WFA-based decomposition?  The latter
   is a substantial dependency and should not block the core VCF modernization
   spec.

8. How will normalization provenance be recorded?  If povu writes raw graph
   records and a later command rewrites them, fields like `ORIGIN`, `TYPE`,
   `LEN`, and original traversal tags are needed to audit the transformation.

9. What thresholds should be defaults versus named profiles?  Upstream examples
   use `vcfbub` parent popping around 10 kb in the README and Cactus docs
   describe a 100 kb default for non-raw VCF output.  These are workflow
   defaults, not mathematical properties.

10. Which downstream consumers are target-critical: bcftools-only tools,
    imputation panels, graph-aware validators, or pangenome browsers?  The best
    choice for nested records and off-reference child coordinates depends on
    this answer.

## Confirmation tasks to file if needed

The following behaviors are useful but should be fixture-confirmed before the
Rust output spec treats them as requirements:

- Run the exact target `vg` version on upstream `test/nesting/*.gfa` for
  top-level, `-a`, `-n`, `-n -R`, and `-u` outputs.
- Run `vcfbub` on a tiny nested VCF where a large `LV=0` parent contains a
  small `LV=1` child, confirming the child-rescue behavior and exact header
  preservation.
- Run `vcfwave` on a multi-ALT repetitive deletion/SNP example to confirm
  genotype nullification around deletion overlap and preservation/update of
  `AT`, `AC`, `AF`, and `AN`.
- Run `bcftools norm -f` after `vg deconstruct` on an indel in a repeat to
  show how POS/REF/ALT shift while `AT` provenance remains graph-specific.

## Source links

- vg deconstruct wiki: https://github.com/vgteam/vg/wiki/VCF-export-with-vg-deconstruct
- vg deconstruct source: https://github.com/vgteam/vg/blob/master/src/subcommand/deconstruct_main.cpp
- vg deconstructor source: https://github.com/vgteam/vg/blob/master/src/deconstructor.cpp
- vg nesting tests: https://github.com/vgteam/vg/tree/master/test/nesting
- vcfbub README: https://github.com/pangenome/vcfbub/blob/main/README.md
- vcfbub source: https://github.com/pangenome/vcfbub/blob/main/src/main.rs
- vcfwave docs: https://github.com/vcflib/vcflib/blob/master/doc/vcfwave.md
- vcfwave source: https://github.com/vcflib/vcflib/blob/master/src/vcfwave.cpp
- Minigraph-Cactus pangenome docs: https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md
- HPRC pangenome resources: https://github.com/human-pangenomics/hpp_pangenome_resources
