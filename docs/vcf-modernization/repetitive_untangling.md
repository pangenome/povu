# Repetitive Untangling and Normalization Strategy

Task: `vcf-modern-untangling-design`

Date: 2026-05-18

## Scope

This document designs the downstream strategy for repetitive-sequence
untangling and normalization after povu has produced the raw graph VCF defined
in `docs/vcf-modernization/rust_vcf_output_spec.md`.

The input to this design is not arbitrary VCF. It is povu's raw graph output:

- one Lean-aligned semantic `VariantCall` per record;
- graph traversal provenance in `AT`;
- flubble hierarchy in `ES` and `LV`;
- hairpin reverse substitutions as `SUBR` without flubble hierarchy fields;
- anchored but not left-normalized `REF`/`ALT`;
- `TANGLED=T` for graph-faithful records that should not be mistaken for
  primitive normalized variants.

The goal is to expose complex and repetitive structure without making the core
flubble/hairpin algorithm responsible for consumer-specific VCF rewriting.

## Decision Summary

povu should keep repetitive-sequence untangling downstream of the core
algorithm. The core algorithm should emit faithful graph sites and hierarchy.
Downstream profiles may then pop large parents, left-normalize alleles,
decompose complex alleles, or adapt tags for external tools.

The safe split is:

| Layer | Responsibility | In core? |
| --- | --- | --- |
| Raw graph emission | Emit deterministic semantic records with `AT`, `ES`/`LV`, anchored alleles, genotypes, and `TANGLED`. | Yes |
| Profile filtering | Select records by level or policy without rewriting sequence alleles. | Downstream |
| vcfbub-like popping | Drop large parent records and rescue eligible children of popped parents. | Downstream |
| Reference normalization | Left-align and normalize `POS`/`REF`/`ALT` against a FASTA. | Downstream |
| vcfwave-like decomposition | Split complex alleles with WFA/BiWFA-style alignment and update genotypes/counts/provenance. | Downstream |
| External compatibility tags | Add `PS`, origin tags, or tool-specific metadata needed by vg/vcfbub/vcfwave-style pipelines. | Downstream or adapter |

This means the default povu claim remains narrow and auditable: raw records are
graph-faithful and semantically checked; downstream records are useful VCF
rewrites with explicit provenance, not additional facts proved by the flubble
detector.

## Why Untangling Is Downstream

Repetitive-sequence untangling is a representation policy, not a property of
the flubble boundary itself. A large repeat-spanning parent may be a valid
flubble and still be poor input for imputation, panel construction, or a
bcftools-only workflow. Conversely, preserving that parent may be exactly what
a graph-aware validator or browser needs.

Untangling also depends on workflow-specific thresholds. A vcfbub-like profile
needs maximum level, maximum allele length, maximum reference length, and a
child-rescue policy. A vcfwave-like profile needs a maximum alignment length,
rules for inversions, and overlap behavior for genotypes. These are not stable
mathematical invariants of povu's core graph algorithm.

Most importantly, normalization rewrites the VCF. `bcftools norm`-style
normalization can change `POS`, `REF`, and `ALT`. vcfwave-like decomposition can
change the record count, split multi-allelic records, alter genotype allele
indexes, update `AC`/`AF`/`AN`, and mark overlapping calls missing. Those
rewrites must not be hidden inside the formally checked emission boundary.

The Lean/Rust guardrail is strongest at the semantic raw-record boundary:
verified flubble/hairpin sources, well-formed `VariantCall`s, deterministic
ordering, and valid graph-derived VCF fields. Downstream untangling should be a
separately tested transformation over those records.

## Core Non-Goals

The povu core algorithm and default Rust emitter should not:

- left-shift indels across tandem repeats;
- trim shared prefixes or suffixes beyond the semantic allele construction
  supplied by the raw call;
- decompose one complex flubble into several SNP/indel records;
- run WFA/BiWFA or require a pairwise alignment library for VCF emission;
- cluster nearly equivalent graph traversals unless an explicit lossy profile
  asks for it;
- fabricate coordinates for a nested child that has no coordinate on the
  selected reference path;
- silently convert nested non-traversal into VCF `*` alleles;
- discard large parent sites merely because they are inconvenient for a
  consumer;
- assume external `vcfbub` can consume default povu VCF before a compatibility
  adapter emits the fields it requires;
- claim normalized primitive VCF semantics for records still marked
  `TANGLED=T`.

The core may expose enough metadata for downstream tools to make those choices:
source ids, parent ids when available, levels, traversal strings, allele
lengths, reference spans, genotype semantics, and stable origin identifiers.

## Staged Pipeline

### Stage 0: Raw Graph VCF

The Rust emitter writes the `raw-graph` profile from the output spec.

Required properties:

- preserve `AT` for REF and every ALT;
- preserve `ES` and `LV` for flubble records;
- keep `SUBR` as graph-native hairpin output without `ES`/`LV`;
- emit anchored `REF`/`ALT`, but do not left-normalize;
- mark complex graph-faithful sites with `TANGLED=T`;
- fail with a structured error rather than inventing coordinates for
  unsupported nested children.

This stage is the only stage that should be considered part of the initial
Lean-aligned core VCF contract.

### Stage 1: Annotation and Compatibility Preparation

A downstream profile runner should first derive explicit transformation
metadata from raw records without changing sequence alleles.

Recommended derived fields or sidecar values:

- stable raw record id, for example `OID` or sidecar `raw_id`;
- parent id for nested records, derived from the semantic hierarchy;
- top-level enclosing id when available;
- `raw_pos`, `raw_ref`, `raw_alt`, and raw traversal strings in a sidecar
  manifest if the VCF fields may later be rewritten;
- maximum allele length and reference allele length;
- profile name and command/configuration.

For external `vcfbub` compatibility, this is where povu can add a `PS` alias
from the semantic parent id. The default raw VCF uses `ES` plus `LV`, so a
compatibility adapter is needed before assuming ordinary `vcfbub` can process
povu records unchanged.

### Stage 2: Simple Filters

The first downstream filters should be explicit about whether they rewrite
records.

Recommended profiles:

- `raw-graph`: no filtering or rewriting.
- `top-level-only`: keep `LV=0` flubble records and all `SUBR` records. This is
  a simple filter and is not vcfbub popping.
- `nested-only` or `audit-nested`: optional diagnostic filter for inspecting
  child records and off-reference cases.

These filters are useful for debugging and consumer selection, but they should
not be described as normalization.

### Stage 3: vcfbub-Like Popping

The `popped` profile should be a two-pass downstream transformation inspired by
vcfbub:

1. Mark each raw record as kept or popped according to configuration:
   `max_level`, `max_ref_length`, `max_allele_length`, and any future
   `max_record_span` threshold.
2. Remember records popped by size as popped parents.
3. Rescue eligible child records whose parent was popped, provided the child
   itself does not violate the configured thresholds.
4. Emit a record set plus provenance describing which parents were popped and
   which child records were rescued.

The important user-facing rule is that `popped` is not the same as
`top-level-only`. With child rescue, a `popped` output configured with
`max_level=0` may still contain `LV>0` records if their large parent was popped.
The profile name and metadata must make that clear.

Implementation can start as an internal Rust pass over semantic records or an
adapter that emits vg-style `LV`/`PS` tags and calls external `vcfbub`. The
internal pass is preferable for preserving povu-specific provenance, but it
needs fixtures before implementation.

### Stage 4: Reference-Based Left Normalization

Left normalization should be a separate `left-normalized` or `normalized`
profile that requires a reference FASTA matching the VCF contigs.

Rules:

- run only after raw graph emission and optional popping/filtering;
- reject or quarantine records whose `CHROM` cannot be matched to the FASTA;
- preserve raw origin in an `ORIGIN`-style field or sidecar manifest;
- treat changes to `POS`, `REF`, and `ALT` as expected rewrite output, not as
  raw povu semantics;
- keep `AT` or raw traversal provenance attached when the normalizer preserves
  it, and otherwise store it in the sidecar manifest.

This stage is where `bcftools norm -f`-like behavior belongs. It is not part of
flubble detection, hairpin detection, or semantic raw emission.

### Stage 5: vcfwave-Like Allele Decomposition

A `decomposed` profile may run a vcfwave-like pass after raw emission and,
usually, after popping/filtering.

Responsibilities:

- align each ALT allele against REF for records selected by policy;
- split complex records into simpler SNP, MNP, insertion, and deletion records;
- update genotype allele indexes;
- recompute or update `AC`, `AF`, `AN`, and `NS`;
- preserve or update `AT` when possible;
- attach an origin field linking every emitted primitive record back to the raw
  record and allele;
- keep long inversions or `SUBR` records graph-native unless an explicit
  inversion-normalization profile is chosen.

The pass must have a maximum allele length and a pass-through policy. If a
record is too long or too ambiguous to decompose safely, the output should keep
the raw record with `TANGLED=T` or place it in a rejected/sidecar stream. It
should not emit primitive-looking records that imply a decomposition was
successful when it was not.

### Stage 6: Final Audit

Every non-raw profile should emit enough metadata to answer:

- which raw record and allele produced this output record;
- whether `POS`, `REF`, `ALT`, record count, or genotypes changed;
- which thresholds and tools were used;
- whether any records were dropped, popped, rescued, left unchanged, or
  rejected;
- whether graph traversal provenance remains in VCF fields or only in a
  sidecar manifest.

For command-line UX, the pipeline should prefer named profiles over many
loosely coupled booleans:

| Profile | Intended consumer | Rewrites VCF alleles? |
| --- | --- | --- |
| `raw-graph` | graph-aware analysis, auditing, conformance | No |
| `top-level-only` | simple inspection of top-level calls | No |
| `popped` | panel/imputation-friendly nested output | No sequence rewrite, but drops/rescues records |
| `left-normalized` | bcftools-compatible indel normalization | Yes |
| `decomposed` | primitive-ish SNP/indel consumers | Yes |
| `audit-all` | development and validation | No, plus sidecar summaries |

## Handling Difficult Cases Without Misleading VCF

### Tandem Repeats

Raw output should anchor indels according to the semantic call. A downstream
reference-normalization step may left-shift them across the repeat if a matching
FASTA is supplied. Both coordinates must be auditable because either spelling
can be valid VCF while only the raw spelling directly reflects povu's graph
source.

### Large Parent with Small Children

A large repeat or segmental-duplication parent may be unusable for some
consumers. The `popped` profile may remove the parent and rescue smaller child
records. The raw output must remain available because dropping the parent
changes the interpretation of the record set.

### Nested Child Inside an Insertion

If a child site has no coordinate on the selected reference path, core povu
must not invent a top-level `CHROM`/`POS`. Acceptable policies are:

- reject the child for default VCF output with a structured unsupported error;
- omit the child from default VCF while preserving it in graph-aware output;
- permit an explicit raw mode with inserted-path contigs;
- add a future compatibility mode once fixtures define the exact coordinate
  contract.

Any of those is preferable to silently placing the child at a misleading
coordinate.

### Missing, Non-Traversal, and Star Alleles

The raw semantic model should distinguish:

- no sample data;
- no traversal through a nested child;
- allele absent because a deletion spans the child;
- explicit VCF `*` spanning deletion allele.

The first raw output contract serializes nested non-traversal as `.`. A future
`star-spanning` profile may introduce `*`, but it should be explicit and
fixture-backed because it changes ALT indexing and genotype interpretation.

### Hairpins and Inversions

`SUBR` should remain a graph-native inversion-like record in raw output.
vcfwave-like decomposition should not blindly turn long reverse-complement
alleles into a field of substitutions. If a downstream profile emits
`INV=YES`, symbolic SV tags, or primitive decomposition, it must retain the raw
`SUBR` origin.

### Lossy Traversal Clustering

Traversal clustering can reduce redundant alleles in repeats, but it is lossy.
The default should preserve distinct traversals. If povu later supports a
clustered profile, it needs explicit similarity thresholds, cluster provenance,
and tests showing that allele counts and genotypes remain interpretable.

## Provenance Contract

Downstream profiles should treat raw povu records as immutable source facts.
Rewritten records need origin metadata.

Minimum provenance for rewritten VCF:

- `ORIGIN` or equivalent source record id;
- original allele index when a raw ALT is decomposed;
- original `AT` traversal string or sidecar pointer;
- profile name and tool/version;
- flags such as `POPPED_PARENT`, `RESCUED_CHILD`, `LEFT_NORMALIZED`,
  `DECOMPOSED`, or `PASSTHROUGH`;
- old-to-new record mapping in a machine-readable sidecar when one raw record
  becomes many output records.

If VCF INFO fields become too crowded, a sidecar JSONL/TSV manifest is safer
than overloading `AT`, `ES`, or `LV` with rewrite state.

## Fixture Classes Needed Before Implementation

Before implementing the downstream profiles, the conformance corpus should add
fixtures for these classes. The expected result for each fixture should include
raw output and, where applicable, one downstream profile output.

| Fixture class | What it protects |
| --- | --- |
| Tandem-repeat insertion/deletion | Raw anchored alleles stay unchanged in `raw-graph`; left normalization shifts only in the normalization profile. |
| Homopolymer repeat with multiple valid VCF spellings | Normalized output is FASTA-dependent and provenance preserves raw coordinates. |
| Large `LV=0` parent with small `LV=1` child | vcfbub-like popping removes the parent and rescues eligible children; `top-level-only` does something different. |
| Child exceeding popping thresholds | Child rescue does not override the child's own size/level rejection. |
| Multi-ALT tangled repeat record | vcfwave-like decomposition splits or passes through according to max-length policy and updates genotype allele indexes. |
| Deletion overlapping SNP/MNP | Downstream decomposition handles overlapping genotypes explicitly instead of treating structural non-traversal as ordinary no-call. |
| Nested child inside insertion | Default policy rejects/omits or explicitly coordinates the child; no fabricated primary-reference coordinate. |
| Hairpin/SUBR inside repetitive sequence | Raw `SUBR` remains graph-native; optional inversion/decomposition output preserves `SUBR` origin. |
| Off-reference inserted-path contig | If an explicit raw mode allows inserted-path `CHROM`, contig/header and consumer warnings are deterministic. |
| External vcfbub compatibility tags | Adapter emits `LV`/`PS` or the chosen equivalent correctly from povu hierarchy. |
| Reference FASTA mismatch | left-normalization rejects missing/mismatched contigs instead of producing shifted alleles against the wrong reference. |
| WFA max-length pass-through | Records above decomposition length limits remain marked `TANGLED=T` or are quarantined, not falsely simplified. |
| Multi-sample phased genotypes | Popping and decomposition preserve phase separators, missing phases, `NS`, `AN`, `AC`, and `AF`. |
| Provenance one-to-many mapping | One raw record split into multiple records has complete origin mapping for each emitted record. |

These fixtures should be created before any implementation task claims that the
downstream profiles are ready for users. They can be smaller than biological
benchmarks; each should isolate one policy decision.

## Filed Follow-Up Work

This design implies downstream work rather than immediate implementation in the
core emitter. Follow-up task `add-downstream-repetitive` was filed to turn the
fixture classes above into concrete conformance assets before synthesis.

Remaining follow-up decisions for synthesis and later implementation:

1. Decide whether the first implementation uses an internal Rust popping pass,
   an external `vcfbub` adapter, or both.
2. Define the provenance sidecar format for record rewrites.
3. Add an integration/synthesis check that verifies final user-facing claims
   distinguish raw graph VCF from popped, normalized, and decomposed VCF.

None of this should block the core Rust VCF emitter from producing raw graph
records. It should block any claim that povu already emits consumer-normalized
primitive VCF in repetitive regions.

## Synthesis Guidance

For `vcf-modern-synthesis-check`, the safe final claim is:

> povu can emit or target a Lean-aligned raw graph VCF contract for supported
> flubble/hairpin calls, and it has a designed downstream path for repetitive
> popping and normalization.

The unsafe claims are:

- povu core solves repetitive-region normalization;
- povu raw output is already left-normalized across repeats;
- povu raw `TANGLED=T` records are primitive normalized VCF records;
- external vcfbub/vcfwave behavior is available without compatibility tags,
  fixtures, and provenance;
- normalized downstream records remain inside the same proof boundary as raw
  `VariantCall` emission.

The synthesis task should therefore report the modernization as ready for raw
VCF implementation and fixture design, but not as ready to advertise
consumer-normalized repetitive VCF until the downstream fixture and profile
tasks land.
