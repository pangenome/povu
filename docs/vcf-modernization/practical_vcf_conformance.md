# Practical VCF Conformance Corpus

Task: `vcf-modern-conformance-corpus`

Date: 2026-05-18

## Scope

This corpus expands the Rust-run practical VCF conformance harness around the
modernization output spec and the Lean/Rust guardrail.  It is intentionally a
semantic row corpus, not a byte-for-byte header corpus:

- the Rust harness in `tests/lean4_conformance/src/main.rs` drives each GFA
  fixture through `povu gfa2vcf`;
- Lean emits the expected semantic VCF rows from explicit
  `VCF.VariantCall`s in `tests/lean4_conformance/lean_reference.lean`;
- the Rust harness parses both VCF streams and compares normalized record
  semantics, including `CHROM`, `POS`, `ID`, `REF`, `ALT`, `INFO`, `FORMAT`,
  samples, and genotypes;
- current C++ header quirks such as duplicate `GT` metadata and dynamic
  `fileDate` are outside this comparison;
- unsupported input fixtures are accepted only as controlled povu failures.

The trusted boundary is the same as
`docs/vcf-modernization/lean_rust_guardrail.md`: final semantic VCF rows from
accepted flubble/hairpin sources.  This corpus does not claim that current Rust
bindings export VCF yet, and it does not widen the proof boundary to raw PVST,
cycle-class, or hairpin-scan internals.

## Conformance Command

Run the complete expanded corpus from the repository root:

```bash
cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- --repo-root .
```

During local iteration after `bin/povu` has already been built, this equivalent
command keeps the Lean build but skips the CMake rebuild:

```bash
cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- --repo-root . --skip-build
```

Run one fixture by id:

```bash
cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- --repo-root . --fixture repeat-anchor-deletion
```

The command currently covers thirteen fixtures: eleven positive semantic VCF
cases and two controlled rejection cases.

## Common Record Defaults

Every positive record below has:

- `QUAL=60`
- `FILTER=PASS`
- `FORMAT=GT`

INFO fields are listed in semantic order:

```text
AC;AF;AN;NS;AT;VARTYPE;TANGLED;ES;LV
```

`ES` and `LV` are absent for `SUBR` hairpin records.  `linear-no-variant`
expects a valid header and sample columns with no data records.

## Positive Fixtures

### `minimal-substitution`

Fixture file: `tests/lean4_conformance/fixtures/minimal_substitution.gfa`

Expected semantic row:

```text
CHROM=HG1#1#chr1 POS=2 ID=>0>3 REF=C ALT=G
INFO=AC=1;AF=0.5;AN=2;NS=2;AT=>1,>2;VARTYPE=SUB;TANGLED=F;ES=>0>3;LV=0
GT: HG1=0 HG2=1
```

Why it matters: this is the smallest ordinary top-level flubble substitution.
It proves the harness still checks the base row shape used by later Rust
emitter tests.

External-tool note: this aligns with raw graph VCF behavior from
`vg deconstruct`: preserve graph traversal provenance and do not invent a
normalization pass.

### `insertion-flubble`

Fixture file: `tests/lean4_conformance/fixtures/insertion_flubble.gfa`

Expected semantic row:

```text
CHROM=HG1#1#chr1 POS=1 ID=>0>1 REF=A ALT=AG
INFO=AC=1;AF=0.5;AN=2;NS=2;AT=>0,>0>2;VARTYPE=INS;TANGLED=F;ES=>0>1;LV=0
GT: HG1=0 HG2=1
```

Why it matters: checks VCF anchoring for insertions, including a nonempty
reference allele and an alternate traversal that retains the inserted graph
segment.

External-tool note: this is raw anchored output.  `bcftools norm` or a
vcfwave-style normalizer may later rewrite an indel against a FASTA, but that
rewrite is outside the proved core.

### `deletion-flubble`

Fixture file: `tests/lean4_conformance/fixtures/deletion_flubble.gfa`

Expected semantic row:

```text
CHROM=HG1#1#chr1 POS=1 ID=>0>1 REF=AG ALT=A
INFO=AC=1;AF=0.5;AN=2;NS=2;AT=>0>2,>0;VARTYPE=DEL;TANGLED=F;ES=>0>1;LV=0
GT: HG1=0 HG2=1
```

Why it matters: checks VCF anchoring for deletions, including a nonempty ALT
and graph traversal provenance for the deleted segment.

External-tool note: downstream normalization may left-align equivalent indels,
but povu's core output preserves the graph-local anchor.

### `nested-deletion`

Fixture file: `tests/lean4_conformance/fixtures/nested_deletion.gfa`

Expected semantic row:

```text
CHROM=HG1#1#chr1 POS=2 ID=>1>4 REF=CT ALT=C
INFO=AC=1;AF=0.5;AN=2;NS=2;AT=>1>3,>1;VARTYPE=DEL;TANGLED=F;ES=>1>4;LV=1
GT: HG1=0 HG2=1 HG3=.
```

Why it matters: protects nested flubble presentation, especially `LV=1` and
the current missing genotype convention for a sample that takes the outer
sibling path and does not traverse the child site.

External-tool note: this is the povu/Lean equivalent of graph-aware nesting
metadata from `vg deconstruct -n`.  A vcfbub-like pass may later pop parents or
rescue children, but the raw corpus keeps the child record explicit.

### `nested-substitution-missing-outer`

Fixture file:
`tests/lean4_conformance/fixtures/nested_substitution_missing_outer.gfa`

Expected semantic row:

```text
CHROM=HG1#1#chr1 POS=3 ID=>1>4 REF=T ALT=G
INFO=AC=1;AF=0.5;AN=2;NS=2;AT=>3,>6;VARTYPE=SUB;TANGLED=F;ES=>1>4;LV=1
GT: HG1=0 HG2=1 HG3=.
```

Why it matters: extends the nested coverage beyond deletion.  The child site is
a substitution, and the outer sibling sample must still serialize as missing
for this record rather than as reference or an invented spanning allele.

External-tool note: a vcfbub-style popping profile would reason over parent and
child levels after raw emission.  This fixture keeps povu responsible only for
the raw child row and its `ES`/`LV` metadata.

### `repeat-anchor-deletion`

Fixture file: `tests/lean4_conformance/fixtures/repeat_anchor_deletion.gfa`

Expected semantic row:

```text
CHROM=HG1#1#chr1 POS=1 ID=>0>2 REF=AA ALT=A
INFO=AC=1;AF=0.5;AN=2;NS=2;AT=>0>1,>0;VARTYPE=DEL;TANGLED=F;ES=>0>2;LV=0
GT: HG1=0 HG2=1
```

Why it matters: this is the explicit repetitive/ambiguous sequence fixture.
The reference path contains an `AAA` homopolymer before a shared suffix, and the
alternate omits one graph segment.  Multiple VCF spellings are equivalent after
left/right alignment, but the conformance expectation is the raw graph
boundary: povu emits this anchored `DEL`, delegates left-normalization to
downstream tools, and does not reject the fixture.

External-tool note: `bcftools norm` and vcfwave-like normalization may shift or
rewrite this representation when a matching FASTA is supplied.  The core povu
row must not silently perform that rewrite because it would change `POS`,
`REF`, `ALT`, and potentially graph provenance.

### `repeat-anchor-insertion`

Fixture file: `tests/lean4_conformance/fixtures/repeat_anchor_insertion.gfa`

Expected semantic row:

```text
CHROM=HG1#1#chr1 POS=1 ID=>0>2 REF=A ALT=AA
INFO=AC=1;AF=0.5;AN=2;NS=2;AT=>0,>0>1;VARTYPE=INS;TANGLED=F;ES=>0>2;LV=0
GT: HG1=0 HG2=1
```

Why it matters: complements the repeat deletion with a repeat insertion at the
same graph-local anchor.  It protects the first Rust contract's rule that
repeat normalization is delegated and raw `AT` traversal stays stable.

External-tool note: a normalized profile may shift the inserted `A` within the
homopolymer.  The raw graph profile keeps the `>0` anchor and `>0>1` alternate
traversal.

### `complex-substitution-span`

Fixture file: `tests/lean4_conformance/fixtures/complex_substitution_span.gfa`

Expected semantic row:

```text
CHROM=HG1#1#chr1 POS=2 ID=>0>3 REF=CG ALT=TA
INFO=AC=1;AF=0.5;AN=2;NS=2;AT=>1>2,>4>5;VARTYPE=SUB;TANGLED=F;ES=>0>3;LV=0
GT: HG1=0 HG2=1
```

Why it matters: stresses decomposition-sensitive alleles.  The graph difference
could be represented as two primitive substitutions after sequence alignment,
but core povu emits one graph-faithful `SUB` record over the whole flubble
span.

External-tool note: this is where a vcfwave-like WFA decomposition pass is
useful downstream.  The corpus intentionally expects no WFA decomposition in
raw povu output.

### `hairpin-inversion-subr`

Fixture file: `tests/lean4_conformance/fixtures/hairpin_inversion_subr.gfa`

Expected semantic row:

```text
CHROM=ref POS=2 ID=>1>5 REF=ACGTA ALT=TACGT
INFO=AC=1;AF=0.5;AN=2;NS=2;AT=>1>2>3>4>5,<5<4<3<2<1;VARTYPE=SUBR;TANGLED=F
GT: ref=0 alt=1
```

Why it matters: protects the graph-native hairpin/SUBR output convention.  The
record must omit `ES` and `LV` because its source is a hairpin, not a flubble.

External-tool note: SV-oriented or vcfwave-like downstream output may annotate
or decompose inversion-like sequence.  Raw povu keeps the `SUBR` record and its
reverse traversal.

### `linear-no-variant`

Fixture file: `tests/lean4_conformance/fixtures/linear_no_variant.gfa`

Expected semantic output:

```text
samples=HG1,HG2
records=[]
```

Why it matters: a graph with VCF-capable reference/sample metadata and no
selected flubbles must produce a valid header-only VCF, not fail and not invent
records.

External-tool note: this is a compatibility baseline for downstream tools:
empty callsets remain valid VCF streams.

### `two-ordered-substitutions`

Fixture file: `tests/lean4_conformance/fixtures/two_ordered_substitutions.gfa`

Expected semantic rows in raw order:

```text
1. CHROM=HG1#1#chr1 POS=2 ID=>0>3 REF=C ALT=G
   INFO=AC=1;AF=0.5;AN=2;NS=2;AT=>1,>2;VARTYPE=SUB;TANGLED=F;ES=>0>3;LV=0
   GT: HG1=0 HG2=1

2. CHROM=HG1#1#chr1 POS=4 ID=>3>6 REF=A ALT=C
   INFO=AC=1;AF=0.5;AN=2;NS=2;AT=>4,>5;VARTYPE=SUB;TANGLED=F;ES=>3>6;LV=0
   GT: HG1=0 HG2=1
```

Why it matters: this fixture opts into strict record-order comparison.  It
guards the deterministic ordering rule that the Rust emitter must follow:
contig order, `POS`, and then source keys.

External-tool note: this is intentionally stronger than historical C++ chunk or
map iteration behavior.  The Rust modernization contract gives semantic order
priority.

## Controlled Rejection Fixtures

### `unsupported-overlap`

Fixture file: `tests/lean4_conformance/fixtures/unsupported_overlap.gfa`

Expected behavior: povu exits nonzero.  The harness treats this as a pass only
because non-zero GFA link overlap is outside the accepted semantic subset used
by the Lean guardrail.

Why it matters: unsupported inputs must reject before semantic VCF comparison;
they must not produce an empty VCF that looks successful.

### `malformed-path-missing-overlaps`

Fixture file:
`tests/lean4_conformance/fixtures/malformed_path_missing_overlaps.gfa`

Expected behavior: povu exits nonzero.  The fixture has a malformed GFA path
line missing the overlaps column.

Why it matters: malformed GFA syntax must remain a boundary failure, not a
semantic VCF case.

## Coverage Map

| Concern | Fixtures |
| --- | --- |
| Ordinary flubble substitution | `minimal-substitution` |
| Anchored insertion/deletion | `insertion-flubble`, `deletion-flubble` |
| Nested/complex flubbles | `nested-deletion`, `nested-substitution-missing-outer` |
| Repetitive ambiguous sequence | `repeat-anchor-deletion`, `repeat-anchor-insertion` |
| Decomposition-sensitive alleles | `complex-substitution-span` |
| Hairpin/SUBR-like structure | `hairpin-inversion-subr` |
| No-variant graphs | `linear-no-variant` |
| Deterministic ordering | `two-ordered-substitutions` |
| Boundary rejection | `unsupported-overlap`, `malformed-path-missing-overlaps` |

## Modernization Guidance

The corpus is deliberately conservative:

- povu emits raw graph records for the repeat fixtures;
- povu delegates repeat left/right normalization and complex-allele
  decomposition to downstream profiles or external tools;
- povu rejects unsupported/malformed GFA before semantic comparison;
- no fixture requires Rust or C++ to export intermediate PVST/flubble/hairpin
  structures, because that remains the separate
  `vcf-modern-structure-export-conformance` dependency;
- no new Rust implementation mismatch was found while expanding this corpus.

Future Rust emitter work should port these exact `VariantCall` expectations
into Rust API tests once the semantic Rust VCF model exists.  If that emitter
disagrees with any row here, the fix belongs in `vcf-modern-rust-emitter` unless
the synthesis task explicitly decides to file a blocking WG task for a changed
contract.
