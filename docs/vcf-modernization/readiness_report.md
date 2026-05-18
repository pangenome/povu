# Practical VCF Modernization Readiness Report

Task: `vcf-modern-synthesis-check`

Date: 2026-05-18

## Decision

povu is not ready to make the broad release claim that the whole practical VCF
modernization batch is fully passing its merged Rust/Lean conformance gate.
The join-point command currently fails because four newer practical VCF
fixtures have final VCF expectations but no Lean structure-reference exports.
That unresolved mismatch is tracked as WG task `fix-extend-lean`.

The safe claim today is narrower:

> The Rust semantic VCF writer now implements the raw graph VCF contract for
> already constructed `VariantCall` data, and its unit tests cover the core
> fixture records and formatting/error-policy decisions from the modernization
> spec. The current CLI still agrees with Lean on the structure-covered
> conformance fixtures, and downstream repetitive profile behavior is executable
> in the conformance harness with explicit provenance. A release claim over the
> full practical corpus remains blocked until every positive practical fixture
> has a passing structure comparison or a documented structure-coverage policy.

This claim intentionally does not say that the public Rust one-shot
`povu::gfa_to_vcf` path performs native GFA-to-VCF conversion, that Rust has a
verified graph-to-`VariantCall` extractor, or that default raw output is already
consumer-normalized primitive VCF.

## Evidence Used

| Area | Predecessor artifact or validation | What it supports |
| --- | --- | --- |
| Rust output contract | `docs/vcf-modernization/rust_vcf_output_spec.md` | Raw graph VCF boundary, Rust semantic data model, deterministic header/record behavior, intentional C++ divergences, and downstream-only normalization. |
| Rust implementation | `povu-rs/src/vcf.rs`, `povu-rs/tests/vcf_emitter_tests.rs`; `vcf-modern-rust-emitter` log: `cargo test --manifest-path povu-rs/Cargo.toml` and full conformance passed at that branch point | The semantic Rust VCF writer and focused Rust tests landed. |
| Lean guardrail | `docs/vcf-modernization/lean_rust_guardrail.md` | Lean proves the semantic VCF boundary, not byte parsing, headers, current implementation internals, or Rust binding VCF generation. |
| Structure guardrail | `vcf-modern-structure-export-conformance` artifacts in `tests/lean4_conformance/**`, `src/mto/to_structure_export.cpp`, and related CLI wiring | Current CLI structure exports can be compared against Lean structure exports for the fixtures that have structure references. |
| Practical corpus | `docs/vcf-modernization/practical_vcf_conformance.md`, `tests/lean4_conformance/src/main.rs`, `tests/lean4_conformance/lean_reference.lean`, fixture README | The corpus covers ordinary, nested, repetitive, decomposition-sensitive, hairpin, no-variant, ordering, and rejection behavior at the VCF row level. |
| Repetitive strategy | `docs/vcf-modernization/repetitive_untangling.md` | Popping, left-normalization, vcfwave-like decomposition, star alleles, inserted-path coordinates, and lossy traversal clustering stay downstream of core raw graph output. |
| Downstream profiles | `docs/vcf-modernization/downstream_repetitive_fixtures.md`, `tests/lean4_conformance/src/downstream_profiles.rs` | The conformance harness can execute raw, left-normalized, top-level-only, popped, decomposed/pass-through, and expected-reject profile cases with provenance. |
| C++ divergence and history | `docs/vcf-modernization/rust_vcf_output_inventory.md`, `fix-sne-subr`, `fix-cleanly-reject` | Current C++ behavior remains user-visible context; SUBR final-chunk and malformed-GFA controlled rejection fixes matter for compatibility expectations. |

## Rust VCF Output Changes

The Rust modernization changed the Rust API surface from "VCF-shaped calls
exist but production output is not available through the binding" to an
explicit semantic writer:

- `povu-rs/src/vcf.rs` defines `VcfDocument`, `VariantCall`,
  `VariantSource`, `AlleleConstruction`, `AlternateAllele`, `GenotypeColumn`,
  semantic `Info`, formatted `Record`, `AlleleFrequency`, and `OrderKey`.
- `VcfDocument::to_vcf_string()` and `VcfDocument::write_path()` serialize the
  same validated document.
- Headers are deterministic: `##fileformat=VCFv4.2`, `##source=povu-rs`,
  optional contigs, INFO metadata matching the spec, exactly one `GT` FORMAT
  metadata line, and no dynamic `fileDate` in deterministic output.
- Records use `QUAL=60`, `FILTER=PASS`, `FORMAT=GT`, ordered INFO fields, graph
  traversal provenance in `AT`, flubble `ES`/`LV`, and no `ES`/`LV` for
  hairpin `SUBR`.
- Insertions and deletions are VCF-anchored by the semantic construction:
  `INS(anchor, inserted)` serializes as `REF=anchor`, `ALT=anchor+inserted`;
  `DEL(anchor, deleted)` serializes as `REF=anchor+deleted`, `ALT=anchor`.
- Record ordering follows the semantic key `(contig_order, POS,
  source_primary, source_secondary)`, not C++ chunk or map iteration order.
- Allele frequency rendering keeps useful deterministic precision; the Rust
  tests include `1/3 -> 0.3333333333333333`, avoiding the historical C++
  one-decimal behavior.
- Fully missing phased columns collapse to `.`, while partial missing phases
  such as `0|.` are preserved.
- Invalid semantic inputs return structured `Error::InvalidVcf` or
  `Error::Unsupported` rather than producing plausible empty VCF or panicking.

The Rust tests cover the spec fixture records for minimal substitution,
insertion, deletion, nested deletion, hairpin/SUBR, header-only no-variant
output, deterministic ordering, AF precision, missing genotype rendering,
path/string parity, duplicate record ids, and source/type validation.

Important limitation: this is a semantic Rust VCF writer, not a complete native
Rust GFA-to-VCF implementation. `povu::gfa_to_vcf` still calls the FFI path, and
`GraphAnalysis::write_vcf` now reports that it needs a verified
graph-to-`VariantCall` extractor and points callers to `VcfDocument::write_path`
for semantic documents.

## Still Downstream

The following behavior is not part of core Rust raw graph emission:

- repeat left/right normalization and reference-FASTA-dependent
  `bcftools norm`-style coordinate shifting;
- vcfbub-like parent popping, threshold filtering, and child rescue;
- vcfwave-like WFA/BiWFA decomposition of complex alleles into primitive rows;
- star-allele spanning deletion semantics such as a future `vg -R` compatible
  mode;
- inserted-path coordinate systems for nested children that lack a selected
  reference coordinate;
- lossy traversal clustering or `vg -u`/`-L`-style untangling fields;
- external `vcfbub` compatibility based on `PS`, because default raw output
  uses the current povu/Lean `ES` plus `LV` convention.

The downstream repetitive runner does implement fixture-defined profile
behavior inside the conformance crate:

- `raw-graph` passthrough or controlled reject;
- `left-normalized` coordinate shifting with sidecar provenance;
- `top-level-only` filtering;
- `popped` large-parent removal with child rescue;
- `decomposed` primitive splitting or pass-through with raw ALT provenance;
- expected-reject behavior for nested child inside insertion and inserted-path
  raw mode;
- `SUBR` preservation through decomposition profiles.

Those downstream records are useful compatibility/profile outputs. They are not
inside the Lean raw `VariantCall` emission proof boundary, and any
coordinate/allele/genotype rewrite needs provenance.

## Lean And Rust Agreement

### Proven Core Guarantee

Lean proves semantic VCF emission correctness at
`PovuLean.Pipeline.semanticGfaToVcf_correct`, assuming accepted semantic GFA
input, certified traversal/flubble/hairpin/hierarchy witnesses, and well-formed
semantic `VCF.VariantCall`s. It proves record well-formedness, semantic order,
source consistency, and supported INFO/genotype semantics for emitted semantic
records.

Lean does not prove byte-level GFA parsing, current C++ traversal-frame
construction, current Rust/FFI internals, `.pvst` round trips, decimal text
formatting, header metadata, or serialized VCF bytes.

### Tested VCF Equivalence

The practical corpus is a Rust-run harness around current `povu gfa2vcf` plus a
Lean reference emitter. It compares normalized semantic VCF rows, including
`CHROM`, `POS`, `ID`, `REF`, `ALT`, INFO values, `FORMAT`, samples, and
genotypes. It intentionally ignores `##` metadata such as duplicate C++ `GT`
headers and dynamic `fileDate`.

Positive VCF row coverage includes:

- ordinary flubble substitution;
- anchored insertion and deletion;
- nested deletion and nested substitution with outer-sibling missing genotype;
- repeat-anchor deletion and insertion;
- complex substitution span kept as one graph-faithful `SUB`;
- hairpin/SUBR reverse traversal;
- header-only no-variant output;
- deterministic ordering.

Expected-failure coverage includes unsupported non-zero overlap and malformed
path records.

### Tested Structure Agreement

The structure-export guardrail added `--structure-export` comparison for
current-povu accepted GFA, boundary/PVST nodes, parent links, and
pre-serialization variant-call fields. The join command currently confirms
structure agreement for the structure-covered fixtures before the failure:

- `minimal-substitution`
- `insertion-flubble`
- `deletion-flubble`
- `nested-deletion`
- `hairpin-inversion-subr`
- `linear-no-variant`
- `two-ordered-substitutions`

It does not currently confirm structure agreement for four newer positive VCF
fixtures because Lean has no structure reference for them yet:

- `nested-substitution-missing-outer`
- `repeat-anchor-deletion`
- `repeat-anchor-insertion`
- `complex-substitution-span`

This is not an observed semantic mismatch in those fixtures; it is a missing
structure oracle at the merged join point. It still blocks a full readiness
claim because the required conformance command exits nonzero.

## C++ Divergence And User Impact

Current C++ remains the user-visible `povu gfa2vcf` implementation used by the
Lean conformance harness. Rust output should be described as semantically
aligned with Lean on the Rust `VcfDocument` surface, not byte-compatible with
current C++ headers or a replacement for the C++ CLI.

Intentional Rust-vs-C++ differences that matter for users:

- Rust emits one `GT` FORMAT metadata line; C++ has historically duplicated it.
- Rust declares `AN` as integer metadata; C++ historically declared it as a
  string.
- Rust deterministic mode omits dynamic `fileDate`; C++ emits dynamic metadata.
- Rust renders `AF` with useful precision; C++ uses one decimal place.
- Rust sorts by the semantic order key; C++ generation can be influenced by
  chunking or map order.
- Rust semantic APIs return structured errors; C++ CLI compatibility is process
  exit status plus diagnostics.

Historical C++ fixes still matter:

- `fix-sne-subr` makes final-partial-chunk `SUBR` output current expected
  behavior. The `hairpin-inversion-subr` fixture protects this convention.
- `fix-cleanly-reject` makes malformed or unsupported GFA a controlled failure
  boundary. The conformance corpus relies on nonzero exits rather than signals
  or successful empty VCF.

These divergences are documented rather than hidden. Users should expect row
semantics for covered raw graph cases to line up, while headers, error surface,
ordering, and precision may differ intentionally.

## Validation Run For This Report

Commands run from the repository root on 2026-05-18:

| Command | Result |
| --- | --- |
| `cargo test --manifest-path povu-rs/Cargo.toml` | Passed. Includes 14 Rust VCF emitter tests plus existing Rust binding tests and doctests. |
| `cargo test --manifest-path tests/lean4_conformance/Cargo.toml` | Passed. Includes harness unit tests, downstream repetitive asset validation, and downstream profile tests. |
| `cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- --repo-root . --downstream-repetitive` | Passed. All selected downstream repetitive profiles passed. |
| `lake build` | Passed. No Lean-facing artifacts were changed by this synthesis task. |
| `cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- --repo-root .` | Failed. VCF/structure-covered fixtures pass up to `two-ordered-substitutions`, and controlled rejection fixtures pass, but four newer positive practical fixtures fail at Lean structure-reference lookup. |

Failure excerpt:

```text
Lean structure reference failed for fixture nested-substitution-missing-outer
unknown Lean conformance structure fixture: nested-substitution-missing-outer

Lean structure reference failed for fixture repeat-anchor-deletion
unknown Lean conformance structure fixture: repeat-anchor-deletion

Lean structure reference failed for fixture repeat-anchor-insertion
unknown Lean conformance structure fixture: repeat-anchor-insertion

Lean structure reference failed for fixture complex-substitution-span
unknown Lean conformance structure fixture: complex-substitution-span
```

## Release-Blocking Follow-Up Tasks

| Task | Blocks | Required outcome |
| --- | --- | --- |
| `fix-extend-lean` | Any release claim that the merged practical VCF corpus passes the Rust/Lean structure guardrail | Add Lean structure-reference exports or an explicit harness policy for the four newer positive fixtures, then make `cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- --repo-root .` pass. |

No additional release blocker is required for the narrower claim that
`povu-rs` has a tested semantic `VcfDocument` writer for already constructed
`VariantCall`s. The broader practical modernization claim stays blocked until
`fix-extend-lean` is done and the full join command passes.
