# Rust VCF Output Inventory

Date: 2026-05-18

Task: `vcf-modern-rust-output-inventory`

This inventory targets Rust modernization.  The Rust surfaces today are the
`povu-rs` binding crate and the Rust Lean conformance harness.  Production
GFA-to-VCF output is still emitted by the C++ CLI, so C++ behavior is included
only as compatibility context and as user expectation inherited by any future
Rust output.

No C++/Rust/Lean implementation files were changed for this inventory.

## Rust Target Summary

There is no implemented Rust VCF emitter today.

The public Rust API advertises two VCF routes:

- `povu::gfa_to_vcf(...)` in `povu-rs/src/lib.rs`.
- `GraphAnalysis::write_vcf(...)` in `povu-rs/src/analysis.rs`.

Both routes currently fail before producing VCF text.  The one-shot API calls
the C FFI function `povu_gfa_to_vcf`, whose implementation returns
`gfa_to_vcf not yet implemented in FFI layer`.  The analysis API returns a
Rust-side `Error::Povu` with `VCF generation not yet implemented` without
calling the FFI VCF handle path.

The existing Rust conformance harness under `tests/lean4_conformance` does not
exercise `povu-rs` output.  It is a Rust test runner that executes the C++ CLI
`povu gfa2vcf`, parses the resulting VCF, normalizes selected textual details,
and compares it to Lean semantic fixture expectations.

## Rust Files and Roles Inspected

| File | Role in current Rust GFA/flubble/output surface |
| --- | --- |
| `povu-rs/src/lib.rs` | Crate root and public one-shot `gfa_to_vcf` function.  Converts paths to C strings, calls `ffi::povu_gfa_to_vcf`, and maps FFI failure to `Error`.  The module docs and feature list still advertise VCF generation even though the FFI path is stubbed. |
| `povu-rs/src/graph.rs` | Main safe wrapper around `PovuGraph`.  Provides `load`, in-memory `new`, `add_vertex`, `add_edge`, stubbed `add_path`, reference selection methods, graph topology queries, and `analyze`.  `analyze` returns `GraphAnalysis` but does not give Rust any allele/VCF record data. |
| `povu-rs/src/analysis.rs` | Wraps detected flubbles as `GraphAnalysis` and exposes `flubble_count`, `pvst_tree`, and `write_vcf`.  `write_vcf` is not implemented and always returns an error. |
| `povu-rs/src/ffi.rs` | Bindgen include for C ABI definitions generated from `povu-rs/povu-ffi/povu_ffi.h`; also defines `PovuError::default`. |
| `povu-rs/src/error.rs` | Converts FFI `PovuError` into Rust `Error::Povu` and frees C-allocated error strings. |
| `povu-rs/src/path.rs` | Defines Rust `Orientation`, `Step`, and `Path`; includes simple PanSN-style helpers (`sample_name`, `haplotype`, `contig`) used by topology callers, not by a VCF emitter. |
| `povu-rs/src/vertex.rs`, `povu-rs/src/edge.rs` | Topology helper types.  Relevant to future allele traversal display but not connected to VCF records today. |
| `povu-rs/build.rs` | Builds C++ libraries through CMake with `POVU_BUILD_FFI=ON`, links `povu_ffi`, `mto`, and `povulib`, then generates Rust bindings with bindgen. |
| `povu-rs/povu-ffi/povu_ffi.h` | C ABI for Rust.  Declares opaque `PovuGraph`, `PovuFlubbles`, `PovuPvstTree`, and `PovuVcfOutput` plus graph, flubble, PVST, and VCF functions.  The VCF function declarations are present but documented as not fully implemented. |
| `povu-rs/povu-ffi/povu_ffi.cpp` | Native backing for the Rust FFI.  Loads GFA via `mto::from_gfa::to_bd`, exposes topology and reference methods, has incomplete flubble/PVST access, and stubs all VCF output paths. |
| `tests/lean4_conformance/src/main.rs` | Rust conformance runner.  Runs `povu gfa2vcf` through `std::process::Command`, parses C++ CLI VCF output into `NormalizedVcf`, normalizes INFO order and usually record order, and compares to Lean fixture output.  This is a Rust harness, not Rust VCF output. |
| `tests/lean4_conformance/lean_reference.lean` | Lean executable used by the Rust harness to emit fixture-specific semantic VCF text after `VCF.emitRecords`. |

## Rust Entry Points and Data Structures

### `povu::gfa_to_vcf`

`povu-rs/src/lib.rs` exposes:

```rust
pub fn gfa_to_vcf(
    gfa_path: impl AsRef<std::path::Path>,
    vcf_path: impl AsRef<std::path::Path>,
    ref_file: Option<impl AsRef<std::path::Path>>,
) -> Result<()>
```

Current behavior:

- Converts `gfa_path`, `vcf_path`, and optional `ref_file` to C strings.
- Calls `ffi::povu_gfa_to_vcf`.
- Returns `Ok(())` only if the FFI call returns `true`.
- Current FFI implementation creates a C++ `core::config`, sets task
  `gfa2vcf`, input GFA, labels, references, and optional reference file, then
  returns `false` with `gfa_to_vcf not yet implemented in FFI layer`.

Modernization implication:

- This is the most obvious Rust one-shot API to make functional.
- A future implementation must decide whether it delegates to the C++ CLI/core,
  ports the output path to Rust, or emits from a Lean-backed semantic adapter.

### `PovuGraph`

`PovuGraph` wraps a raw `*mut ffi::PovuGraph`.  It can be built from GFA or
manually:

- `PovuGraph::load(path)` calls `povu_graph_from_gfa`.
- `PovuGraph::new(vertex_capacity, edge_capacity, path_capacity)` calls
  `povu_graph_new` and panics if the pointer is null.
- `add_vertex` and `add_edge` call FFI graph builders.
- `add_path` calls `povu_graph_add_path`, but the native FFI currently returns
  `false` because path construction is still a TODO.
- `set_references_from_file` and `set_references_from_prefixes` update native
  config only.
- `analyze` calls `povu_graph_find_flubbles` and wraps the result in
  `GraphAnalysis`.

Modernization risks:

- Manual graph construction cannot currently add reference paths, so it cannot
  create VCF-capable in-memory inputs through the Rust API.
- `GraphAnalysis` has no lifetime tied to `PovuGraph` even though the native
  `PovuFlubbles` stores a non-owning graph pointer.  Current `write_vcf` returns
  before using that pointer, but a future Rust VCF implementation must address
  this ownership/lifetime hazard.
- The native `povu_graph_find_flubbles` path builds a spanning tree with the
  right size but has a TODO saying it still needs proper initialization before
  calling `povu::flubbles::find_flubbles`.  That makes current Rust flubble
  counts/PVST data unsuitable as trusted VCF input.

### `GraphAnalysis` and `PvstTree`

`GraphAnalysis` wraps `*mut ffi::PovuFlubbles` and exposes:

- `flubble_count()` via `povu_flubbles_count`.
- `pvst_tree()` via `povu_flubbles_get_pvst_tree`.
- `write_vcf(path)`, which currently returns `Error::Povu { code: 1,
  message: "VCF generation not yet implemented" }`.

`PvstTree` only exposes `vertex_count()`.

Modernization implication:

- The Rust API has no `VcfRecord`, allele, traversal, genotype, INFO, or output
  writer type.  The natural data model still has to be designed.
- If `GraphAnalysis::write_vcf` remains the high-level API, it needs access to
  the graph, reference selection, flubble/PVST hierarchy, and full variant
  records, not only a flubble count and PVST vertex count.

### FFI VCF Handles

The C ABI declares:

- `PovuVcfOutput`
- `povu_flubbles_call_variants`
- `povu_vcf_write_to_file`
- `povu_vcf_to_string`
- `povu_vcf_free`
- `povu_gfa_to_vcf`

Current behavior:

- `povu_flubbles_call_variants` returns null and sets
  `VCF generation not yet implemented in FFI layer`.
- `povu_vcf_write_to_file` and `povu_vcf_to_string` can serialize a
  `PovuVcfOutput`, but no implemented path constructs one.
- `povu_gfa_to_vcf` returns `false` and sets
  `gfa_to_vcf not yet implemented in FFI layer`.

Modernization implication:

- The FFI shape anticipates string and file VCF output, but the actual semantic
  record production and C++ writer integration are missing.
- A Rust-native output spec should not treat `PovuVcfOutput` as established
  behavior beyond its intended ownership and serialization shape.

## Current Output Conventions

### Rust `povu-rs`

Current Rust output convention is failure:

- `povu::gfa_to_vcf` returns `Err(Error::Povu { code: 1, message:
  "gfa_to_vcf not yet implemented in FFI layer" })`.
- `GraphAnalysis::write_vcf` returns `Err(Error::Povu { code: 1, message:
  "VCF generation not yet implemented" })`.
- No VCF file or VCF string is emitted through `povu-rs`.

This is surprising because `povu-rs/src/lib.rs`, `povu-rs/README.md`, and
`povu-rs/IMPLEMENTATION.md` present VCF generation as a public feature or
convenience API, while the implementation and limitations sections identify it
as missing.

### Rust Lean Conformance Harness

The Rust harness output convention is semantic comparison, not production
emission:

- It runs `povu gfa2vcf -i <fixture> -t 1 -P <prefix>` using `Command`.
- It parses the VCF text into `NormalizedVcf { samples, records }`.
- It ignores all `##` metadata header lines.
- It stores INFO fields in a `BTreeMap`, normalizing INFO field order.
- It preserves record order only for fixtures whose `check_record_order` is
  true.  All current success fixtures except `two-ordered-substitutions`
  normalize record order.
- Expected-failure fixtures require controlled nonzero process exit with a
  process status code.  A signal crash is still a harness failure.

The harness currently covers C++ CLI behavior and Lean comparison semantics.  It
does not cover the `povu-rs` API, FFI VCF stubs, or a Rust-native emitter.

### C++ CLI Output Context

The C++ CLI is the behavior users actually see today when running
`povu gfa2vcf` or `povu call`.  Any Rust modernization that claims compatibility
should account for these current conventions:

- `gfa2vcf` is a two-step flow: decompose GFA into temporary PVST files, then run
  `call` on that temporary forest.
- With no output directory, `gfa2vcf`/`call` write one combined VCF stream to
  stdout.  With an output directory, C++ splits output into `<sample>.vcf` files
  by reference prefix.
- Header writer emits `##fileformat=VCFv4.2`, dynamic `##fileDate=<today>`,
  `##source=povu`, INFO metadata, contig lines, and the `#CHROM` column header.
- The C++ header currently emits the `GT` FORMAT metadata line twice.
- `##INFO=<ID=AN,...>` is declared with `Type=String`, although VCF convention
  and the emitted value are numeric.
- Records use fixed `QUAL=60`, `FILTER=PASS`, and `FORMAT=GT`.
- INFO is emitted in this order: `AC`, `AF`, `AN`, `NS`, `AT`, `VARTYPE`,
  `TANGLED`, then `ES` and `LV` for non-`SUBR` records.
- `SUBR` records omit `ES` and `LV`; this is also modeled by Lean.
- `AF` is rendered with one decimal place by C++ (`fmt::format("{:.1f}", v)`).
- All-missing genotype phases in a sample column are collapsed to a single `.`;
  otherwise phases are joined with `|`.
- Record grouping is by reference id in `std::map` order, with per-reference
  vectors retaining generation order from chunked RoV processing and SNE output.

## Concrete Surprises and Output Risks

1. Rust advertises VCF output but cannot produce it.
   - `povu-rs/src/lib.rs` exposes and documents `gfa_to_vcf`.
   - `povu-rs/src/analysis.rs` documents `analysis.write_vcf`.
   - `povu-rs/povu-ffi/povu_ffi.cpp` returns not-implemented errors for both
     `povu_flubbles_call_variants` and `povu_gfa_to_vcf`.
   - No fixture currently asserts this Rust failure mode, so the public API
     mismatch is easy to miss.

2. The Rust flubble-to-output precondition is not reliable yet.
   - `povu-rs/povu-ffi/povu_ffi.cpp` constructs a `povu::spanning_tree::Tree`
     by size in `povu_graph_find_flubbles`, leaves a TODO to initialize it
     properly, then calls `povu::flubbles::find_flubbles`.
   - `povu_flubbles_get` always returns null, so callers cannot inspect
     flubble boundaries, walks, or types.
   - `PovuGraph::add_path` is stubbed, so manually built Rust graphs cannot
     supply the path/reference data VCF output needs.

3. Header conformance is weaker than row conformance.
   - The C++ writer emits duplicate `GT` FORMAT metadata lines.
   - The Rust conformance harness ignores all `##` metadata lines before
     comparison, so duplicate FORMAT metadata, dynamic `fileDate`, `source`,
     contig order, and INFO metadata typing are not checked by the Lean harness.
   - This matters for a Rust emitter because it must decide whether to preserve
     the duplicate line for byte compatibility, drop it for standards hygiene,
     or explicitly version the change.

4. `AF` formatting is not semantically exact.
   - The C++ `VcfRec::get_af` renders one decimal place.
   - Lean models allele frequencies semantically as counts over total alleles
     and leaves decimal rendering outside the verified layer.
   - Current Lean fixtures mainly exercise `1/2 -> 0.5`, so they do not expose
     precision loss such as `1/3 -> 0.3`.

5. Record order is only partly protected.
   - The Rust harness normalizes record order for most success fixtures.
   - Only `two-ordered-substitutions` sets `check_record_order: true`.
   - The C++ writer emits by chunk/ref-map generation order, so a Rust emitter
     must choose whether to match current order exactly or define a new
     deterministic semantic order.

6. Missing genotype output collapses all-missing phases.
   - C++ writes `.` for a sample column when every phase is missing; otherwise
     it joins phases with `|`.
   - Lean represents missing phases semantically as `GenotypeAllele.missing`.
   - The `nested-deletion` fixture protects a one-phase missing sample, but
     varying ploidy and fully missing multi-phase columns need explicit Rust
     specification.

7. `SUBR` metadata policy is intentionally different from flubble records.
   - C++ omits `ES`/`LV` for `SUBR`; Lean requires this same absence for
     hairpin records.
   - The `hairpin-inversion-subr` fixture protects this convention.
   - A Rust emitter must not blindly attach flubble hierarchy metadata to all
     records.

## Lean Semantic VCF Assumptions

The Lean VCF layer is a semantic record layer, not a byte-level VCF writer.

Key assumptions from `PovuLean/VCF/Spec.lean`,
`PovuLean/VCF/Emit.lean`, `PovuLean/VCF/Correctness.lean`, and
`PovuLean/Pipeline.lean`:

- The trusted boundary starts from accepted semantic GFA records and verified
  upstream witnesses.  Byte parsing, C++/Rust graph construction, and serialized
  VCF formatting are outside the theorem.
- `VariantCall` is the input to the semantic VCF emitter.  It already carries
  source variant, chrom, one-based positive position, nonempty id, REF,
  reference traversal, alternate allele constructions, variant type, tangled
  flag, optional `ES`/`LV`, reference allele count, and genotype columns.
- Supported variant types are `DEL`, `INS`, `SUB`, and `SUBR`.
- Flubble-tree sources support `DEL`, `INS`, and `SUB`; hairpin sources support
  only `SUBR`.
- Insertions and deletions are anchored so serialized REF and ALT alleles are
  nonempty.  Substitutions and reverse substitutions use supplied REF/ALT
  sequences directly.
- Semantic records fix `QUAL=60`, `FILTER=PASS`, and `FORMAT=GT`.
- INFO consistency is semantic: `AC`, `AF`, `AN`, `NS`, `AT`, `VARTYPE`,
  `TANGLED`, and `ES`/`LV` match the `VariantCall`.
- `ES` and `LV` are required for non-`SUBR` flubble records and absent for
  `SUBR` hairpin records.
- Records are required to be ordered by Lean's `OrderKey` over contig order,
  position, source primary key, and source secondary key.
- `PovuLean.Pipeline.semanticGfaToVcf_correct` exposes the final theorem:
  accepted GFA plus traversal/class/scan witnesses, supported hierarchy input,
  and a well-formed ordered `ReferenceCallSet` produce semantically correct VCF
  records through `VCF.emitRecords`.

Current Rust gaps relative to Lean:

- `povu-rs` has no `VariantCall` or `Record` equivalent.
- `povu-rs` has no Rust allele construction type for anchored `INS`/`DEL`,
  `SUB`, or `SUBR`.
- `povu-rs` has no Rust genotype column model that preserves ploidy and missing
  phase semantics for VCF.
- `povu-rs` cannot currently extract flubble boundaries, walks, traversal
  strings, nested levels, or hairpin/SUBR sources from `GraphAnalysis`.
- The Rust conformance harness maps C++ CLI text to a normalized comparison
  object, but it does not prove or test a Rust API mapping to Lean `VariantCall`
  assumptions.

## C++ Divergence and Prior Fixes

These notes are included only to set user expectations for future Rust output.
They are not modernization targets for this task.

### Production C++ Emits VCF; Rust Does Not

The most important divergence is feature completeness.  C++ `povu gfa2vcf`
already produces user-visible VCF.  Rust `povu-rs` exposes VCF-shaped APIs but
returns not-implemented errors.

Downstream Rust modernization should decide whether the first supported Rust
behavior is:

- compatibility wrapper around the existing C++ writer,
- Rust-native emitter matching current C++ row conventions,
- Rust-native semantic emitter aligned with Lean, with explicit differences from
  current C++ formatting,
- or a staged combination.

### Prior C++ Fix: SNE/SUBR Final-Chunk Gating

WG task `fix-sne-subr` fixed a prior C++ issue where SNE/SUBR generation was
effectively gated on `chunk_num == CHUNK_SIZE`.  Current C++ code now runs
`ise::sne(...)` exactly once when processing the final RoV chunk, because SNE
consumes the cumulative pin cushion and per-chunk calls would duplicate earlier
pins.

Evidence in current tree:

- `src/ita/genomics/genomics.cpp` computes `final_chunk = end == N` and calls
  `ise::sne(g, pc, to_call_ref_ids)` only when `final_chunk`.
- `tests/integration_tests/genomics_tests.cc` includes
  `GenomicsTest.RunsSneForFinalPartialChunkAndEmitsSingleSubr`, which expects
  one `SUBR` record for a one-chunk inversion fixture.
- `tests/lean4_conformance/fixtures/hairpin_inversion_subr.gfa` and the
  registered `hairpin-inversion-subr` conformance fixture protect the end-to-end
  `SUBR` row convention.

Rust implication:

- Future Rust output should treat final-chunk `SUBR` availability as current
  expected C++ behavior, not as an open C++ bug to work around.

### Prior C++ Fix: Malformed-GFA Crash Handling

WG task `fix-cleanly-reject` fixed prior malformed-GFA crashes in the C++
`gfa2vcf` boundary.  Current C++ code preflights GFA records before handing them
to liteseq, checks liteseq status/pointers, and the `gfa2vcf` command catches
exceptions and exits with `EXIT_FAILURE` after cleanup.

Evidence in current tree:

- `src/mto/from_gfa.cpp` rejects unsupported record types and malformed/empty
  segment sequences before or immediately after liteseq parsing.
- `app/subcommand/gfa2vcf.cpp` catches `std::exception` and unknown exceptions,
  removes the temp directory, reports an error, and exits failure.
- `tests/integration_tests/gfa2vcf_tests.cc` checks unsupported-record and
  missing-sequence inputs exit without a signal and include `Invalid GFA`.
- The Lean conformance corpus has expected-failure fixtures for unsupported
  non-zero overlap and malformed path records.

Rust implication:

- A Rust emitter/runner should report unsupported or malformed GFA as controlled
  failures, not as panics, null-pointer crashes, or successful empty VCF output.
- The Rust API should probably surface these as structured `Result::Err`
  variants rather than process exits.

## Tests and Fixtures Covering Output Today

### `povu-rs`

Current Rust binding tests do not cover successful VCF output.

Observed coverage:

- `povu-rs/tests/builder_tests.rs` covers in-memory graph creation, vertex/edge
  addition, finalization, topology queries, orientation combinations, empty and
  long sequences, self loops, and a comparison against `tests/data/simple.gfa`
  if that file exists.
- `povu-rs/tests/graph_tests.rs` covers GFA loading, topology queries,
  `analyze`, reference prefix setting, reverse complement helpers, edge helper
  behavior, and nonexistent-file errors.
- Many `graph_tests.rs` and one builder comparison test look for
  `tests/data/simple.gfa`, but the current repository only has
  `tests/data/LPA.gfa`, so those tests return early rather than exercising GFA
  loading.
- There are no tests asserting the current `povu::gfa_to_vcf` or
  `GraphAnalysis::write_vcf` not-implemented behavior.
- There are no Rust fixtures for VCF header fields, REF/ALT spelling, INFO
  fields, genotype columns, ordering, or malformed-GFA behavior through
  `povu-rs`.

### Rust Lean Conformance Harness

`tests/lean4_conformance/src/main.rs` is the strongest current Rust test code
related to VCF semantics, but it validates the C++ CLI rather than `povu-rs`.

Registered success fixtures:

- `minimal-substitution`: ordinary top-level `SUB`, `AT`, `ES`, `LV=0`, and
  genotype columns.
- `insertion-flubble`: anchored `INS` (`A` to `AG`) and traversal formatting.
- `deletion-flubble`: anchored `DEL` (`AG` to `A`) and traversal formatting.
- `nested-deletion`: inner deletion with `LV=1` and a missing genotype for the
  outer-only sample.
- `hairpin-inversion-subr`: reverse traversal `SUBR` and absence of `ES`/`LV`.
- `linear-no-variant`: header/sample handling with no records.
- `two-ordered-substitutions`: strict record ordering for two top-level calls.

Registered expected-failure fixtures:

- `unsupported-overlap`: non-zero GFA link overlap is outside the accepted
  subset.
- `malformed-path-missing-overlaps`: malformed path record missing the overlaps
  column is rejected before Lean VCF comparison.

Harness comparison behavior:

- Parses `#CHROM` and records.
- Ignores `##` metadata.
- Normalizes INFO field order.
- Normalizes record order except when `check_record_order` is true.
- Requires expected failures to be controlled nonzero exits, not signals.

### C++ Tests Relevant to User Expectations

Although not Rust target tests, these current C++ tests explain behavior Rust
users will expect if Rust claims compatibility:

- `tests/integration_tests/genomics_tests.cc`:
  `RunsSneForFinalPartialChunkAndEmitsSingleSubr`.
- `tests/integration_tests/gfa2vcf_tests.cc`:
  `RejectsUnsupportedRecordWithoutSignal` and
  `RejectsSegmentMissingSequenceWithoutSignal`.

## Suggested Inputs for `vcf-modern-output-spec`

The downstream Rust output spec should make these decisions explicit:

- Whether Rust's first supported VCF path is file output, string output, or both.
- Whether Rust output is required to be byte-compatible with current C++ or only
  semantically compatible with Lean.
- Header policy: duplicate `GT` FORMAT line, `fileDate`, contig lines, INFO
  metadata types, and metadata order.
- `AF` formatting policy, especially for non-terminating decimals.
- Record order policy: C++ generation order, Lean semantic `OrderKey`, or a new
  documented deterministic order.
- Exact Rust data structures for `VariantCall`, `Record`, `Info`, genotype
  columns, allele constructions, traversal strings, and output writers.
- How Rust obtains verified or conformance-checked flubble/hairpin/PVST data:
  through fixed FFI extraction, a native Rust port, or a Lean-backed adapter.
- Error policy for malformed/unsupported GFA: structured Rust errors should
  replace C++ process exits at the API boundary.
- Test plan requiring `povu-rs` VCF tests, not only the Rust conformance harness
  that runs the C++ CLI.
