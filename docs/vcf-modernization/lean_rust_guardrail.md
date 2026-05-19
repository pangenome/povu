# Lean/Rust Guardrail for VCF Modernization

Task: `vcf-modern-lean-rust-guardrail`

Date: 2026-05-18

## Executive Conclusion

The current guardrail is strong enough to say this:

- Lean proves the semantic reference boundary for accepted GFA semantics,
  certified flubble/hairpin witnesses, supported laminar flubble hierarchy
  inputs, and well-formed semantic `VCF.VariantCall`s.
- The existing Rust conformance harness checks that current `povu gfa2vcf`
  output normalizes to the same semantic VCF records as Lean for the registered
  fixture corpus.

It is not strong enough to say this:

- current Rust or C++ intermediate flubble/PVST/variant structures have been
  proved correct;
- current implementation traversal frames, cycle-class stacks, hairpin scan
  state, or flubble hierarchy nodes have been compared directly to Lean;
- current Rust binding VCF generation agrees with Lean, because the Rust FFI
  VCF path is not implemented.

The practical VCF modernization work can use the existing guardrail for final
semantic VCF row preservation on the covered fixtures. Final synthesis must not
claim implementation-level flubble/PVST structural correctness until a
canonical intermediate export is compared. The dependent task is
`vcf-modern-structure-export-conformance`, now a prerequisite of
`vcf-modern-synthesis-check`.

## Lean-Proved Boundary

The trusted Lean join point is
`PovuLean.Pipeline.semanticGfaToVcf_correct`
(`PovuLean/Pipeline.lean:29-60`). It composes the completed Lean modules and
returns `VCF.EmissionCorrect` for `VCF.emitRecords calls`. Its assumptions are
explicit: accepted semantic GFA input, a traversal frame, a flubble
cycle-class assignment, a hairpin scan assignment, a supported flubble hierarchy
input, and a `VCF.ReferenceCallSet` for the semantic calls
(`PovuLean/Pipeline.lean:29-54`).

### Flubble Boundaries

Lean models a base flubble as a canonical pair of real black tree edges in a
certified candidate stack. A boundary is a pair of `openEdge` and `closeEdge`
(`PovuLean/Algorithms/Flubble/Spec.lean:18-22`). The stack is the ordered filter
of real black tree links from a supplied traversal frame
(`PovuLean/Algorithms/Flubble/InputInvariant.lean:20-41`). The paper-facing
predicate requires both edges to be boundary candidates, distinct, canonical by
edge id, and cycle-equivalent
(`PovuLean/Algorithms/Flubble/Spec.lean:142-149`).

The main theorem `detectFlubbles_correct` proves three facts under
`SupportedInput` and `CycleClassAssignment.Correct`:

- soundness: every emitted boundary satisfies `IsFlubbleBoundary`;
- completeness: every supported canonical flubble boundary is emitted;
- duplicate freedom for emitted boundaries.

See `PovuLean/Algorithms/Flubble/Correctness.lean:168-225`. The GFA-facing
corollary carries the same result for accepted semantic GFA documents
(`PovuLean/Algorithms/Flubble/Correctness.lean:227-242`).

What remains outside this proof is also explicit. The Lean flubble module takes
the traversal frame and class assignment as certified inputs. It does not prove
that current povu's DFS/bracket implementation produces that frame, that its
equivalence classes satisfy `CycleClassAssignment.Correct`, or that its stack
order is the same canonical stack order
(`PovuLean/Algorithms/Flubble/InputInvariant.lean:6-11`,
`PovuLean/Algorithms/Flubble/Correctness.lean:3-10`).

### Flubble Hierarchy

Lean builds hierarchy from detected flubble boundaries using stack spans. A
boundary span is the interval between the boundary endpoints in the candidate
stack (`PovuLean/Algorithms/FlubbleTree/Spec.lean:21-84`). Parenthood means
strict span containment, and supported hierarchy inputs require every boundary
to have a span and every pair to be laminar
(`PovuLean/Algorithms/FlubbleTree/Spec.lean:91-116`). Non-laminar pairs and
missing endpoint spans are named unsupported inputs
(`PovuLean/Algorithms/FlubbleTree/Spec.lean:118-134`).

The theorem `buildHierarchy_correct` proves that, given the flubble detector
proofs and `SupportedHierarchyInput`, the hierarchy is correct:

- hierarchy nodes correspond exactly to detected flubble boundaries;
- boundaries are duplicate-free;
- parent links are valid nesting links or absent when no nesting parent exists;
- detector order and canonical boundary orientation are preserved for
  downstream emitters.

See `PovuLean/Algorithms/FlubbleTree/Spec.lean:152-205` and
`PovuLean/Algorithms/FlubbleTree/Correctness.lean:124-176`.

Outside Lean: current povu must still demonstrate that its concrete PVST forest
matches this hierarchy. The production PVST path serializes only selected
route-visible data (`src/mto/to_pvst.cpp:50-105`) and reconstructs flubble-like
families from file symbols in `from_pvst` (`src/mto/from_pvst.cpp:207-260`).
The pipeline inventory notes that richer metadata such as computed bounds is
not fully represented after PVST round-trip
(`docs/lean4-proof/povu_pipeline_inventory.md:439-458`).

### Hairpin Boundaries

Hairpin support is part of the integrated VCF theorem. Lean specifies a
structural hairpin boundary with two real candidate edges, canonical orientation,
and a graph-structural reverse-stem witness
(`PovuLean/Algorithms/Hairpin/Spec.lean:76-135`). The theorem
`detectHairpins_correct` proves soundness, completeness, and duplicate freedom
under `SupportedInput` plus `HairpinScanAssignment.Correct`
(`PovuLean/Algorithms/Hairpin/Correctness.lean:109-164`).

Outside Lean: current povu still must supply or compare the reverse scan
certificate. The hairpin module says the concrete reverse DFS/bracket state is a
later conformance obligation
(`PovuLean/Algorithms/Hairpin/Spec.lean:11-15`,
`PovuLean/Algorithms/Hairpin/Spec.lean:181-197`).

### Semantic VCF Emission

Lean's VCF layer starts after graph variant sources and allele spellings have
already been supplied. `VariantCall` carries the source, coordinates, REF/ALT
construction, traversal strings, counts, genotype columns, and source-specific
INFO fields (`PovuLean/VCF/Spec.lean:343-414`). The supported subset includes
`DEL`, `INS`, `SUB`, and `SUBR`, fixed `QUAL`, `FILTER`, and `FORMAT`, semantic
`AC`/`AF`/`AN`/`NS`/`AT`/`VARTYPE`/`TANGLED`, and optional `ES`/`LV`
(`PovuLean/VCF/Spec.lean:16-31`).

`recordOfCall_wellFormed`, `emitRecords_ordered`, and `emitRecords_correct`
prove that semantic records emitted from well-formed calls are well formed,
ordered, and derived from verified flubble-tree or hairpin sources
(`PovuLean/VCF/Emit.lean:31-110`,
`PovuLean/VCF/Correctness.lean:45-73`,
`PovuLean/VCF/Correctness.lean:126-156`). The GFA-facing theorem
`emitVcfRecords_correct_for_gfa` composes this with flubble-tree and hairpin
correctness (`PovuLean/VCF/Correctness.lean:157-202`).

Outside Lean: byte-level GFA parsing, extraction of DNA strings from current
graph slices, decimal formatting for `AF`, header rendering, and serialized VCF
text are not proved. The VCF spec calls those downstream conformance obligations
(`PovuLean/VCF/Spec.lean:4-14`), and the pure emitter says tabular text
conversion, `fileDate`, and decimal rendering are outside the verified layer
(`PovuLean/VCF/Emit.lean:3-10`).

## Current Rust-vs-Lean Conformance

The current conformance binary is Rust, in
`tests/lean4_conformance/src/main.rs`. It is a harness around the current
`povu` CLI, not a proof of Rust implementation internals.

The harness:

- runs `lake build` unless skipped (`tests/lean4_conformance/src/main.rs:65-70`);
- builds the `povu` CLI with CMake unless skipped
  (`tests/lean4_conformance/src/main.rs:77-79`,
  `tests/lean4_conformance/src/main.rs:309-329`);
- runs each VCF fixture through `povu gfa2vcf -t 1 -P ...`
  (`tests/lean4_conformance/src/main.rs:468-482`);
- asks Lean to emit the expected fixture VCF with
  `lake env lean --run tests/lean4_conformance/lean_reference.lean <fixture>`
  (`tests/lean4_conformance/src/main.rs:456-466`);
- parses both outputs into normalized semantic VCF records
  (`tests/lean4_conformance/src/main.rs:396-410`,
  `tests/lean4_conformance/src/main.rs:530-599`);
- sorts records before comparison unless a fixture opts into strict order
  (`tests/lean4_conformance/src/main.rs:601-617`);
- treats unsupported or malformed boundary fixtures as expected controlled
  povu failures rather than Lean semantic VCF cases
  (`tests/lean4_conformance/src/main.rs:425-453`).

The positive fixture set currently covers:

- minimal top-level substitution;
- top-level insertion;
- top-level deletion;
- nested deletion with `LV=1` and a missing genotype column;
- hairpin/SUBR inversion output;
- no-variant/header-only output;
- two ordered substitutions with strict record-order checking.

These registrations are in `tests/lean4_conformance/src/main.rs:207-291`; the
fixture notes mirror the same scope in
`tests/lean4_conformance/fixtures/README.md:6-16`. The Lean reference constructs
the expected `VCF.VariantCall`s directly and checks their well-formedness with
`native_decide` examples before emitting `VCF.emitRecords`
(`tests/lean4_conformance/lean_reference.lean:116-154`,
`tests/lean4_conformance/lean_reference.lean:162-245`,
`tests/lean4_conformance/lean_reference.lean:253-325`).

The harness comparison is meaningful for final semantic VCF rows. It does not
normalize away allele spelling, positions, IDs, INFO values, sample names,
genotypes, `QUAL`, `FILTER`, or `FORMAT`; the original conformance document
states that boundary explicitly (`docs/lean4-proof/conformance.md:94-115`).

## What Current Tests Do Not Exercise

The current conformance suite does not compare any of these implementation
structures directly against Lean:

- parsed GFA as a `GFA.Document` plus `Document.Accepted`;
- `Core.TraversalFrame` and ordered real-black candidate stack;
- `CycleClassAssignment.Correct` evidence or raw class ids;
- canonical flubble boundaries before PVST construction;
- `FlubbleTree.Hierarchy` nodes, parent links, spans, or detector order;
- `HairpinScanAssignment` state and hairpin boundary candidates;
- pre-serialization `VCF.VariantCall` objects;
- PVST serialization/round-trip fidelity against Lean hierarchy semantics;
- subflubble families (`tiny`, `parallel`, `concealed`, `midi`, `smothered`);
- Rust binding VCF behavior.

The last point matters for wording. The Rust binding exposes flubble counts and
a minimal PVST vertex count, but detailed flubble access returns `nullptr`
(`povu-rs/povu-ffi/povu_ffi.cpp:397-410`) and VCF generation reports "not yet
implemented" (`povu-rs/povu-ffi/povu_ffi.cpp:440-456`,
`povu-rs/povu-ffi/povu_ffi.cpp:504-521`,
`povu-rs/src/analysis.rs:26-48`). The binding implementation notes list
detailed flubble iteration and VCF generation as TODOs
(`povu-rs/IMPLEMENTATION.md:121-130`,
`povu-rs/IMPLEMENTATION.md:202-214`).

## Required Structure Export

To compare Rust/current-povu intermediate structures against Lean instead of
only final VCF semantics, add a canonical deterministic export with at least
these fields:

- fixture/input identity and normalized accepted-input status;
- ordered real-black candidate stack entries with edge ids, orientation, color,
  provenance, and implementation cycle-class id;
- detected base flubble boundaries with open/close ids, canonical orientation,
  stack positions, and detector order;
- hierarchy/PVST nodes with node id, family, boundary endpoints, route
  direction, parent id, child order, level, and stable ordering key;
- hairpin boundaries and enough scan state to map them to the Lean hairpin
  boundary predicate;
- pre-serialization variant-call fields: source kind and source id, chrom,
  contig order, one-based position, record id, REF, reference traversal,
  alternate construction/traversal/counts, variant type, tangled flag,
  `ES`/`LV`, allele counts/frequencies as semantic ratios, and genotype columns.

The comparison target should be one of:

- a Lean reference structure export produced before `VCF.emitRecords`; or
- an explicitly documented canonical semantic export shared by Lean and the
  implementation, with VCF text formatting excluded from the comparison.

This is deeper than a small fixture addition because current intermediate data
are not exposed through the Rust binding or conformance harness. The dependent
WG task is `vcf-modern-structure-export-conformance`.

## Downstream Guardrail Rules

For `vcf-modern-output-spec`:

- It may rely on Lean for the semantic record contract and on the current harness
  for fixture-level final VCF semantic agreement.
- It must not describe current Rust binding VCF output as implemented.
- It must not claim that current implementation PVST/flubble hierarchy
  structures are Lean-conformant unless the intermediate export task lands.

For `vcf-modern-conformance-corpus`:

- It should add final-VCF fixtures in the current harness when the behavior under
  test is allele spelling, coordinates, genotype columns, INFO values, order, or
  controlled rejection.
- It should not try to infer structural correctness solely from final VCF rows;
  structure-sensitive cases should wait for or depend on
  `vcf-modern-structure-export-conformance`.

For `vcf-modern-synthesis-check`:

- Safe final claim today: "Lean proves the semantic reference boundary, and the
  current CLI agrees with Lean final VCF semantics on the maintained conformance
  fixtures."
- Unsafe final claim today: "Rust/current povu intermediate flubble or variant
  structures are formally correct."
- Synthesis is now explicitly blocked on
  `vcf-modern-structure-export-conformance` for any structure-level claim.

## Validation Record

No Lean or Rust conformance fixtures were added in this task; the missing
intermediate comparison is not small fixture scope. Validation performed for
this artifact:

- reviewed Lean flubble, flubble-tree, hairpin, VCF, and pipeline theorem
  boundaries;
- reviewed the Rust conformance harness and Lean reference emitter;
- confirmed the trusted Lean source tree contains no `sorry` or `admit` matches
  with `rg -n '\b(sorry|admit)\b' PovuLean tests/lean4_conformance/lean_reference.lean`;
- `lake build` passed;
- `cargo test --manifest-path tests/lean4_conformance/Cargo.toml` passed;
- `cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- --repo-root .`
  passed all registered fixtures;
- filed `vcf-modern-structure-export-conformance` as the dependent intermediate
  structure comparison task and added it before `vcf-modern-synthesis-check`.
