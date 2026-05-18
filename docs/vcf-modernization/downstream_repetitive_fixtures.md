# Downstream Repetitive VCF Fixture Assets

Task: `add-downstream-repetitive`

Date: 2026-05-18

## Scope

The downstream repetitive corpus lives at
`tests/lean4_conformance/fixtures/downstream_repetitive/`.  It converts the
fixture classes from `docs/vcf-modernization/repetitive_untangling.md` into
concrete inputs, raw expectations, downstream profile expectations, and
provenance sidecars.

These assets are intentionally not registered as current `povu gfa2vcf`
conformance cases.  The current core contract remains `raw-graph` emission as
described in `docs/vcf-modernization/rust_vcf_output_spec.md`; downstream
popping, left-normalization, and decomposition profiles are not implemented in
the core emitter by this task.

## Validation Command

The fixture corpus is checked by a Rust integration test in the existing
conformance crate:

```bash
cargo test --manifest-path tests/lean4_conformance/Cargo.toml downstream_repetitive_fixture_assets_are_complete
```

The test verifies that:

- every fixture has a raw expectation and a separate downstream expectation;
- VCF expectations contain parseable row structure and required core INFO
  fields;
- expected-reject assets name the rejected profile and required behavior;
- coordinate-shifting rewrites have old/new coordinate provenance;
- one-to-many decomposition has raw ALT index and output record mapping;
- the corpus covers tandem repeat normalization, vcfbub-like child rescue,
  vcfwave-like decomposition/pass-through, nested child inside insertion, and
  SUBR preservation.

## Fixture Matrix

| Fixture | Raw graph contract | Downstream contract |
| --- | --- | --- |
| `tandem-repeat-left-normalization` | Two homopolymer indels remain at the graph-local right anchor, `POS=4`. | `left-normalized` shifts both spellings to `POS=1` against `reference.fa`; `provenance.jsonl` records raw/output coordinates, allele strings, `AT`, and profile configuration. |
| `popped-parent-child-rescue` | Raw output contains a large `LV=0` parent and a small `LV=1` child. | `top-level-only` keeps the parent only; `popped` removes the parent and rescues the child because the child satisfies size thresholds. |
| `vcfwave-complex-decomposition` | One tangled multi-ALT `SUB` preserves graph traversal for both alternates. | `decomposed` maps raw ALT 1 to two SNP rows and keeps raw ALT 2 as a `TANGLED=T` pass-through because it exceeds the configured max allele length. |
| `nested-child-inside-insertion` | Default raw VCF is an expected reject because the child exists only inside an inserted alternate path and has no selected-reference coordinate. | Future inserted-path raw mode is also an expected reject until contig naming, header policy, and origin semantics are defined. |
| `subr-inversion-preservation` | Raw output keeps one graph-native `SUBR` record with reverse traversal and no `ES`/`LV`. | `decomposed` passes the record through with `SUBR_ORIGIN=T`; the fixture forbids symbolic `<INV>` replacement or primitive decomposition without explicit future provenance. |

## Provenance Expectations

Coordinate-shifting and one-to-many profile outputs are downstream rewrites,
not additional raw graph facts.  Their sidecar records must preserve:

- raw record id and raw ALT index;
- raw `CHROM`, `POS`, `REF`, `ALT`, and `AT`;
- output record id or ids;
- output `CHROM`, `POS`, `REF`, and `ALT` where applicable;
- profile name and profile configuration;
- flags such as `LEFT_NORMALIZED`, `DECOMPOSED`, `PASSTHROUGH`,
  `POPPED_PARENT`, and `RESCUED_CHILD`.

The decomposition fixture also documents genotype policy: samples carrying a
different overlapping raw ALT are serialized as missing (`.`) on primitive rows
rather than being silently converted to reference.

## Expected Rejects and Pass-Throughs

`nested-child-inside-insertion` remains expected-reject for both default raw VCF
and future inserted-path raw mode.  This is deliberate: the core emitter must
not fabricate coordinates for a child with no selected-reference coordinate,
and no inserted-path coordinate contract exists yet.

`vcfwave-complex-decomposition` and `subr-inversion-preservation` include
pass-through output rows.  A pass-through is still a downstream profile result;
it carries `ORIGIN`, `RAW_ALT_INDEX`, `PROFILE`, and `PASSTHROUGH` so consumers
can distinguish policy-retained records from raw graph output.

## Implementation Status

These files are conformance assets for downstream implementation.  They do not
change povu's Rust emitter, C++ implementation, or current raw VCF proof
boundary.  A follow-up implementation task should make the downstream profile
runner consume this corpus before any synthesis task claims repetitive VCF
normalization is implemented.
