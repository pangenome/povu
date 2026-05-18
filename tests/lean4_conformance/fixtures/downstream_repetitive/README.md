# Downstream Repetitive VCF Fixture Corpus

These assets instantiate the fixture classes from
`docs/vcf-modernization/repetitive_untangling.md` without changing the core
Rust emitter or the C++ implementation.  They are golden conformance inputs and
expected outputs for future downstream profiles.

The corpus deliberately separates two contracts:

- `raw-graph` files are the graph-faithful output shape aligned with
  `docs/vcf-modernization/rust_vcf_output_spec.md`.
- non-raw profile files are downstream rewrites, filters, pass-through records,
  or expected rejects.  These are not claims about current core povu behavior.

`manifest.json` is the machine-readable index.  The Rust integration test
`downstream_repetitive_assets.rs` verifies that every fixture has a raw
expectation, a non-raw profile expectation, concrete input assets, valid VCF or
expected-reject files, and provenance for coordinate-shifting or one-to-many
rewrites.

## Fixture Summary

| Fixture | Raw expectation | Downstream expectation |
| --- | --- | --- |
| `tandem-repeat-left-normalization` | Two right-anchored homopolymer indels remain at `POS=4`. | `left-normalized` shifts both records to `POS=1` and records old/new coordinates in JSONL provenance. |
| `popped-parent-child-rescue` | Raw output contains a large `LV=0` parent and a small `LV=1` child. | `top-level-only` keeps only the parent; `popped` drops the parent and rescues the eligible child with explicit provenance. |
| `vcfwave-complex-decomposition` | One tangled multi-ALT `SUB` record preserves the graph span. | `decomposed` splits ALT 1 into two SNP records and passes ALT 2 through because of the max-length policy. |
| `nested-child-inside-insertion` | Default raw VCF is an expected reject because the child has no selected-reference coordinate. | Future inserted-path raw mode is also an expected reject until the coordinate/header contract exists. |
| `subr-inversion-preservation` | One graph-native `SUBR` record omits `ES`/`LV`. | `decomposed` passes it through with `SUBR` origin metadata rather than emitting `<INV>` or primitive rows. |

## Provenance Rules Captured Here

Coordinate-shifting rewrites must preserve at least the raw record id, raw ALT
index, raw `CHROM/POS/REF/ALT`, output record id, output `CHROM/POS/REF/ALT`,
raw `AT` traversal, profile name, and flags such as `LEFT_NORMALIZED`.

One-to-many decomposition must preserve the raw record id, raw ALT index, raw
allele, all output record ids, profile configuration, and whether any other
raw ALT was intentionally passed through.  Genotypes in decomposed rows must not
silently turn an unrelated raw ALT into reference; the fixtures use `.` for the
sample carrying the non-decomposed overlapping ALT.

Pass-through cases are still profile decisions.  They keep `ORIGIN`,
`RAW_ALT_INDEX`, and `PASSTHROUGH` fields so downstream consumers can
distinguish "left unchanged by policy" from raw graph output.
