# Lean4 End-to-End Fixture Notes

These fixtures are registered by `tests/lean4_conformance/src/main.rs` and
documented in `docs/lean4-proof/e2e_validation.md`.

| Fixture | Protects |
| --- | --- |
| `minimal_substitution.gfa` | Ordinary flubble substitution from two diverging haploid paths. |
| `insertion_flubble.gfa` | Anchored `INS` allele construction and traversal strings. |
| `deletion_flubble.gfa` | Anchored `DEL` allele construction and traversal strings. |
| `nested_deletion.gfa` | Nested flubble-tree VCF fields, including `LV=1` and a missing genotype column. |
| `nested_substitution_missing_outer.gfa` | Nested `SUB` output with `LV=1` and missing genotype semantics for the outer sibling path. |
| `repeat_anchor_deletion.gfa` | Repetitive `AAA` deletion emitted at the raw graph anchor; downstream normalization is deliberately not applied. |
| `repeat_anchor_insertion.gfa` | Repetitive `AAA` insertion emitted at the raw graph anchor; downstream normalization is deliberately not applied. |
| `complex_substitution_span.gfa` | Two-base `SUB` allele (`CG` -> `TA`) kept as one decomposition-sensitive graph record. |
| `hairpin_inversion_subr.gfa` | End-to-end SUBR hairpin inversion output derived from the SNE/SUBR regression shape. |
| `linear_no_variant.gfa` | Header-only VCF output for a graph with no variant records. |
| `two_ordered_substitutions.gfa` | Deterministic ordering of multiple VCF records on one contig. |
| `unsupported_overlap.gfa` | Unsupported non-zero-overlap GFA boundary is rejected instead of compared to Lean semantics. |
| `malformed_path_missing_overlaps.gfa` | Malformed path record missing the GFA overlaps column is rejected instead of compared to Lean semantics. |
