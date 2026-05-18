# Lean4 End-to-End Validation Corpus

Date: 2026-05-15
Task: `lean4-e2e-validation-corpus`

This corpus expands the Lean4 conformance harness from one smoke fixture into a
small, reproducible boundary suite for current povu GFA-to-VCF behavior.  The
suite is intentionally fixture-sized: each input is a minimal GFA whose current
`povu gfa2vcf` output is compared against a Lean semantic VCF expectation, or
whose unsupported boundary is expected to fail before any Lean VCF comparison.

The documented conformance command remains:

```bash
cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- --repo-root .
```

For local iteration after `lake build` and `bin/povu` already exist:

```bash
cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- \
  --repo-root . --skip-build --skip-lean-build --povu-bin bin/povu
```

## Corpus

| Fixture id | GFA | Expected boundary | Protected behavior |
| --- | --- | --- | --- |
| `minimal-substitution` | `tests/lean4_conformance/fixtures/minimal_substitution.gfa` | Lean semantic VCF | Ordinary top-level flubble where two haploid paths diverge for one base and rejoin; protects `SUB`, `AT`, `ES`, `LV=0`, and genotype columns. |
| `insertion-flubble` | `tests/lean4_conformance/fixtures/insertion_flubble.gfa` | Lean semantic VCF | Top-level flubble where the alternate path inserts segment `2` after anchor segment `0`; protects anchored `INS` spelling (`A` to `AG`) and insertion traversal formatting. |
| `deletion-flubble` | `tests/lean4_conformance/fixtures/deletion_flubble.gfa` | Lean semantic VCF | Top-level flubble where the alternate path skips reference segment `2`; protects anchored `DEL` spelling (`AG` to `A`) and deletion traversal formatting. |
| `nested-deletion` | `tests/lean4_conformance/fixtures/nested_deletion.gfa` | Lean semantic VCF | Small outer/inner flubble shape derived from the paper's nested-flubble obligation; current povu emits the inner deletion with `LV=1` and a missing genotype for the outer-only sample. |
| `nested-substitution-missing-outer` | `tests/lean4_conformance/fixtures/nested_substitution_missing_outer.gfa` | Lean semantic VCF | Substitution analogue of `nested-deletion`; protects `LV=1` child substitution output and missing genotype semantics for the outer sibling sample. |
| `repeat-anchor-deletion` | `tests/lean4_conformance/fixtures/repeat_anchor_deletion.gfa` | Lean semantic VCF | Homopolymer deletion emitted at the raw graph anchor; protects the decision to delegate repeat left/right normalization downstream. |
| `repeat-anchor-insertion` | `tests/lean4_conformance/fixtures/repeat_anchor_insertion.gfa` | Lean semantic VCF | Homopolymer insertion emitted at the raw graph anchor; complements repeat deletion and protects stable `AT` provenance through ambiguous sequence. |
| `complex-substitution-span` | `tests/lean4_conformance/fixtures/complex_substitution_span.gfa` | Lean semantic VCF | Two-base `CG` to `TA` allele emitted as one graph-faithful `SUB`; protects against accidental primitive decomposition in raw output. |
| `hairpin-inversion-subr` | `tests/lean4_conformance/fixtures/hairpin_inversion_subr.gfa` | Lean semantic VCF | Current SNE/SUBR regression shape through the full `gfa2vcf` path; protects reverse traversal, `SUBR`, and the absence of `ES`/`LV` for hairpin records. |
| `linear-no-variant` | `tests/lean4_conformance/fixtures/linear_no_variant.gfa` | Lean semantic VCF with no records | Two identical paths over a linear graph; protects header/sample handling when no VCF records are emitted. |
| `two-ordered-substitutions` | `tests/lean4_conformance/fixtures/two_ordered_substitutions.gfa` | Lean semantic VCF with strict record order | Two independent top-level flubbles on one contig; protects deterministic record ordering in addition to semantic row equality. |
| `unsupported-overlap` | `tests/lean4_conformance/fixtures/unsupported_overlap.gfa` | Expected controlled povu failure | Non-zero GFA link overlap is outside the accepted corpus subset; protects the unsupported-input boundary by requiring controlled non-success instead of assigning Lean VCF semantics. |
| `malformed-path-missing-overlaps` | `tests/lean4_conformance/fixtures/malformed_path_missing_overlaps.gfa` | Expected controlled povu failure | A `P` line missing the required overlaps column is malformed GFA; protects the parser/input boundary by requiring controlled non-success instead of assigning Lean VCF semantics. |

The same short notes are mirrored near the fixtures in
`tests/lean4_conformance/fixtures/README.md` so future fixture additions have a
nearby convention to follow.  Exact semantic rows and modernization-specific
tool-boundary notes for the expanded practical corpus are documented in
`docs/vcf-modernization/practical_vcf_conformance.md`.

## Harness Registration

The fixture list lives in `tests/lean4_conformance/src/main.rs`.  Successful
fixtures use `ExpectedOutcome::Vcf` and receive a Lean expected output from
`tests/lean4_conformance/lean_reference.lean`.

Most VCF fixtures continue to normalize record order before semantic
comparison, matching the original conformance-harness behavior.  The
`two-ordered-substitutions` fixture sets `check_record_order: true`, so the raw
record sequence emitted by povu must match the Lean call order.  INFO field
ordering and metadata header ordering remain normalized away.

The unsupported and malformed-boundary fixtures use
`ExpectedOutcome::PovuFailure`.  They are not sent to the Lean reference emitter
because the trusted Lean VCF boundary starts after accepted GFA semantics; these
fixtures only check that current povu does not silently produce a VCF for inputs
outside the supported subset.  The harness requires these failures to have a
controlled process exit code so a crash is still reported as a conformance
failure.

## Discovered Follow-Up

During probing, some malformed GFA inputs that should be clean parser-boundary
failures were found to exit with status 139.  A WG follow-up was filed:
`fix-cleanly-reject` ("Fix: cleanly reject malformed GFA in gfa2vcf").

Reproducer captured in that task:

```bash
tmp=$(mktemp --suffix=.gfa)
printf "H\tVN:Z:1.0\nZ\tnot\tvalid\nS\t0\tA\nP\tHG1#1#chr1\t0+\t*\n" > "$tmp"
bin/povu gfa2vcf -i "$tmp" -t 1 -P HG1
```

This corpus includes one malformed and one unsupported non-crashing boundary
fixture now, and leaves the crashing malformed-GFA rejection cases to the
follow-up task.
