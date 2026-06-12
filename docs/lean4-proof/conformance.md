# Lean4 Conformance Harness

Date: 2026-05-15
Task: `lean4-conformance-harness`

This document describes the first deterministic Rust/Lean conformance harness
for the povu GFA-to-VCF path.  It is intentionally small: it proves the harness
shape, normalization, diagnostics, and trusted boundary without becoming the
broad validation corpus owned by `lean4-e2e-validation-corpus`.

## Command

From a clean checkout with Rust, CMake, a C++17 compiler, and Lean/Lake
available:

```bash
cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- --repo-root .
```

The command runs `lake build`, configures CMake in
`build/lean4-conformance`, builds the `povu` CLI target, runs each fixture
through `povu gfa2vcf`, asks Lean to emit the expected semantic VCF records, and
compares normalized VCF semantics.

Useful development variants:

```bash
cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- \
  --repo-root . --fixture minimal-substitution

cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- \
  --repo-root . --skip-build --skip-lean-build --povu-bin bin/povu

cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- \
  --repo-root . --fixture minimal-substitution \
  --checked-translator build/lean4-conformance/minimal-substitution-witness.json
```

The second form is only for local iteration after `lake build` and `bin/povu`
already exist.

The `--checked-translator` form runs the same C++ CLI and Lean reference paths,
then writes a machine-checkable JSON artifact instead of only printing pass/fail
text.  The artifact schema is `povu.lean4.checked-translator.v1`.  Each
supported fixture result contains:

- the Lean theorem boundary being targeted,
  `PovuLean.Pipeline.semanticGfaToVcf_correct`;
- every trusted assumption still needed before the artifact can be read as a
  formal runtime proof;
- the accepted-GFA structure exported by current povu;
- exported C++ decomposition state visible through `--structure-export`,
  including boundary candidates, PVST nodes, and variant calls;
- the matching Lean structure reference and normalized VCF witness;
- explicit boolean checks showing that the povu VCF and structure export matched
  the Lean semantic references.

Expected rejection fixtures are represented with status
`unsupported_diagnostic`, the precise expected reason, exit code, stdout, and
stderr.  This lets downstream tooling distinguish "checked semantic witness"
from "unsupported by the accepted subset" without scraping harness prose.

## Fixture Coverage

The initial fixture is
`tests/lean4_conformance/fixtures/minimal_substitution.gfa`.

It contains two haploid paths over a four-segment graph:

- `HG1#1#chr1`: `0+,1+,3+`
- `HG2#1#chr1`: `0+,2+,3+`

The harness runs current povu as:

```bash
bin/povu gfa2vcf -i tests/lean4_conformance/fixtures/minimal_substitution.gfa \
  -t 1 -P HG1
```

Current povu emits one `SUB` record with the `HG1` path as reference and the
`HG2` path as alternate.  The Lean side defines the matching
`VCF.VariantCall`, emits it through `VCF.emitRecords`, and serializes only the
semantic row fields needed by the harness.

## Trusted Boundary

The proof-owned semantic boundary is exposed by:

```lean
PovuLean.Pipeline.semanticGfaToVcf_correct
```

That theorem proves semantic VCF emission correctness after the caller supplies:
an accepted semantic `GFA.Document`, traversal and scan witnesses, flubble-tree
support assumptions, and a `VCF.ReferenceCallSet` for semantic
`VCF.VariantCall`s.

The conformance harness does not change those proof obligations.  Its Lean
reference emitter, `tests/lean4_conformance/lean_reference.lean`, imports the
pipeline and VCF emitter modules, constructs the fixture's semantic
`VCF.VariantCall`, checks it with `native_decide`, and calls `VCF.emitRecords`.

## What The Proof Checks

Lean checks that the semantic record emitted from the fixture call is well
formed under the supported VCF subset:

- nonempty `CHROM`, `ID`, `REF`, and alternate allele data;
- positive one-based `POS`;
- fixed `QUAL = 60`, `FILTER = PASS`, and `FORMAT = GT`;
- `SUB` is supported for the flubble source;
- `AC`, semantic `AF`, `AN`, `NS`, `AT`, `VARTYPE`, `TANGLED`, `ES`, and `LV`
  are derived from the semantic call fields;
- genotype allele indexes are valid for the alternate count.

The integrated theorem additionally covers the semantic path from accepted GFA
records and certified algorithm witnesses to semantically correct VCF records.

## What The Test Checks

The Rust harness checks the current implementation boundary that Lean does not
yet prove:

- the byte-level GFA fixture is accepted by current povu;
- the actual `gfa2vcf` command exercises decomposition plus variant calling;
- the serialized VCF row normalizes to the same semantic record as the Lean
  reference;
- sample columns and genotype calls match the Lean fixture expectation.
- when requested, the checked translator artifact can serialize those successful
  comparisons into stable JSON for downstream bridge checks.

The harness intentionally normalizes away byte details that are outside the
trusted theorem and are not meaningful for this fixture:

- `##fileDate` and other metadata header ordering;
- duplicate or reordered metadata lines;
- INFO field ordering within a record;
- record ordering before semantic comparison, except for fixtures that
  explicitly opt into strict record-order checking.

It does not normalize allele spelling, positions, record IDs, INFO values,
sample names, genotype values, `QUAL`, `FILTER`, or `FORMAT`.

The checked translator artifact is not a new proof obligation and does not hide
the remaining trusted base.  Its `trusted_assumptions` array names the current
external assumptions: byte-parser refinement to Lean `GFA.Document`, faithful
serialization of C++ state, unexported traversal-frame and cycle-class
correctness, hierarchy laminarity/parenthood correctness, variant-call
extraction refinement, and missing runtime cost witnesses.  Downstream bridge
tasks can retire those assumptions one at a time by replacing them with direct
checks or proofs.

## Diagnostics

On mismatch, the harness prints:

- fixture id and description;
- the exact `povu gfa2vcf` command;
- the input GFA contents;
- normalized expected semantic VCF from Lean;
- normalized actual semantic VCF from povu;
- raw povu stdout and stderr.

This makes differences actionable without requiring a developer to manually
decode VCF header noise before seeing the semantic disagreement.

## Adding Fixtures

Keep this harness small.  Add only deterministic fixtures needed to exercise a
specific conformance boundary.  Larger corpus expansion belongs to
`lean4-e2e-validation-corpus`.

To add one fixture:

1. Add a GFA file under `tests/lean4_conformance/fixtures/`.
2. Add a semantic call or call list to
   `tests/lean4_conformance/lean_reference.lean` for accepted-GFA VCF fixtures.
3. Add a fixture entry in `tests/lean4_conformance/src/main.rs` with the GFA
   path, reference prefix arguments passed to `povu gfa2vcf`, and expected
   outcome. Unsupported or malformed boundary fixtures should use an expected
   povu failure instead of assigning Lean VCF semantics; crashes still fail the
   harness.
4. Run:

```bash
cargo test --manifest-path tests/lean4_conformance/Cargo.toml
cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- --repo-root .
lake build
```

If the new fixture reveals a real implementation mismatch, keep the diagnostic
output in the downstream issue or WG task so the failing semantic fields are
visible.

The expanded end-to-end corpus and per-fixture notes are maintained in
`docs/lean4-proof/e2e_validation.md`.
