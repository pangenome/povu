# Checked Translator Bridge Artifact

Date: 2026-06-12
Task: `bridge-checked-translator`

The checked translator artifact is an executable bridge record produced by the
Rust Lean conformance harness:

```bash
cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- \
  --repo-root . --fixture minimal-substitution \
  --checked-translator build/lean4-conformance/minimal-substitution-witness.json
```

The output JSON schema is `povu.lean4.checked-translator.v1`.  It is intended to
be consumed by downstream bridge checks before any runtime-proof claim is made.
It records the exact fixture, current povu command, Lean theorem boundary,
trusted assumptions, exported current-povu structure, Lean reference structure,
and normalized semantic VCF witness.

## Result Statuses

`checked_semantic_witness` means the fixture was accepted by current povu, the
VCF row semantics matched the Lean reference, the `--structure-export` JSON was
well formed, and the exported C++ decomposition/variant state exactly matched
the Lean structure reference.

`unsupported_diagnostic` means the fixture belongs to the accepted negative
corpus for this bridge.  The artifact names the expected unsupported reason and
records the process exit code plus stdout/stderr, so callers can distinguish a
controlled unsupported boundary from a crash or unexpected success.

## Trusted Assumptions

Every artifact contains a `trusted_assumptions` array.  These assumptions are
part of the machine-readable artifact because the current export is still a
conformance bridge, not a formal runtime proof.  The assumptions currently named
are:

- `accepted_gfa_bytes_refine_lean_document`;
- `structure_export_preserves_cpp_state`;
- `traversal_frame_and_cycle_classes_correct`;
- `hierarchy_laminarity_and_parenthood_correct`;
- `variant_call_extraction_refines_lean_reference`;
- `cost_counters_match_lean_stage_contracts`.

Downstream bridge tasks should retire these by adding lower-level exports,
independent checkers, or Lean proofs.  Until then, a checked translator artifact
is precise fixture evidence that current povu and Lean agree at the exported
semantic boundary, with the remaining trusted base explicitly visible.

## Lean Boundary

The artifact targets:

```lean
PovuLean.Pipeline.semanticGfaToVcf_correct
```

The bridge does not synthesize a Lean proof term.  Instead, it packages the
implementation-side state that currently maps to the theorem inputs:
accepted-GFA records, boundary candidates, PVST hierarchy nodes, semantic
variant calls, and normalized VCF records.  Equality against the Lean reference
structure is the machine-checkable witness that the fixture reaches the same
semantic objects as the proof-owned model.
