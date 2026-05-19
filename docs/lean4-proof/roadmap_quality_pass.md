# Lean4 Proof Roadmap Quality Pass

Date: 2026-05-15
Task: `lean4-proof-quality-pass`

This pass reviewed the Lean4 proof roadmap tasks for the objective of porting
povu's GFA-to-VCF pipeline to Lean4, proving the algorithmic boundary, and
using the proof-backed reference as the validation basis for future language
targets.

The WG graph metadata is the source of truth for task dependencies and
descriptions. This document records the quality decisions made so downstream
workers can see the intended order and ownership model without reverse
engineering the graph.

## Dependency Sequence

The roadmap now follows this staged order:

1. `lean4-proof-quality-pass`
2. `lean4-paper-map` and `lean4-povu-pipeline-inventory`
3. `lean4-formalization-architecture`
4. `lean4-scaffold`
5. `lean4-core-graph-model`
6. `lean4-gfa-spec`
7. `lean4-flubble-correctness`
8. `lean4-hairpin-correctness` and `lean4-flubble-tree-correctness`
9. `lean4-vcf-semantics`
10. `lean4-proof-synthesis-check`
11. `lean4-conformance-harness`
12. `lean4-e2e-validation-corpus`
13. `lean4-language-targeting-plan`

The most important correction was replacing scaffold-only readiness for proof
tasks with stable interface gates. GFA waits for the core graph model; flubble
waits for GFA; hairpin and flubble-tree work wait for the flubble interfaces;
VCF waits for hairpin and flubble-tree; conformance waits for the synthesis
check rather than individual proof modules.

## File Ownership

Task descriptions were tightened with explicit file ownership so parallel work
does not require editing the same files:

- `lean4-paper-map` owns `docs/lean4-proof/paper_theorem_map.md`.
- `lean4-povu-pipeline-inventory` owns
  `docs/lean4-proof/povu_pipeline_inventory.md`.
- `lean4-formalization-architecture` owns
  `docs/lean4-proof/architecture.md`.
- `lean4-scaffold` owns Lean build metadata and top-level import scaffolding
  such as `lean-toolchain`, `lakefile.lean`, `lake-manifest.json`,
  `PovuLean.lean`, and minimal source-tree setup.
- `lean4-core-graph-model` owns `PovuLean/Core/**`.
- `lean4-gfa-spec` owns `PovuLean/GFA/**`.
- `lean4-flubble-correctness` owns `PovuLean/Algorithms/Flubble/**`.
- `lean4-hairpin-correctness` owns `PovuLean/Algorithms/Hairpin/**`.
- `lean4-flubble-tree-correctness` owns
  `PovuLean/Algorithms/FlubbleTree/**`.
- `lean4-vcf-semantics` owns `PovuLean/VCF/**`.
- `lean4-proof-synthesis-check` owns the cross-module integration point and
  `docs/lean4-proof/proof_obligation_check.md`.
- `lean4-conformance-harness` owns the conformance harness and
  `docs/lean4-proof/conformance.md`.
- `lean4-e2e-validation-corpus` owns end-to-end fixtures and
  `docs/lean4-proof/e2e_validation.md`.
- `lean4-language-targeting-plan` owns
  `docs/lean4-proof/language_targets.md`.

Tasks that need to change files owned by an upstream module should add a
blocking WG follow-up or document a named obligation instead of silently
refactoring another task's scope.

## Validation Check

Every downstream task in the batch now has:

- a clear scope statement;
- explicit predecessor expectations in the WG graph and description;
- a concrete `## Validation` section;
- deliverable or file ownership guidance;
- a rule for handling cross-module changes without creating parallel file
  conflicts.

No task was split during this pass. Broad tasks were tightened by scope,
ownership, and dependency gates rather than decomposed, because the roadmap
already has the necessary implementation, synthesis, conformance, corpus, and
language-targeting stages.
