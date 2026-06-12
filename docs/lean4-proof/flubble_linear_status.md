# Flubble Linear-Time Proof Status

Date: 2026-06-12
Task: `synthesis-flubble-linear`

This report summarizes the flubble linear-time proof batch after
`bridge-povu-implementation`.  It is a user-facing status note, not a new proof.
The Lean theorem boundary is semantic: the current formal results prove the
Lean reference definitions and conditional cost compositions.  Current povu C++
and Rust behavior has conformance evidence, but there is not yet a
machine-checkable implementation refinement or runtime proof for the production
implementation.

## Bottom Line

The batch proves the semantic flubble detector, hierarchy output-size bounds,
indexed-detector equivalence, exact indexed lookup/update counts, and
flubble-stage linearity in Lean.  The algorithms-facing theorem is that the
indexed flubble detector/extractor plus flat hierarchy construction is linear
in the candidate stack under the standard word-RAM class-index convention; once
the candidate stack is bounded linearly by graph size, the stage is linear in
`|V| + edgeCount`.  The full decomposition theorem is conditional on named
stage contracts and intermediate-size bounds.  The current C++ and Rust
implementations are validated against fixture and port-level conformance tests,
but the project must still build a concrete cost model plus a checked bridge
from implementation state to Lean witnesses before claiming that current povu
itself has a formal linear-runtime proof.

## Formally Proved In Lean

These items build in the trusted Lean modules under `PovuLean`.

| Area | Module path | Important names | Status |
| --- | --- | --- | --- |
| Flubble semantic definitions | `PovuLean.Algorithms.Flubble.Spec` (`PovuLean/Algorithms/Flubble/Spec.lean`) | `Boundary`, `Boundary.CanonicalIds`, `Boundary.ordered`, `firstSameClass?`, `closeAfterGap?`, `CanonicalClassBoundary`, `IsPaperBoundary`, `IsFlubbleBoundary`, `NoDuplicateBoundaries` | Defines the formal semantic notion of a canonical flubble boundary over a traversal candidate stack. |
| Supported flubble inputs and class correctness | `PovuLean.Algorithms.Flubble.InputInvariant` (`PovuLean/Algorithms/Flubble/InputInvariant.lean`) | `IsBoundaryCandidate`, `candidateStack`, `SupportedInput`, `CycleClassAssignment`, `CycleClassAssignment.Correct`, `SupportedGFAInput`, `toSupportedInput` | Defines the proof-side assumptions for real black tree-edge candidates and sound/complete cycle-class assignments. |
| Reference detector output bounds | `PovuLean.Algorithms.Flubble.Detect` (`PovuLean/Algorithms/Flubble/Detect.lean`) | `detectFlubbles_length_le_candidateStack_length`, `candidateStack_length_le_treeLinks_length`, `detectFlubbles_length_le_treeLinks_length`, `treeLinks_length_le_graph_links_length`, `detectFlubbles_length_le_graph_links_length`, `detectFlubbles_length_le_graph_edgeCount`, `detectFlubbles_length_le_supportedInput_graph_edgeCount`, `detectFlubbles_noDuplicates` | Proves the detector emits at most one deduplicated boundary per stack candidate and bounds output by tree links, graph links, and `edgeCount` under the stated input invariants. |
| Reference detector correctness | `PovuLean.Algorithms.Flubble.Correctness` (`PovuLean/Algorithms/Flubble/Correctness.lean`) | `detectFlubbles_sound`, `detectFlubbles_complete`, `detectFlubbles_canonical_noDuplicates`, `detectFlubbles_correct`, `detectFlubbles_correct_for_gfa` | Proves soundness, completeness, duplicate freedom, and GFA semantic-input corollary for the Lean reference detector. |
| Indexed detector semantic equivalence | `PovuLean.Algorithms.Flubble.Indexed` (`PovuLean/Algorithms/Flubble/Indexed.lean`) | `lookup_indexedEntriesFrom_eq_firstSameClass?`, `indexedBoundaryAt?_eq_closeAfterGap?`, `detectStackIndexedRaw_eq_detectStackRaw`, `detectStackIndexed_eq_detectStack`, `detectFlubblesIndexed_eq_detectFlubbles`, `detectFlubblesIndexed_iff_isFlubbleBoundary`, `detectFlubblesIndexed_noDuplicates` | Proves the one-pass indexed detector returns the same boundaries as the reference detector. |
| Indexed detector operation counts | `PovuLean.Algorithms.Flubble.Indexed` (`PovuLean/Algorithms/Flubble/Indexed.lean`) | `scanStackIndexedAux_visited`, `scanStackIndexedAux_lookups`, `scanStackIndexedAux_updates`, `detectStackIndexedRaw_one_lookup_and_update_per_candidate` | Proves exactly one lookup and one update per candidate for the indexed scan shape. |
| Hierarchy semantics and correctness | `PovuLean.Algorithms.FlubbleTree.Spec`, `PovuLean.Algorithms.FlubbleTree.Correctness` (`PovuLean/Algorithms/FlubbleTree/Spec.lean`, `PovuLean/Algorithms/FlubbleTree/Correctness.lean`) | `SupportedHierarchyInput`, `UnsupportedHierarchyInput`, `IsForestHierarchy`, `IsCorrectHierarchy`, `buildHierarchyFrom_correct`, `buildHierarchy_correct`, `buildHierarchy_rejects_unsupported` | Proves hierarchy correctness for supported laminar boundary spans and rejects treating unsupported/non-laminar inputs as supported. |
| Hierarchy output-size bounds | `PovuLean.Algorithms.FlubbleTree.Correctness` (`PovuLean/Algorithms/FlubbleTree/Correctness.lean`) | `buildHierarchy_nodes_length`, `buildHierarchy_nodes_length_le_detectFlubbles_length`, `buildHierarchy_parentReferenceSlots_length`, `buildHierarchy_boundaryReferenceSlots_length`, `buildHierarchy_flatMetadataSlots_length`, `buildHierarchy_nodes_length_le_supportedInput_graph_edgeCount`, `buildHierarchy_flatMetadataSlots_le_supportedInput_graph_edgeCount_twice` | Proves hierarchy node count equals detector output count and flat parent/boundary metadata is linear in detector output and graph edge count under supported input. |
| Semantic VCF pipeline correctness | `PovuLean.Pipeline` (`PovuLean/Pipeline.lean`) | `semanticGfaToVcf_correct` | Proves semantic VCF emission correctness when accepted semantic GFA, traversal/class/scan witnesses, supported hierarchy input, and semantic calls are supplied. |

## Conditionally Proved In Lean

These are theorem statements proved by Lean, but their hypotheses are explicit
contracts.  They should not be read as facts about the current executable
implementation until a later bridge supplies those contracts for that code.

| Area | Module path | Important names | Conditions |
| --- | --- | --- | --- |
| Cost vocabulary | `PovuLean.Complexity.Flubble` (`PovuLean/Complexity/Flubble.lean`) | `GraphSizes`, `DecompositionSizes`, `StageCost`, `LinearStageWith`, `LinearStage`, `ClassIndexContract`, `IndexedDetectorOperationContract`, `TraversalFrameSizeContract`, `OutputSizeReuseContract` | Defines proof-side sizes and cost contracts; it does not assert implementation costs. |
| Indexed operation linearity | `PovuLean.Complexity.Flubble` (`PovuLean/Complexity/Flubble.lean`) | `IndexedDetectorOperationContract.operationCost_linear_in_candidates` | Requires a `ClassIndexContract` with constant lookup/update costs and exact lookup/update counts. |
| Output-size reuse | `PovuLean.Complexity.Flubble` (`PovuLean/Complexity/Flubble.lean`) | `OutputSizeReuseContract.flubbles_linear_in_input`, `OutputSizeReuseContract.hierarchy_linear_in_input` | Requires output-size bounds and an output-linear-in-input hypothesis. |
| Boundary and hierarchy output-size composition | `PovuLean.Complexity.Flubble` (`PovuLean/Complexity/Flubble.lean`) | `boundaryOutputSize_linear_in_candidateStack`, `hierarchyConstructionOutputSize_linear_in_candidateStack` | The boundary result uses the indexed detector count bound.  The hierarchy result additionally requires `FlubbleTree.SupportedHierarchyInput`. |
| Indexed flubble-stage linearity | `PovuLean.Complexity.Flubble` (`PovuLean/Complexity/Flubble.lean`) | `ExtractionPipelineCosts`, `ConditionalExtractionCostContract`, `conditional_extraction_linear_in_candidateStack`, `indexed_flubble_stage_linear_in_candidateStack`, `indexed_flubble_stage_linear_in_graphSize` | Requires constant-time class-index operations, scan bookkeeping linear in candidates, boundary emission linear in emitted boundaries, hierarchy construction linear in hierarchy output size, supported hierarchy input, and, for the graph-size corollary, a candidate-stack linear-size bound.  The conclusion is linear in `candidateStack` length and then in graph size once that size bound is supplied. |
| Abstract decomposition-stage arithmetic | `PovuLean.Complexity.Flubble` (`PovuLean/Complexity/Flubble.lean`) | `DecompositionStageCosts.detectorPrefix_linear_with`, `DecompositionStageCosts.hierarchySuffix_linear_with`, `DecompositionStageCosts.total_linear_with` | Composes supplied per-stage linear bounds; does not supply those bounds for an implementation. |
| Conditional full semantic decomposition | `PovuLean.Complexity.Decomposition` (`PovuLean/Complexity/Decomposition.lean`) | `SemanticDecomposeIntermediateSizes`, `SemanticDecomposeIntermediateSizeBounds`, `SemanticDecomposePipelineCosts`, `SemanticDecomposeStageContracts`, `conditional_semantic_decompose_linear` | Requires accepted semantic GFA, supported flubble/hierarchy/hairpin witnesses, semantic call witnesses, intermediate-size bounds, and linear contracts for component decomposition, tip/dummy augmentation, traversal-frame construction, cycle-class assignment, hairpin scan, and flubble extraction/hierarchy.  The conclusion is a proof-side `LinearBound costs.total sizes.inputGraphSize` plus semantic VCF correctness. |

The external stage names required by
`PovuLean.Complexity.Decomposition.conditional_semantic_decompose_linear` are:
`componentDecompositionLinear`, `tipDummyAugmentationLinear`,
`traversalFrameConstructionLinear`, `cycleClassAssignmentLinear`,
`hairpinScanLinear`, and `flubbleExtractionAndHierarchy`.  The first five are
textbook-linear graph/tree traversal or stack/indexing facts that are not yet
mechanized in Lean.  The sixth is discharged by the indexed flubble-stage
theorem above.

## Validated By C++/Rust Conformance Tests

Conformance tests are evidence that current implementation behavior matches
selected Lean semantic expectations.  They are not Lean theorems about C++ or
Rust runtime.

| Evidence | Path or command | Current status |
| --- | --- | --- |
| Documented Lean/C++ conformance harness | `cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- --repo-root .` | Passed on 2026-06-12. The run built Lean with `lake build`, configured/built the C++ `povu` CLI in `build/lean4-conformance`, then passed all 13 fixtures: `minimal-substitution`, `insertion-flubble`, `deletion-flubble`, `nested-deletion`, `nested-substitution-missing-outer`, `repeat-anchor-deletion`, `repeat-anchor-insertion`, `complex-substitution-span`, `hairpin-inversion-subr`, `linear-no-variant`, `two-ordered-substitutions`, `unsupported-overlap`, and `malformed-path-missing-overlaps`. |
| Lean conformance crate tests | `cargo test --manifest-path tests/lean4_conformance/Cargo.toml` | Passed on 2026-06-12: harness unit tests, downstream repetitive asset validation, and downstream profile tests all passed. |
| Rust binding/native semantic tests | `cargo test --manifest-path povu-rs/Cargo.toml` | Passed on 2026-06-12. Notable flubble tests in `povu-rs/tests/native_gfa_vcf_tests.rs` include `flubble_stack_port_matches_lean_close_after_gap_rule`, `flubble_stack_port_preserves_lean_first_later_same_class_rule`, and `flubble_stack_port_canonicalizes_and_deduplicates_like_lean_detector`. |
| Bridge audit artifact | `docs/lean4-proof/bridge_povu_implementation.md` | Records the implementation-stage map for component decomposition, tip/dummy augmentation, traversal-frame construction, bracket/cycle-class assignment, flubble extraction, hierarchy construction, and VCF/structure export. It explicitly treats implementation evidence as conformance only. |

The upstream bridge task also logged successful `cargo test --manifest-path
povu-rs/Cargo.toml`, `cargo test --manifest-path
tests/lean4_conformance/Cargo.toml`, the full Lean/Rust conformance harness,
and C++ `test_povu`/`ctest` validation before commit `248b1a5`.

## Still Open

The following blockers prevent a formal runtime claim for current povu:

- No checked translator relates byte-level GFA parsing and current C++ component
  decomposition to Lean `GFA.Document.Accepted` and `GFA.Document.toGraph`.
- Current C++ spanning-tree, tip/dummy augmentation, and traversal state are not
  exported as a Lean `TraversalFrame` witness with proof that the C++ stack order
  equals `Flubble.candidateStack`.
- Current C++ bracket classes are not proved to satisfy
  `Flubble.CycleClassAssignment.Correct`; soundness and completeness against
  `CycleEquivalentEdges` remain external.
- C++ flubble extraction is not machine-checked against Lean
  `closeAfterGap?`, `firstSameClass?`, `Boundary.ordered`, and
  `detectFlubbles` for arbitrary supported inputs.  The Rust semantic port has
  focused conformance tests, but this is not an implementation proof.
- C++ hierarchy parent edges are not proved to match Lean
  `FlubbleTree.SupportedHierarchyInput` and canonical strict span-containment
  parent semantics for arbitrary supported laminar inputs.
- No concrete runtime cost model connects C++/Rust loops, maps, stacks,
  allocations, FFI boundaries, or serialization to the Lean `StageCost`,
  `ExtractionPipelineCosts`, and `SemanticDecomposePipelineCosts` records.
- Container-operation assumptions for `std::vector`, `std::unordered_map`,
  `std::unordered_set`, `std::map`, `std::list`, Rust collections, and
  serialization buffers are not named in a checked implementation cost artifact.
- Existing conformance fixtures are finite regression evidence; they do not
  quantify over all supported semantic graphs or all accepted byte-level inputs.

The bridge task created follow-up WG tasks for these gaps:
`bridge-export-flubble`, `bridge-check-cycle`, `bridge-check-hierarchy`,
`bridge-cost-instrumentation`, and `bridge-checked-translator`.

## Validation Recorded For This Report

Commands rerun in this worktree on 2026-06-12:

- `lake build`: passed; 39 Lean jobs built successfully.
- `rg -n '\b(sorry|admit|axiom|opaque|unsafe)\b' PovuLean`: passed with no
  matches.  The command exits with status 1 when `rg` finds no matches.
- `cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- --repo-root .`:
  passed; all Lean/C++ conformance fixtures passed.
- `cargo test --manifest-path tests/lean4_conformance/Cargo.toml`: passed;
  conformance harness tests and downstream profile tests passed.
- `cargo test --manifest-path povu-rs/Cargo.toml`: passed; Rust binding,
  native GFA/VCF, and VCF emitter tests passed.

## Claim Boundary

It is accurate to say that Lean proves conditional semantic linearity for the
flubble extraction/hierarchy portion and for the full semantic decomposition
pipeline under named contracts.  It is also accurate to say that current C++ and
Rust behavior is conformance-tested against the Lean semantic boundary on the
maintained fixtures.

It is not yet accurate to say that current povu has a formal linear-runtime
proof.  That claim requires a concrete implementation cost model and a
machine-checkable bridge proving that current implementation executions supply
the Lean semantic witnesses, size bounds, and stage-cost contracts consumed by
`PovuLean.Complexity.Decomposition.conditional_semantic_decompose_linear`.
