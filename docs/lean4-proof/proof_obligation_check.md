# Lean4 Proof Obligation Check

Date: 2026-05-15
Task: `lean4-proof-synthesis-check`

This report is the synthesis checkpoint for the Lean proof modules completed by
the core graph, GFA, flubble, hairpin, flubble-tree, and VCF tasks.  The checked
trusted Lean join point is `PovuLean.Pipeline`, imported from the root
`PovuLean.lean` aggregator.

## Status Summary

`lake build` succeeds for the integrated Lean project after adding
`PovuLean/Pipeline.lean`.

The trusted proof path contains no `sorry` or `admit` tokens under `PovuLean`.
No draft Lean module with placeholders is imported by the trusted path.

The Lean reference now covers the semantic GFA-to-semantic VCF boundary stated
by `PovuLean.Pipeline.semanticGfaToVcf_correct`: accepted semantic GFA records,
certified traversal/class/scan witnesses, supported flubble hierarchy inputs,
and well-formed pipeline-reported variant calls produce semantically correct
VCF records.  It does not yet prove that current povu byte parsing,
DFS/bracket implementation, allele extraction, or serialized VCF text conform
to those semantic witnesses; those are blocking follow-ups for
`lean4-conformance-harness`, with broader fixture expansion left to
`lean4-e2e-validation-corpus`.

## Integrated Modules

| Module family | Responsibility in the integrated path | Status |
| --- | --- | --- |
| `PovuLean.Core.Basic` | Core segment, oriented segment, link, graph, color, provenance, graph well-formedness, and bidirected closure predicates. | Integrated. |
| `PovuLean.Core.Walk` | Walks, paths, cycles, reachability, component views, edge cuts, and traversal-frame hooks consumed by algorithm modules. | Integrated; `Core.TODO.TraversalFrameSpanningForest` is a named future strengthening, not a placeholder proof. |
| `PovuLean.GFA.Basic` | Normalized semantic GFA records, accepted subset checks, graph construction, path support, and accepted-document construction theorem. | Integrated. |
| `PovuLean.Algorithms.Flubble.*` | Supported input contract, cycle-class assignment contract, canonical flubble detector, soundness, completeness, and no-duplicate theorem. | Integrated. |
| `PovuLean.Algorithms.Hairpin.*` | Hairpin detector built on the flubble candidate-stack interface, certified scan assignment, soundness, completeness, and no-duplicate theorem. | Integrated. |
| `PovuLean.Algorithms.FlubbleTree.*` | Stack-span hierarchy specification, supported laminar input contract, hierarchy builder, and hierarchy correctness theorem. | Integrated. |
| `PovuLean.VCF.*` | Semantic VCF subset, variant calls, semantic emitter, record well-formedness, source derivation, ordering, and GFA-facing emission theorem. | Integrated. |
| `PovuLean.Pipeline` | Synthesis-owned join point exposing `semanticGfaToVcf_correct` over the completed modules. | Added by this task. |
| `*/Examples.lean` modules | Build-checked examples for local coverage of each module family. | Built by Lake; not a Rust/Lean conformance corpus and not extended here. |

## Theorem Obligations

| Area | Proved in the trusted path | Delegated or externally supplied | Still open |
| --- | --- | --- | --- |
| Core graph and walks | `Graph.BidirectedWellFormed`, `IsWalk.reverse`, `Walk.reverse`, and reachability composition are available to downstream modules. | Concrete DFS/spanning-forest construction is represented by `TraversalFrame` data supplied to algorithms. | A full traversal-frame spanning-forest theorem remains a future strengthening. |
| GFA semantics | `Document.accepted_toGraph_bidirectedWellFormed` proves accepted semantic documents construct bidirected well-formed core graphs; `toNamedPath_walk_vertices` preserves path support. | Byte-level `liteseq` parsing into `GFA.Document` and evidence of `Document.Accepted` are external to Lean for now. | Parser conformance and rejected-byte diagnostics belong to `lean4-conformance-harness`. |
| Flubble detection | `detectFlubbles_sound`, `detectFlubbles_complete`, `detectFlubbles_canonical_noDuplicates`, `detectFlubbles_correct`, and `detectFlubbles_correct_for_gfa` prove the Lean reference detector against `IsFlubbleBoundary`. | The current implementation must supply a `CycleClassAssignment.Correct` witness matching povu's bracket/cycle-equivalence state and canonical candidate order. | Complexity/count bounds such as the paper-level `flubble_count_le_numEdges` are not part of the current semantic emission theorem. |
| Hairpin detection | `detectHairpins_sound`, `detectHairpins_complete`, `detectHairpins_canonical_noDuplicates`, `detectHairpins_correct`, and `detectHairpins_correct_for_gfa` prove the Lean reference detector against `IsPaperHairpinBoundary`. | The implementation must supply `HairpinScanAssignment.Correct` from the reverse DFS/bracket scan. | Sequence-level proof that current inversion/SNE extraction yields the semantic `SUBR` allele calls remains conformance work. |
| Flubble hierarchy | `buildHierarchy_correct` proves the builder yields `IsCorrectHierarchy` under duplicate-free detector output and `SupportedHierarchyInput`; `buildHierarchy_rejects_unsupported` isolates non-laminar or missing-span inputs. | The integrated theorem receives `SupportedHierarchyInput` for the detected boundaries as an explicit assumption. | Proving that every povu-supported input produces laminar hierarchy inputs is not yet connected to the current implementation. |
| VCF semantics | `recordOfCall_wellFormed`, `emitRecords_ordered`, `emitRecords_correct`, and `emitVcfRecords_correct_for_gfa` prove emitted semantic records are well formed, ordered, and derived from verified flubble-tree or hairpin sources. | Semantic allele strings, traversal strings, genotype columns, and `ReferenceCallSet` membership are supplied as `VariantCall` data. | Byte-level VCF formatting, header details, decimal `AF` rendering, and current output-order normalization remain outside Lean. |
| Pipeline synthesis | `Pipeline.semanticGfaToVcf_correct` exposes the final semantic theorem boundary using the completed component interfaces. | The conformance harness must instantiate or compare the theorem assumptions against current povu behavior. | No Lean placeholder remains in the trusted path; remaining blockers are external/conformance tasks. |

## Interface Reconciliation

The architecture expected a synthesis-owned `PovuLean/Pipeline.lean` module, but
the completed component work only provided area aggregators.  This task added
`PovuLean/Pipeline.lean` and imported it from `PovuLean.lean`.  The new theorem
is a narrow wrapper around `VCF.emitVcfRecords_correct_for_gfa`, preserving the
upstream ownership boundaries while giving downstream tasks one stable theorem
name: `Pipeline.semanticGfaToVcf_correct`.

No duplicate or conflicting Lean definitions were found when importing the
core/path, GFA, flubble, hairpin, flubble-tree, VCF, and pipeline modules
together.  Existing names that intentionally overlap, such as flubble and
hairpin boundary/candidate concepts, remain separated by namespace or are
explicit `abbrev` reuses.

The hairpin proof correctly reuses the flubble candidate-stack and supported
input interface through abbreviations in `Algorithms/Hairpin/InputInvariant.lean`
instead of duplicating a competing traversal contract.

The flubble-tree layer correctly consumes `Flubble.Boundary` as an abbreviation
and adds only hierarchy-specific span, laminarity, parent, and ordering
contracts.

The VCF layer correctly consumes the verified flubble-tree and hairpin source
contracts through `VariantSource.Derived`, `ReferenceCallSet`, and
`EmissionCorrect`.  The synthesis layer did not add allele-extraction or Rust
comparison fixtures, because that would cross into the conformance task's file
ownership.

## Trusted Boundary

The trusted Lean theorem boundary is:

```lean
PovuLean.Pipeline.semanticGfaToVcf_correct
```

It proves semantic VCF emission correctness from:

- `doc : GFA.Document`;
- `accepted : doc.Accepted`;
- `frame : Core.TraversalFrame doc.toGraph`;
- `classes : Algorithms.Flubble.CycleClassAssignment doc.toGraph`;
- `scan : Algorithms.Hairpin.HairpinScanAssignment doc.toGraph`;
- `Algorithms.Flubble.SupportedGFAInput doc accepted frame`;
- `Algorithms.Flubble.CycleClassAssignment.Correct frame classes`;
- `Algorithms.Hairpin.HairpinScanAssignment.Correct frame scan`;
- `Algorithms.FlubbleTree.SupportedHierarchyInput` for the detected flubbles;
- `VCF.ReferenceCallSet` for the semantic variant calls.

The result is `VCF.EmissionCorrect` for `VCF.emitRecords calls`, which includes
emission identity, record well-formedness, derivation from verified flubble or
hairpin sources, and deterministic semantic ordering.

## Remaining Blockers for Conformance Testing

The next blocking WG task is `lean4-conformance-harness`.  It must provide or
check the bridge between current povu behavior and the Lean semantic theorem.
The concrete blockers are:

- translate actual GFA parser output into `GFA.Document` values and validate
  `Document.Accepted`;
- derive or compare traversal frames against `Core.TraversalFrame` and the
  `Flubble.SupportedGFAInput` contract;
- derive `CycleClassAssignment.Correct` from povu's bracket/cycle-equivalence
  implementation;
- derive `HairpinScanAssignment.Correct` from povu's reverse scan state;
- establish `FlubbleTree.SupportedHierarchyInput` for detected flubble spans or
  document unsupported non-laminar cases;
- extract semantic `VCF.VariantCall` data and prove/check `VCF.ReferenceCallSet`
  for flubble-tree and hairpin sources;
- compare serialized povu VCF output to `VCF.Record` semantics, including
  documented normalization for header fields, decimal formatting, and output
  ordering.

After the harness exists, `lean4-e2e-validation-corpus` should expand concrete
GFA-to-VCF fixtures.  This synthesis task intentionally did not add harness
fixtures, Rust comparison infrastructure, or end-to-end corpus files.

## Validation Record

- `lake build` succeeded with the integrated module set.
- `rg -n '\b(sorry|admit)\b' PovuLean` returned no trusted-path matches.
- The imported graph/path, GFA, flubble, hairpin, flubble-tree, VCF, and
  pipeline modules build together without duplicate-definition or conflicting
  import failures.
