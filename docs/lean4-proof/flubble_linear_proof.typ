#set document(
  title: "Flubble Linear-Time Algorithm Proof",
  author: "povu Lean proof notes",
)
#set page(paper: "us-letter", margin: (x: 0.75in, y: 0.7in))
#set text(font: "Libertinus Serif", size: 10pt)
#set par(justify: true, leading: 0.55em)
#set heading(numbering: "1.")

#show raw.where(block: true): it => block(
  fill: rgb("#f6f8fa"),
  inset: 8pt,
  radius: 2pt,
  width: 100%,
  it,
)

#align(center)[
  #text(size: 18pt, weight: "bold")[Flubble Linear-Time Algorithm Proof]

  #text(size: 9pt)[Rendered proof note for `revise-flubble-linear`, 2026-06-12]
]

#outline(title: [Contents], indent: auto)

= Main Claim

The fast indexed flubble detector/extractor is a linear-time algorithm over the
candidate stack under the ordinary word-RAM convention that class-index lookup
and update are constant-time operations. Lean now exposes this statement as:

```lean
PovuLean.Complexity.Flubble.indexed_flubble_stage_linear_in_candidateStack
```

Lean also now exposes the Nadia-facing count theorem:

```lean
PovuLean.Algorithms.Flubble.flubble_count_le_numEdges
```

For a supported traversal frame and certified cycle-class assignment, this
states that the canonical detector output has length at most `edgeCount(G)`.
That is the formal linear-number-of-flubbles statement for the current Lean
definition: it counts deduplicated canonical `detectFlubbles` boundaries, not
all pairwise cycle-equivalent candidate-edge pairs.

The theorem covers the indexed detector, boundary emission, and flat hierarchy
construction. It uses the mechanically checked indexed detector facts that the
scan performs exactly one class-index lookup and one class-index update per
candidate, plus the checked output-size facts that emitted flubble boundaries
and hierarchy metadata are linear in the candidate stack.

Once the candidate-stack size is bounded linearly by the graph-size measure

$ n(G) = |V| + "edgeCount"(G), $

the same stage is linear in graph size. Lean exposes that corollary as:

```lean
PovuLean.Complexity.Flubble.indexed_flubble_stage_linear_in_graphSize
```

This is the algorithms-level message: after the traversal frame and
cycle-class assignment provide a candidate stack of linear size, the indexed
flubble detector/extractor and hierarchy builder add only linear work.

= Algorithm Model

Let `G` be the Lean graph, with segment list `V` and directed black/grey link
measure `edgeCount(G)`. A traversal frame `F : TraversalFrame G` contains an
ordered tree-link list. The flubble stage scans only the candidate stack

$ S(F) = "candidateStack"(F) = [e_0, e_1, dots, e_(m-1)], $

the ordered subsequence of real black tree edges eligible to be flubble
endpoints. A cycle-class assignment `C` labels candidates by the cycle
equivalence classes needed by the flubble characterization.

For candidate `e_i`, the detector asks for the first later same-class
candidate after a nonempty gap:

$ j = min { k | i + 1 < k < m and C(e_i) = C(e_k) }. $

When such `j` exists, the detector emits the canonical ordered boundary

$ b = "Boundary.ordered"(e_i, e_j). $

The reference detector defines this relation directly. The fast detector keeps
an index from class id to the next later candidate already seen in the reverse
scan. Under the standard word-RAM model, table lookup and table update are
constant-cost operations; the scan therefore does one constant amount of index
work per stack candidate.

== Hierarchy Construction

Each emitted boundary has a span in the candidate stack. If

$ b = (e_i, e_j) " with " i < j, $

then `boundarySpan? S b = some { start := i, finish := j }`. Supported
hierarchy inputs are laminar: for any two distinct spans, one strictly contains
the other or the spans are disjoint. Parenthood is nearest strict containment.
The flat hierarchy records one node per detected flubble boundary and parent /
boundary metadata slots for those nodes.

Lean proves that supported laminar inputs produce the canonical forest
hierarchy and that the flat hierarchy output size is at most three slots per
candidate-stack entry.

= Mechanically Checked Lean Facts

This section lists the theorem surface used by the algorithms claim. The names
are checked by `lake build`.

== Flubble Semantics

Module: `PovuLean.Algorithms.Flubble.Spec`

Key names:

+ `PovuLean.Algorithms.Flubble.Boundary`
+ `PovuLean.Algorithms.Flubble.Boundary.CanonicalIds`
+ `PovuLean.Algorithms.Flubble.Boundary.ordered`
+ `PovuLean.Algorithms.Flubble.CanonicalClassBoundary`
+ `PovuLean.Algorithms.Flubble.IsPaperBoundary`
+ `PovuLean.Algorithms.Flubble.IsFlubbleBoundary`
+ `PovuLean.Algorithms.Flubble.NoDuplicateBoundaries`

Module: `PovuLean.Algorithms.Flubble.InputInvariant`

Key names:

+ `PovuLean.Algorithms.Flubble.IsBoundaryCandidate`
+ `PovuLean.Algorithms.Flubble.candidateStack`
+ `PovuLean.Algorithms.Flubble.SupportedInput`
+ `PovuLean.Algorithms.Flubble.CycleClassAssignment`
+ `PovuLean.Algorithms.Flubble.CycleClassAssignment.Correct`
+ `PovuLean.Algorithms.Flubble.SupportedGFAInput`
+ `PovuLean.Algorithms.Flubble.SupportedGFAInput.toSupportedInput`

The `CycleClassAssignment.Correct` obligation is a substantive correctness
question, independent of the linear-time accounting. It states that equal class
ids are sound and complete for cycle-equivalent boundary candidates. The
flubble linear-time theorem assumes this semantic classification is already
available; it does not make the class-assignment proof disappear.

== Detector Correctness And Size Bounds

Module: `PovuLean.Algorithms.Flubble.Detect`

Key names:

+ `PovuLean.Algorithms.Flubble.detectFlubbles`
+ `PovuLean.Algorithms.Flubble.detectFlubbles_length_le_candidateStack_length`
+ `PovuLean.Algorithms.Flubble.detectFlubbles_length_le_treeLinks_length`
+ `PovuLean.Algorithms.Flubble.detectFlubbles_length_le_graph_links_length`
+ `PovuLean.Algorithms.Flubble.detectFlubbles_length_le_graph_edgeCount`
+ `PovuLean.Algorithms.Flubble.detectFlubbles_length_le_supportedInput_graph_edgeCount`
+ `PovuLean.Algorithms.Flubble.detectFlubbles_noDuplicates`

Module: `PovuLean.Algorithms.Flubble.Correctness`

Key names:

+ `PovuLean.Algorithms.Flubble.detectFlubbles_sound`
+ `PovuLean.Algorithms.Flubble.detectFlubbles_complete`
+ `PovuLean.Algorithms.Flubble.detectFlubbles_canonical_noDuplicates`
+ `PovuLean.Algorithms.Flubble.detectFlubbles_correct`
+ `PovuLean.Algorithms.Flubble.flubble_count_le_numEdges`
+ `PovuLean.Algorithms.Flubble.detectFlubbles_correct_for_gfa`

These theorems prove that the Lean reference detector is sound and complete for
the flubble-boundary relation, emits canonical duplicate-free boundaries, and
has output size bounded by the candidate stack and graph edge count. The named
`flubble_count_le_numEdges` theorem packages the graph-edge-count bound as the
paper-facing canonical flubble count theorem.

== Indexed Detector

Module: `PovuLean.Algorithms.Flubble.Indexed`

Semantic-equivalence names:

+ `PovuLean.Algorithms.Flubble.lookup_indexedEntriesFrom_eq_firstSameClass?`
+ `PovuLean.Algorithms.Flubble.indexedBoundaryAt?_eq_closeAfterGap?`
+ `PovuLean.Algorithms.Flubble.detectStackIndexedRaw_eq_detectStackRaw`
+ `PovuLean.Algorithms.Flubble.detectStackIndexed_eq_detectStack`
+ `PovuLean.Algorithms.Flubble.detectFlubblesIndexed_eq_detectFlubbles`
+ `PovuLean.Algorithms.Flubble.detectFlubblesIndexed_iff_isFlubbleBoundary`
+ `PovuLean.Algorithms.Flubble.detectFlubblesIndexed_noDuplicates`

Operation-count names:

+ `PovuLean.Algorithms.Flubble.scanStackIndexedAux_visited`
+ `PovuLean.Algorithms.Flubble.scanStackIndexedAux_lookups`
+ `PovuLean.Algorithms.Flubble.scanStackIndexedAux_updates`
+ `PovuLean.Algorithms.Flubble.detectStackIndexedRaw_one_lookup_and_update_per_candidate`

These theorems justify using the indexed scan as the fast detector: it returns
exactly the same boundaries as the reference detector and exposes the
one-lookup / one-update per-candidate accounting used by the cost theorem.

== Hierarchy Correctness And Size Bounds

Module: `PovuLean.Algorithms.FlubbleTree.Spec`

Key names:

+ `PovuLean.Algorithms.FlubbleTree.Span`
+ `PovuLean.Algorithms.FlubbleTree.boundarySpan?`
+ `PovuLean.Algorithms.FlubbleTree.SupportedHierarchyInput`
+ `PovuLean.Algorithms.FlubbleTree.UnsupportedHierarchyInput`
+ `PovuLean.Algorithms.FlubbleTree.IsForestHierarchy`
+ `PovuLean.Algorithms.FlubbleTree.IsCorrectHierarchy`

Module: `PovuLean.Algorithms.FlubbleTree.Correctness`

Key names:

+ `PovuLean.Algorithms.FlubbleTree.buildHierarchyFrom_correct`
+ `PovuLean.Algorithms.FlubbleTree.buildHierarchy_correct`
+ `PovuLean.Algorithms.FlubbleTree.buildHierarchy_rejects_unsupported`
+ `PovuLean.Algorithms.FlubbleTree.buildHierarchy_nodes_length`
+ `PovuLean.Algorithms.FlubbleTree.buildHierarchy_nodes_length_le_detectFlubbles_length`
+ `PovuLean.Algorithms.FlubbleTree.buildHierarchy_nodes_length_le_supportedInput_graph_edgeCount`
+ `PovuLean.Algorithms.FlubbleTree.buildHierarchy_flatMetadataSlots_le_supportedInput_graph_edgeCount_twice`

These theorems prove canonical hierarchy construction for supported laminar
boundary spans, and prove the hierarchy output-size bounds used by the
linear-time statement.

== Cost And Composition Theorems

Module: `PovuLean.Complexity.Flubble`

Key names:

+ `PovuLean.Complexity.Flubble.GraphSizes`
+ `PovuLean.Complexity.Flubble.DecompositionSizes`
+ `PovuLean.Complexity.Flubble.StageCost`
+ `PovuLean.Complexity.Flubble.ClassIndexContract`
+ `PovuLean.Complexity.Flubble.IndexedDetectorOperationContract`
+ `PovuLean.Complexity.Flubble.IndexedDetectorOperationContract.operationCost_linear_in_candidates`
+ `PovuLean.Complexity.Flubble.ConditionalExtractionCostContract`
+ `PovuLean.Complexity.Flubble.boundaryOutputSize_linear_in_candidateStack`
+ `PovuLean.Complexity.Flubble.hierarchyConstructionOutputSize_linear_in_candidateStack`
+ `PovuLean.Complexity.Flubble.conditional_extraction_linear_in_candidateStack`
+ `PovuLean.Complexity.Flubble.indexed_flubble_stage_linear_in_candidateStack`
+ `PovuLean.Complexity.Flubble.indexed_flubble_stage_linear_in_graphSize`
+ `PovuLean.Complexity.Flubble.DecompositionStageCosts.total_linear_with`

The two `indexed_flubble_stage_*` theorems are the short algorithms-facing
entry points. They package the indexed operation-count theorem, output-size
bounds, hierarchy-size bounds, and candidate-stack size composition.

= Full Pipeline Perspective

The complete semantic decomposition theorem is:

```lean
PovuLean.Complexity.Decomposition.conditional_semantic_decompose_linear
```

Its conclusion combines a proof-side linear cost bound with semantic VCF
correctness:

```lean
LinearBound costs.total sizes.inputGraphSize
  /\ VCF.EmissionCorrect frame classes scan
       (Algorithms.FlubbleTree.buildHierarchy frame classes)
       calls
       (VCF.emitRecords calls)
```

The theorem uses `SemanticDecomposeIntermediateSizeBounds` to compose all
intermediate measures back to the graph size. The candidate-stack bound in that
record is the bridge from the flubble-stage theorem to the whole graph-size
measure.

Five upstream full-pipeline stage linearities are textbook linear algorithmic
facts but are not mechanized in Lean yet:

+ component decomposition,
+ tip/dummy augmentation,
+ traversal-frame construction,
+ cycle-class assignment, and
+ hairpin scan.

They are standard graph/tree traversal or stack/indexing facts under the same
word-RAM model, but the current Lean theorem accepts them as named stage
contracts. This is a mechanization boundary, not evidence of an algorithmic
unknown. The separate deployed-code refinement question is whether current C++,
Rust, parser, export, and serialization paths instantiate the semantic
witnesses and cost slots for arbitrary supported inputs.

= Implementation Conformance Context

The implementation material validates and audits conformance to the algorithmic
model. It is not the main theorem.

Relevant documentation and surfaces:

+ `docs/lean4-proof/bridge_cost_instrumentation.md`
+ `docs/lean4-proof/bridge_povu_implementation.md`
+ `docs/lean4-proof/checked_translator.md`
+ `src/povu/common/stage_cost.cpp`
+ `include/povu/common/stage_cost.hpp`
+ `tests/lean4_conformance/src/main.rs`
+ `src/mto/to_structure_export.cpp`
+ `povu-rs/src/native_gfa.rs`
+ `povu-rs/tests/native_gfa_vcf_tests.rs`

== C++ Stage Counters

`POVU_STAGE_COST_TRACE=1` emits stage-counter lines with stable `contract=`
identifiers aligned with Lean stage slots:

+ `componentDecomposition`
+ `tipDummyAugmentation`
+ `traversalFrameConstruction`
+ `cycleClassAssignment`
+ `boundaryEmission`
+ `hierarchyConstruction`

These counters are useful for audits and regression checks. They are not Lean
proof objects and they do not by themselves prove a C++ runtime theorem.

== Rust And Fixture Evidence

The Rust semantic port in `povu-rs/src/native_gfa.rs` contains
`detect_flubble_stack`, with tests for close-after-gap,
first-later-same-class, canonicalization, and duplicate removal. The C++
conformance path exposes stack, class, and hierarchy information through
`--structure-export` and `src/mto/to_structure_export.cpp`.

The Lean conformance harness can produce checked-translator artifacts such as:

```bash
cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- \
  --repo-root . --fixture minimal-substitution \
  --checked-translator build/lean4-conformance/minimal-substitution-witness.json
```

These artifacts connect selected fixtures to the semantic theorem boundary
`PovuLean.Pipeline.semanticGfaToVcf_correct`. They provide finite conformance
evidence for maintained examples; they do not synthesize a Lean proof term for
all current implementation executions.

== Cycle-Class And Hierarchy Audits

The cycle-class checker compares exported C++ class ids against an independent
fixture-sized cycle-equivalence oracle in both directions: same class id implies
same oracle group, and same oracle group implies same class id. This targets
the nontrivial `CycleClassAssignment.Correct` obligation.

The hierarchy span checker reconstructs flubble boundaries from exported stack
entries using the Lean first-later same-class nonempty-gap rule, computes stack
spans, derives nearest strict span-containment parents, and compares those
parents with exported PVST parent edges.

Both checks are valuable conformance evidence. They remain fixture-level audits,
not arbitrary-input refinement proofs.

= Scope Summary

The algorithm theorem is:

#quote[
Under the standard word-RAM convention for constant-cost class-index
lookup/update, the indexed flubble detector/extractor and flat hierarchy
construction are linear in the candidate stack. Once the candidate-stack size
is linearly bounded by `|V| + edgeCount`, the flubble stage is linear in graph
size.
]

The Lean mechanization currently checks:

+ flubble boundary semantics,
+ detector soundness, completeness, canonicalization, and duplicate freedom,
+ the canonical flubble count bound `flubble_count_le_numEdges`,
+ indexed detector equivalence to the reference detector,
+ exact indexed scan lookup/update counts,
+ detector and hierarchy output-size bounds,
+ the candidate-stack and graph-size linearity corollaries for the flubble
  stage, and
+ conditional full semantic decomposition composition.

The implementation evidence currently covers:

+ C++ stage counters and structure export,
+ fixture-level C++ class and hierarchy audits,
+ Rust semantic detector tests, and
+ checked-translator artifacts for selected semantic VCF fixtures.

The deployed-code refinement scope still includes parser refinement, arbitrary
C++ traversal-frame and class-assignment witnesses, arbitrary hierarchy
refinement, serialization/FFI accounting, and a concrete runtime cost model for
the implementation containers. Those are implementation-refinement obligations,
not caveats to the standard algorithms statement above.

= Build And Reference Check

The rendered PDF for this source is a generated build artifact. Render for
validation to an output path outside git:

```text
/tmp/flubble_linear_proof.pdf
```

Build command:

```bash
typst compile docs/lean4-proof/flubble_linear_proof.typ \
  /tmp/flubble_linear_proof.pdf
```

The Lean names listed in this document are checked by `lake build`. A focused
temporary `#check` file may import:

```lean
import PovuLean.Algorithms.Flubble.Correctness
import PovuLean.Algorithms.Flubble.Indexed
import PovuLean.Algorithms.FlubbleTree.Correctness
import PovuLean.Complexity.Flubble
import PovuLean.Complexity.Decomposition
import PovuLean.Pipeline
```
