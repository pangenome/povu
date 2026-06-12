#set document(
  title: "Flubble Linear Decomposition Proof Boundary",
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
  #text(size: 18pt, weight: "bold")[Flubble Linear Decomposition Proof Boundary]

  #text(size: 9pt)[Rendered proof note for `render-flubble-linear`, 2026-06-12]
]

#outline(title: [Contents], indent: auto)

= Purpose

This note is a distributable mathematical description of the current flubble
decomposition proof story in povu. It deliberately separates four layers:

+ facts formally proved in Lean,
+ conditional linear-time theorems and their stage contracts,
+ current C++/Rust bridge and conformance evidence, and
+ the remaining open gap to a full arbitrary C++ runtime refinement proof.

The formal boundary is semantic. Lean proves properties of proof-side graph,
walk, detector, hierarchy, cost, and VCF emission models. Current C++ and Rust
behavior is supported by conformance tests and bridge artifacts, but there is
not yet a Lean theorem stating that every current povu execution refines the
semantic model or satisfies the proof-side runtime-cost records.

= Mathematical Model

Let a graph be represented by a Lean `Graph` with segment list `V` and directed
link list `E`. The proof-side graph-size measure used by the cost layer is

$ n(G) = |V| + "edgeCount"(G). $

A traversal frame `F : TraversalFrame G` carries an ordered list of tree links.
The flubble detector does not scan all links directly. It scans the candidate
stack

$ S(F) = "candidateStack"(F) = [e_0, e_1, dots, e_(m-1)], $

the ordered subsequence of real black tree edges eligible to be flubble
endpoints. A cycle-class assignment `C` classifies candidate links. Its formal
correctness contract says class equality is sound and complete for the Lean
cycle-equivalence relation on boundary candidates.

For candidate `e_i`, the base detector looks for the first later same-class
candidate after a nonempty gap:

$ j = min { k | i + 1 < k < m and C(e_i) = C(e_k) }. $

If such `j` exists, the candidate pair is canonicalized by edge id:

$ b = "Boundary.ordered"(e_i, e_j). $

The decomposition relation for base flubbles is therefore:

$ "IsFlubbleBoundary"(F, C, b) <=> b " is a paper boundary and "
  b " arises from the canonical close-after-gap scan over " S(F). $

The emitted boundary list is deduplicated and deterministic. The indexed
detector introduces a class-index state so the same relation can be computed by
one pass over `S(F)`: one lookup and one update per candidate.

== Hierarchy Semantics

Each boundary has a span in the candidate stack. If

$ b = (e_i, e_j) " with " i < j, $

then `boundarySpan? S b = some { start := i, finish := j }`. A hierarchy input
is supported when every emitted boundary has a span and all spans are laminar:
for every two distinct spans, one strictly contains the other or they are
disjoint. Parent semantics is nearest strict span containment:

$ "parent"(b) = "arg min"_(p) "width"(p) $

among spans `p` that strictly contain the span for `b`; if no such span exists,
`b` is a root of the flubble forest. Lean represents the hierarchy as a flat
list of nodes with optional parent boundaries and proves it is a forest that
preserves detector order, duplicate freedom, and canonical endpoint orientation.

== Cost Model

The flubble extraction cost model is proof-side, not an observed timing model.
It has four slots:

+ `indexedOperations`,
+ `scanBookkeeping`,
+ `boundaryEmission`, and
+ `hierarchyConstruction`.

The indexed operation contract supplies constant lookup/update costs and the
exact count theorem that the indexed scan performs one lookup and one update per
candidate. Boundary output size is bounded by `m = |S(F)|`. Under supported
hierarchy input, flat hierarchy output size is bounded by `3m`: one node slot
plus two reference metadata slots per detected boundary.

Thus Lean proves the conditional statement:

$ T_"extract"(F, C) in O(|S(F)|), $

assuming the explicit operation, bookkeeping, boundary-emission, and hierarchy
construction cost contracts.

For the full semantic decomposition pipeline, Lean composes named stage costs:

+ component decomposition,
+ tip/dummy augmentation,
+ traversal-frame construction,
+ cycle-class assignment,
+ hairpin scan, and
+ flubble extraction plus hierarchy construction.

Given intermediate-size bounds back to `n(G)`, the composed theorem concludes:

$ T_"decompose"(G) in O(n(G)). $

This theorem is conditional on the named stage contracts. It is not a statement
that current C++ povu has already been formally proved linear.

= Formally Proved In Lean

The following names are Lean definitions or theorems in the `PovuLean` library.
They are the formal proof boundary this document relies on.

== Flubble Definition And Input Contracts

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

These definitions establish the semantic input model: a supported frame, a
candidate stack of eligible real black tree links, and a correct cycle-class
assignment.

== Detector Correctness And Bounds

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
+ `PovuLean.Algorithms.Flubble.detectFlubbles_correct_for_gfa`

These theorems prove that the Lean reference detector is sound and complete for
the flubble boundary relation, emits canonical duplicate-free boundaries, and
has output size bounded by the candidate stack and graph edge count.

== Indexed Detector And Equivalence

Module: `PovuLean.Algorithms.Flubble.Indexed`

Key semantic-equivalence names:

+ `PovuLean.Algorithms.Flubble.lookup_indexedEntriesFrom_eq_firstSameClass?`
+ `PovuLean.Algorithms.Flubble.indexedBoundaryAt?_eq_closeAfterGap?`
+ `PovuLean.Algorithms.Flubble.detectStackIndexedRaw_eq_detectStackRaw`
+ `PovuLean.Algorithms.Flubble.detectStackIndexed_eq_detectStack`
+ `PovuLean.Algorithms.Flubble.detectFlubblesIndexed_eq_detectFlubbles`
+ `PovuLean.Algorithms.Flubble.detectFlubblesIndexed_iff_isFlubbleBoundary`
+ `PovuLean.Algorithms.Flubble.detectFlubblesIndexed_noDuplicates`

Key operation-count names:

+ `PovuLean.Algorithms.Flubble.scanStackIndexedAux_visited`
+ `PovuLean.Algorithms.Flubble.scanStackIndexedAux_lookups`
+ `PovuLean.Algorithms.Flubble.scanStackIndexedAux_updates`
+ `PovuLean.Algorithms.Flubble.detectStackIndexedRaw_one_lookup_and_update_per_candidate`

These theorems justify replacing the reference scan by an indexed scan without
changing the emitted boundaries, while exposing the per-candidate operation
counts needed by the linear-time cost theorem.

== Hierarchy Correctness

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

These theorems prove that supported laminar boundary spans produce the canonical
flat hierarchy and that hierarchy output metadata is linear in the detector
output and, under supported input, graph edge count.

== Cost And Linearity Theorems

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
+ `PovuLean.Complexity.Flubble.DecompositionStageCosts.total_linear_with`

Module: `PovuLean.Complexity.Decomposition`

Key names:

+ `PovuLean.Complexity.Decomposition.SemanticDecomposeIntermediateSizes`
+ `PovuLean.Complexity.Decomposition.SemanticDecomposeIntermediateSizeBounds`
+ `PovuLean.Complexity.Decomposition.SemanticDecomposePipelineCosts`
+ `PovuLean.Complexity.Decomposition.SemanticDecomposeStageContracts`
+ `PovuLean.Complexity.Decomposition.conditional_semantic_decompose_linear`

The `conditional_semantic_decompose_linear` theorem composes stage linearity,
intermediate-size bounds, flubble extraction linearity, and semantic VCF
correctness into one proof-side statement.

== Semantic Pipeline Result

Module: `PovuLean.Pipeline`

Key name:

+ `PovuLean.Pipeline.semanticGfaToVcf_correct`

This theorem is the semantic GFA-to-VCF join point. Given accepted GFA input,
supported flubble input, correct cycle classes, correct hairpin scan, supported
hierarchy input, and semantic variant-call witnesses, emitted VCF records are
correct.

= Conditional Linear-Time Theorem

The strongest current full-pipeline theorem is:

```lean
PovuLean.Complexity.Decomposition.conditional_semantic_decompose_linear
```

Its conclusion is:

```lean
LinearBound costs.total sizes.inputGraphSize
  /\ VCF.EmissionCorrect frame classes scan
       (Algorithms.FlubbleTree.buildHierarchy frame classes)
       calls
       (VCF.emitRecords calls)
```

The theorem receives these stage contracts through
`SemanticDecomposeStageContracts`:

+ `componentDecompositionLinear`,
+ `tipDummyAugmentationLinear`,
+ `traversalFrameConstructionLinear`,
+ `cycleClassAssignmentLinear`,
+ `hairpinScanLinear`, and
+ `flubbleExtractionAndHierarchy`.

It also receives `SemanticDecomposeIntermediateSizeBounds`, which is the proof
that component graph size, augmented graph size, traversal-frame size,
candidate-stack size, and hairpin-boundary count are each linear in the input
graph-size measure.

The flubble portion is delegated to:

```lean
PovuLean.Complexity.Flubble.conditional_extraction_linear_in_candidateStack
```

That theorem assumes:

+ constant class-index lookup/update costs,
+ exact lookup/update counts from the indexed detector theorem,
+ scan bookkeeping linear in the candidate stack,
+ boundary emission linear in emitted boundaries,
+ hierarchy construction linear in hierarchy output size, and
+ supported laminar hierarchy input.

Therefore the correct claim is: Lean proves conditional semantic linearity under
named contracts. The incorrect claim would be: current arbitrary C++ povu
executions have already been formally proved to run in linear time.

= C++ Bridge And Conformance Evidence

The implementation bridge is evidence, not a formal refinement theorem.
Relevant documentation and surfaces are:

+ `docs/lean4-proof/bridge_cost_instrumentation.md`
+ `docs/lean4-proof/bridge_povu_implementation.md`
+ `docs/lean4-proof/checked_translator.md`
+ `src/povu/common/stage_cost.cpp`
+ `include/povu/common/stage_cost.hpp`
+ `tests/lean4_conformance/src/main.rs`
+ `src/mto/to_structure_export.cpp`
+ `povu-rs/src/native_gfa.rs`
+ `povu-rs/tests/native_gfa_vcf_tests.rs`

== Cost Instrumentation

The opt-in environment variable `POVU_STAGE_COST_TRACE=1` emits stage-counter
lines with stable `contract=` identifiers. These identifiers align implementation
trace points with the Lean stage slots:

+ `componentDecomposition`
+ `tipDummyAugmentation`
+ `traversalFrameConstruction`
+ `cycleClassAssignment`
+ `boundaryEmission`
+ `hierarchyConstruction`

The counters report calls, size proxies, outputs, and elapsed nanoseconds. They
are useful bridge artifacts, but elapsed time and container behavior are not
Lean cost witnesses.

== Flubble Export And Structure Evidence

The C++ conformance path exposes canonical structure information through
`--structure-export` in the CLI and `src/mto/to_structure_export.cpp`. The
bridge documentation describes exported stack/class/hierarchy information used
by fixture-level checks. The Rust semantic port in `povu-rs/src/native_gfa.rs`
contains `detect_flubble_stack`, with tests for close-after-gap,
first-later-same-class, canonicalization, and duplicate removal.

This is conformance evidence over maintained fixtures and port-level tests. It
does not quantify over all supported graphs.

== Checked Translator

The checked-translator artifact is produced by the Lean conformance harness:

```bash
cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- \
  --repo-root . --fixture minimal-substitution \
  --checked-translator build/lean4-conformance/minimal-substitution-witness.json
```

Its schema is `povu.lean4.checked-translator.v1`. It records the fixture,
current povu command, Lean theorem boundary, trusted assumptions, exported C++
structure, Lean reference structure, and normalized semantic VCF witness. The
artifact targets `PovuLean.Pipeline.semanticGfaToVcf_correct`; it does not
synthesize a Lean proof term.

== Cycle Class Checker

The bridge audit records a direct fixture checker for C++ class ids against an
independent cycle-equivalence oracle. It tests both directions of the Lean
`CycleClassAssignment.Correct` obligation on selected small supported structure
fixtures:

+ same exported C++ class id implies same oracle group, and
+ same oracle group implies same exported C++ class id.

This is strong finite conformance evidence. It is not a theorem for arbitrary
supported inputs.

== Hierarchy Span Checker

The bridge audit also records a PVST parent-edge checker. It reconstructs
flubble boundaries from exported stack entries using the Lean first-later
same-class nonempty-gap rule, computes stack spans, derives nearest strict
span-containment parents, and compares them with exported PVST parents.

Again, the result is fixture-level conformance evidence for exported stacks,
not a general proof of C++ hierarchy construction.

= Remaining Open Gap

A full arbitrary C++ runtime refinement proof still needs a machine-checkable
bridge from implementation executions to the Lean witnesses and cost contracts.
The missing items are:

+ parser refinement from accepted byte-level GFA records to Lean
  `GFA.Document.Accepted` and `GFA.Document.toGraph`;
+ component-decomposition refinement showing C++ components preserve Lean graph
  semantics and have total size linear in the input;
+ traversal-frame refinement exporting the C++ spanning tree, tip/dummy
  augmentation, ordered tree links, real/black flags, and candidate-stack order
  as a Lean `TraversalFrame`;
+ class refinement proving C++ bracket classes satisfy
  `CycleClassAssignment.Correct` for arbitrary supported inputs;
+ extraction refinement proving C++ `next_seen` and flubble emission match Lean
  `closeAfterGap?`, `firstSameClass?`, `Boundary.ordered`, and
  `detectFlubbles`;
+ hierarchy refinement proving C++ PVST parent edges match Lean nearest strict
  stack-span containment for arbitrary supported laminar inputs;
+ a concrete runtime cost model for loops, allocations, serialization, FFI, and
  container operations such as `std::vector`, `std::map`,
  `std::unordered_map`, `std::unordered_set`, `std::set`, `std::list`, Rust
  collections, and serialization buffers; and
+ proof that observed or exported stage counters instantiate the Lean
  `StageCost`, `ExtractionPipelineCosts`, and
  `SemanticDecomposePipelineCosts` records.

Until those obligations are discharged, the precise project status is:

#quote[
Lean proves semantic detector correctness, hierarchy correctness, indexed
detector equivalence, output-size bounds, and conditional linear-time
composition. Current povu has bridge artifacts and conformance evidence for
selected implementation surfaces. A formal arbitrary C++ runtime refinement
proof remains open.
]

= Build And Reference Check

The rendered PDF for this source is:

```text
docs/lean4-proof/flubble_linear_proof.pdf
```

Build command:

```bash
typst compile docs/lean4-proof/flubble_linear_proof.typ \
  docs/lean4-proof/flubble_linear_proof.pdf
```

The Lean names listed above are intended to be checked with a temporary
`#check` file importing:

```lean
import PovuLean.Algorithms.Flubble.Correctness
import PovuLean.Algorithms.Flubble.Indexed
import PovuLean.Algorithms.FlubbleTree.Correctness
import PovuLean.Complexity.Flubble
import PovuLean.Complexity.Decomposition
import PovuLean.Pipeline
```
