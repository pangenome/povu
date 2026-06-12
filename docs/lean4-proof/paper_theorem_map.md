# Paper Theorem Map for Lean4 Formalization

Task: `lean4-paper-map`

Paper analyzed: Njagi Mwaniki, Erik Garrison, Nadia Pisanti, "Popping
Bubbles in Pangenome Graphs", arXiv:2410.20932v1.

Version/date used: arXiv v1, submitted 2024-10-28 11:26:06 UTC. The
analysis uses the arXiv abstract page and the v1 TeX/PDF source available from
https://arxiv.org/abs/2410.20932 and https://arxiv.org/e-print/2410.20932.

Scope note: this document is a theorem/definition map only. It does not create
Lean modules or assert final names for the architecture task. Candidate Lean
names below are stable handles for downstream discussion.

## Executive Summary

The paper gives a compact proof outline rather than a fully formal sequence of
lemmas. For Lean, the formalization should separate four layers:

1. A finite biedged/variation graph model with edge identities, black/grey edge
   colors, valid alternating walks, labels, connected components, tips, and
   compactness.
2. A DFS-rooted intermediate structure `T'_s`: a spanning tree plus non-tree
   back edges, ancestor order, brackets, and cycle-equivalence classes of tree
   edges.
3. Algorithmic correctness for reported structures: hairpin inversions,
   flubble boundary pairs, and the flubble forest/tree.
4. Resource proofs: node/edge count bounds, output-size bounds, and a cost
   model showing linear time and space.

The highest-risk proof gaps are:

- the paper's "back-edges do not overlap" property needs a precise statement
  and may require stronger graph/DFS assumptions than the paper states;
- cycle equivalence must be defined over either all graph cycles or only valid
  alternating biedged cycles;
- the "suitably modified" `T'_s` with virtual backedges for separating tips,
  hairpins, and flubbles is not specified in enough detail for Lean;
- the hairpin observation is necessary, not by itself a full soundness and
  completeness proof for a detector;
- canonical selection of exactly two black boundary edges from a cycle
  equivalence class is under-specified.

## Paper Claims Coverage Matrix

| Paper claim or definition | Candidate Lean obligation | Status for Lean |
| --- | --- | --- |
| Graphs, walks, cycles, spanning trees, DFS trees | `Walk`, `Cycle`, `SimpleCycle`, `SpanningTree`, `DFSSpanningTree` | In scope; foundational definitions. |
| Bidirected, sequence, biedged, and variation graphs | `BidirectedGraph`, `BiedgedGraph`, `VariationGraph`, `ValidBiedgedWalk` | In scope; requires a concrete edge-id representation. |
| Edge and node cycle equivalence; CEEP/CEVP | `EdgeCycleEquivalent`, `NodeCycleEquivalent`, `edgeCycleEquivalent_equivalence`, `black_ceep_iff_bidirected_cevp` | In scope; core equivalence layer. |
| `k`-edge and `k`-vertex disconnectability | `KEdgeDisconnectable`, `KVertexDisconnectable`, `KBlackEdgeDisconnectable` | In scope; black-edge variant must be added explicitly. |
| Compact input assumption and component decomposition | `CompactVariationGraph`, `compact_collapse_preserves_obligations`, `components_independent_deconstruct` | In scope; either prove preprocessing or make compactness an input invariant. |
| Tip dummy-node augmentation | `Tip`, `addDummyRoot`, `dummyRoot_connected`, `dummyRoot_preserves_real_flubbles` | In scope; preservation proof needed. |
| DFS ancestor order and non-tree edges are ancestor/descendant | `dfsAncestor_partialOrder`, `dfs_nonTreeEdge_comparable` | In scope, with assumptions about undirected traversal of the biedged carrier. |
| Bracket sets characterize cycle equivalence in `T'_s` | `treeEdge_cycleEquivalent_iff_bracketSet_eq` | In scope; central theorem, proof not supplied in detail by paper. |
| Hairpin definition and hairpin observation | `IsHairpin`, `HairpinLoop`, `HairpinStem`, `hairpin_has_unbracketed_treeEdge` | In scope; paper states only the necessary observation. |
| Hairpin detector reports exactly hairpins | `detectHairpins_sound`, `detectHairpins_complete` | In scope for implementation correctness, but not fully proven in the paper. |
| Flubble as cycle-equivalent black tree-edge boundary pair | `IsFlubbleBoundary`, `Flubble`, `flubble_boundary_is_ceep` | In scope; requires canonical boundary convention. |
| Reverse-DFS boundary detection | `reverseDfs_flubbleBoundaries_sound`, `reverseDfs_flubbleBoundaries_complete` | In scope; algorithm proof obligation. |
| At most `|E|` flubbles | `flubble_count_le_numEdges` | In scope; maps to paper Theorem 1. |
| Flubble tree/forest is computable in linear time | `buildFlubbleForest_linear_time`, `buildFlubbleForest_linear_space` | In scope; maps to paper Corollary 1 plus cost model. |
| Leaves are minimal flubbles; depth/width interpretation | `flubbleTree_leaf_iff_minimal`, `flubbleTree_parent_contains_child` | Leaf/minimal and nesting are in scope. Biological interpretation is out of scope. |
| Post-order popping/decomposition applications | `popFlubble_preserves_language` or similar | Out of scope for this batch unless downstream tasks choose to formalize decomposition. |
| Experimental runtime comparisons with `vg` and `BubbleGun` | None | Out of scope for Lean proof; empirical claim only. |

## Foundational Graph Definitions

### Plain Finite Graph

Paper source: Section 2.

Candidate Lean definitions:

- `structure FiniteGraph`
  - finite vertex type or finite set of vertex ids;
  - finite edge id type or finite set of edge ids;
  - endpoint function from edge id to one or two vertices, depending on the
    directed/undirected carrier chosen by the architecture.
- `Walk G a b`
- `ValidWalk G w`
- `Cycle G c`
- `SimpleCycle G c`
- `SpanningTree G T`
- `DFSSpanningTree G root T`

Statement sketches:

- `walk_edges_incident`: every consecutive vertex pair in a walk is joined by
  an edge.
- `cycle_is_closed_walk`: a cycle is a walk whose first and last vertices
  coincide.
- `simpleCycle_no_internal_repeats`: a simple cycle has no repeated vertices
  except the closing vertex.
- `spanningTree_connected_acyclic`: a spanning tree contains all vertices and
  is connected/acyclic.

Lean-risk notes:

- The paper starts with edges as pairs, but GFA/biedged graphs need edge
  identities. Lean should not model `E` as a set of endpoint pairs if the
  implementation may need parallel grey edges, edge labels, or stable output
  ids.
- The paper's cycle definition says "path" but earlier defines only walks. Lean
  must choose whether cycles are valid walks, edge-simple walks, or
  vertex-simple cycles depending on the theorem.

### Bidirected, Sequence, Biedged, and Variation Graphs

Paper source: Section 2 and Section 6.

Candidate Lean definitions:

- `Side := Head | Tail`
- `BidirectedVertexSide := VertexId x Side`
- `structure BidirectedGraph`
- `structure BiedgedGraph`
  - vertices are side nodes;
  - every edge has a color `Black` or `Grey`;
  - each vertex is incident to at most one black edge;
  - black edges pair opposite sides of an original segment/node;
  - labels live on black edges for the GFA/biedged view.
- `VariationGraph := BiedgedGraph` plus any GFA/sequence-label invariants.
- `ValidBiedgedWalk G w`: consecutive edges alternate black and grey.
- `collapseBiedgedToBidirected`
- `expandBidirectedToBiedged`

Statement sketches:

- `black_edge_unique_at_vertex`: each side node has zero or one incident black
  edge.
- `validBiedgedWalk_alternates`: no two adjacent traversed edges have the same
  color.
- `biedged_bidirected_roundtrip`: expanding then collapsing preserves the
  bidirected graph up to isomorphism, under the chosen edge-id convention.
- `gfa_segments_are_black_edges`: GFA segment lines map to labelled black
  edges.
- `gfa_links_are_grey_edges`: GFA link lines map to grey edges.

Lean-risk notes:

- The paper alternates between vertex labels in sequence graphs and black-edge
  labels in biedged graphs. For povu/GFA conformance, labels should probably
  live on black edge ids.
- The paper does not say whether cycle equivalence in a biedged graph ranges
  over all carrier cycles or only valid alternating cycles. For genomic
  semantics, valid alternating cycles are the meaningful objects; for graph DFS
  theorems, the undirected carrier may be easier. This choice must be explicit.
- If Lean uses both a carrier graph and a valid-walk predicate, every theorem
  must state which layer it quantifies over.

## Cycle Equivalence and Disconnectability

### Edge and Node Cycle Equivalence

Paper source: Definition 1.

Candidate Lean definitions:

- `EdgeCycleEquivalent G e1 e2`
- `NodeCycleEquivalent G v1 v2`
- `CEEP G e1 e2 := EdgeCycleEquivalent G e1 e2`
- `CEVP G v1 v2 := NodeCycleEquivalent G v1 v2`
- `CycleEquivalenceClass G e`

Statement sketches:

- `edgeCycleEquivalent_refl`: every edge is cycle equivalent to itself.
- `edgeCycleEquivalent_symm`: cycle equivalence is symmetric.
- `edgeCycleEquivalent_trans`: cycle equivalence is transitive.
- `edgeCycleEquivalent_equivalence`: edge cycle equivalence is an equivalence
  relation.
- `nodeCycleEquivalent_equivalence`: node cycle equivalence is an equivalence
  relation.
- `black_ceep_iff_bidirected_cevp`: for a black edge connecting opposite sides
  in the biedged graph, edge cycle equivalence corresponds to node cycle
  equivalence after collapsing to the bidirected graph.

Dependencies:

- finite graph/cycle definitions;
- choice of valid-cycle semantics;
- biedged/bidirected isomorphism.

Proof-risk notes:

- The paper's footnote says cycle-equivalence classes are indeed equivalence
  classes. Lean must prove this, including transitivity. The proof is direct
  from the "same set of cycles contains both" characterization.
- The CEEP/CEVP correspondence is only a remark in the paper; it becomes a
  separate theorem if downstream proofs move between bidirected and biedged
  representations.

### Disconnectability

Paper source: Definition 2 and Section 4.1/4.2.

Candidate Lean definitions:

- `DisconnectsSubgraph G cut H`
- `KEdgeDisconnectable G H k`
- `KVertexDisconnectable G H k`
- `KBlackEdgeDisconnectable G H k`
- `TwoBlackEdgeConnectedSubgraph G H`

Statement sketches:

- `kEdgeDisconnectable_minimal`: if `KEdgeDisconnectable G H k`, then some
  `k`-edge cut disconnects `H`, and no smaller cut does.
- `kBlackEdgeDisconnectable_ignores_grey`: grey edges are not removable in the
  black-edge cut variant.
- `twoBlackEdgeDisconnectable_boundaryPair`: a proper subgraph with exactly two
  black boundary edges is 2-black-edge-disconnectable under the chosen cut
  convention.

Dependencies:

- subgraph representation;
- graph deletion/removal semantics;
- connectedness of the remaining carrier graph.

Proof-risk notes:

- The paper defines general `k`-edge-disconnectability, but later relies on
  "blackedge" variants. Lean should not overload these: black-edge cuts should
  be a distinct definition.
- Need a convention for whether the boundary edges themselves belong to the
  subgraph, to the complement, or to a boundary relation between them. This
  affects hairpin and flubble boundary theorems.

## povu Intermediate Structure

Paper source: Section 4.

### Compactness and Connected Components

Candidate Lean definitions:

- `CompactVariationGraph G`
- `LinearChainOfBlackEdges G chain`
- `collapseLinearBlackChain G chain`
- `ConnectedComponent G C`
- `deconstructComponents G`

Statement sketches:

- `compact_collapse_preserves_valid_walk_labels`: collapsing a linear chain of
  black edges preserves spelled sequences of valid walks.
- `compact_collapse_preserves_cycle_equivalence`: after quotienting a linear
  black chain, cycle-equivalence obligations transport across the quotient.
- `components_independent_flubbles`: flubbles/hairpins of a disconnected input
  are the disjoint union of those in its connected components.
- `components_independent_complexity`: total work across components is bounded
  by the work on the original graph.

Dependencies:

- graph quotient/isomorphism;
- valid-walk label semantics;
- finite component decomposition.

Proof-risk notes:

- The paper says compactness is assumed without loss of generality because
  linear chains can be collapsed. Lean can either prove the collapse preserves
  the relevant language and structures, or make compactness a precondition of
  all later theorems. The first path is stronger but increases scope.

### Tips, Dummy Root, and Augmented Input

Candidate Lean definitions:

- `Tip G v`: a side node not incident with a grey edge.
- `addDummyRootForTips G : BiedgedGraph`
- `DummyRoot s`
- `RealEdge G' e`: edge was present before augmentation.
- `VirtualEdge G' e`: dummy or virtual edge added for algorithmic separation.

Statement sketches:

- `dummyRoot_connected`: after adding a dummy root and dummy grey edges to all
  tips, the augmented component used for DFS is connected.
- `dummyRoot_noTips_or_arbitraryRoot`: if there are no tips, any chosen real
  node can serve as DFS root.
- `dummyRoot_preserves_real_cycles`: cycles using only real edges are preserved
  by augmentation.
- `dummyRoot_preserves_real_flubbles`: reported flubbles over real black edges
  are preserved after removing dummy artifacts, provided virtual edges are
  filtered.
- `dummyRoot_tip_artifacts_filtered`: structures whose only support is a tip
  augmentation are not reported as real flubbles.

Dependencies:

- tip definition;
- edge provenance;
- connectedness after augmentation.

Proof-risk notes:

- The paper later mentions adding further virtual backedges to distinguish
  tips, hairpins, and flubbles, but does not define this construction. Lean
  needs an explicit `VirtualEdge` provenance discipline before correctness can
  be stated.

### DFS Tree, Ancestor Order, Back-Edges, and Brackets

Candidate Lean definitions:

- `DFSTree G root T`
- `Ancestor T u v`
- `TreeLt T u v`: `u` and `v` lie on one root-leaf path and `u` is closer to
  the root.
- `BackEdge G T e`: an edge in `G` but not in `T`, with endpoints comparable
  by `TreeLt`.
- `AugmentedTree G T`: `T'_s`, the DFS tree plus back edges and ancestor order.
- `BracketOf T backEdge treeEdge`
- `BracketSet T e : Finset EdgeId`

Statement sketches:

- `treeLt_irrefl`: no vertex is closer to the root than itself.
- `treeLt_trans`: ancestor order is transitive along a root-leaf path.
- `treeLt_antisymm`: two distinct vertices cannot be ancestors of each other.
- `dfsAncestor_partialOrder`: `TreeLt` plus equality gives a partial order.
- `dfs_nonTreeEdge_comparable`: every non-tree edge in the DFS carrier connects
  ancestor/descendant endpoints.
- `backEdge_creates_fundamentalCycle`: every back edge plus the tree path
  between its endpoints forms a cycle.
- `treeEdge_in_fundamentalCycle_iff_bracket`: a tree edge lies on the
  fundamental cycle of a back edge exactly when the back edge brackets it.
- `treeEdge_cycleEquivalent_iff_bracketSet_eq`: two tree edges are cycle
  equivalent in `T'_s` iff their bracket sets are equal.

Dependencies:

- DFS theorem for the graph type chosen by architecture;
- unique tree path;
- cycle/fundamental-cycle construction;
- finite sets of back edges.

Proof-risk notes:

- The theorem `dfs_nonTreeEdge_comparable` is standard for undirected DFS. The
  paper's graph model includes oriented sides, but its algorithm appears to use
  the underlying biedged carrier. Lean should state the carrier explicitly.
- The paper says DFS ensures back-edges do not overlap. A Lean theorem needs a
  formal predicate such as `LaminarBackEdges T backEdges`. This property is not
  true for every arbitrary undirected graph without qualification, so downstream
  tasks must either prove the needed laminarity from extra biedged/compact
  invariants or weaken the algorithm proof to avoid relying on that phrase.
- `treeEdge_cycleEquivalent_iff_bracketSet_eq` is the central bridge from graph
  cycles to the reverse-DFS algorithm. It should be proven before flubble
  correctness.

## Hairpin Inversions

Paper source: Section 4.1 and experimental Section 6.

### Definitions

Candidate Lean definitions:

- `IsHairpin G H`
- `HairpinLoop G H loop`
- `HairpinStem G H stem`
- `HairpinBoundary G H e1 e2`
- `HairpinInversion G H`: optional biological interpretation layer.

Statement sketch:

```lean
def IsHairpin (G : BiedgedGraph) (H : Subgraph G) : Prop :=
  KBlackEdgeDisconnectable G H 1
  /\ ExistsSimpleCycleIn H
  /\ ExistsLoopStemDecomposition G H
```

The loop/stem decomposition should include:

- loop is a 1-black-edge-disconnectable simple cycle;
- stem is cycle-free and 2-black-edge-disconnectable, or a single connector
  edge;
- stem attaches the loop to the rest of the graph through the chosen boundary
  convention.

### Paper Observation 1

Paper claim:

- A hairpin in `G` remains a hairpin in `T'_s`.
- If `G` contains a hairpin, at least one tree edge lacks a bracket in `T'_s`.

Candidate Lean obligations:

- `hairpin_lifts_to_augmentedTree`
  - If `IsHairpin G H` and `T'` is the augmented tree view of the same carrier,
    then `IsHairpin T' H`.
- `hairpin_has_unbracketed_treeEdge`
  - If `IsHairpin G H`, `T` is the DFS tree used by povu, and `T'` is the
    augmented tree, then there exists a tree edge in the stem or connector
    whose `BracketSet T e` is empty.

Extra algorithm obligations not fully supplied by the paper:

- `detectHairpins_complete`
  - Every real hairpin satisfying the formal input assumptions is reported by
    the detector.
- `detectHairpins_sound`
  - Every reported real boundary pair corresponds to an `IsHairpin` subgraph.
- `reportedHairpin_boundaries_are_stem_ends`
  - The two output black-edge ids are exactly the stem boundary edges under the
    chosen convention.
- `hairpin_virtualBackedges_separate_from_flubbles`
  - Virtual backedges added for classification prevent hairpins from being
    emitted as flubbles.

Dependencies:

- black-edge disconnectability;
- DFS/bracket sets;
- bridge/articulation facts if the proof follows the paper footnote;
- virtual-edge filtering.

Proof-risk notes:

- The paper's observation is one-way. An unbracketed tree edge may indicate a
  bridge/tip artifact rather than a hairpin. Lean needs additional predicates
  for the detector's candidate construction.
- The paper treats hairpins as inversions biologically, but a formal inversion
  theorem would require sequence labels and reverse-complement semantics. That
  is not necessary for flubble-tree correctness unless downstream tasks want to
  certify the biological interpretation.

## Flubbles

Paper source: Section 4.2.

### Definitions and Boundary Selection

Candidate Lean definitions:

- `CycleEquivalenceClassInAugmentedTree G T e`
- `OpensCycleClass T e classId`
- `ClosesCycleClass T e classId`
- `IsFlubbleBoundary G T e1 e2`
- `structure Flubble`
  - `leftBoundary : EdgeId`
  - `rightBoundary : EdgeId`
  - proof boundaries are distinct real black tree edges;
  - proof boundaries are cycle equivalent in `T'_s`;
  - proof boundaries are canonical open/close edges for one class.

Statement sketch:

```lean
def IsFlubbleBoundary (G : BiedgedGraph) (T : DFSTree G root)
    (e1 e2 : EdgeId) : Prop :=
  e1 <> e2
  /\ IsRealBlackTreeEdge G T e1
  /\ IsRealBlackTreeEdge G T e2
  /\ EdgeCycleEquivalent (AugmentedTree G T) e1 e2
  /\ CanonicalClassBoundary G T e1 e2
```

Candidate obligations:

- `flubbleBoundary_same_rootLeaf_path`
  - Boundary edges of a flubble lie on the same root-leaf path of `T`.
- `flubbleBoundary_has_multiple_paths`
  - There are at least two valid paths connecting the boundary context, using
    the paper's intended notion of "multiple paths".
- `flubbleBoundary_twoBlackEdgeDisconnectable`
  - The boundary pair makes the represented subgraph 2-black-edge-
    disconnectable.
- `flubble_boundary_is_ceep`
  - Every flubble boundary pair is a CEEP in the real graph or in the
    augmented-tree graph, depending on the chosen theorem target.
- `cycleClass_boundary_unique`
  - Each cycle-equivalence class has a unique canonical open/close boundary
    pair after tie-breaking.

Dependencies:

- bracket-set characterization;
- canonical ordering of tree edges by reverse DFS;
- black tree-edge filter;
- compact graph assumption.

Proof-risk notes:

- The paper identifies a flubble with a pair of black tree edges, but cycle
  equivalence classes can contain more than two edges. Lean needs the
  "open/close" rule as an explicit canonicalization theorem.
- The paper says compactness and DFS imply multiple paths and same root-leaf
  path. The same-root-leaf part follows from ancestor comparability; the
  multiple-path claim needs a graph-theoretic construction.
- The statement "each flubble is a CEEP" may be definitional if flubbles are
  defined by cycle equivalence, but it is still useful as a named theorem for
  downstream modules.

### Reverse-DFS Detection

Paper claim:

- Flubble boundaries are detected in linear time by traversing the DFS tree in
  reverse DFS order and taking the two black edges that open and close a cycle
  equivalence class.

Candidate Lean obligations:

- `reverseDfs_classState_invariant`
  - During reverse DFS, the maintained class state equals the bracket-set
    equivalence information for the processed subtree boundary.
- `reverseDfs_detects_openClose_edges`
  - The algorithm's emitted pair for a class is exactly the canonical open and
    close pair.
- `reverseDfs_flubbleBoundaries_sound`
  - Every emitted pair satisfies `IsFlubbleBoundary`.
- `reverseDfs_flubbleBoundaries_complete`
  - Every `IsFlubbleBoundary` pair is emitted once.
- `reverseDfs_no_duplicate_flubbles`
  - No canonical flubble is emitted more than once.

Dependencies:

- finite DFS postorder/reverse-preorder traversal;
- bracket-set equality theorem;
- class-state data structure specification;
- virtual-edge filtering.

Proof-risk notes:

- The paper does not give pseudocode for the reverse-DFS state. If the Rust
  implementation has a concrete state machine, the Lean reference should
  formalize that state machine directly rather than reconstructing it from the
  prose.

## Flubble Tree and Forest

Paper source: Sections 4.2, 5, and 6.

### Definitions

Candidate Lean definitions:

- `FlubbleContains G Fparent Fchild`
- `StrictFlubbleContains`
- `ImmediateFlubbleParent`
- `FlubbleTree G component`
- `FlubbleForest G`
- `MinimalFlubble G F`

Suggested containment convention:

- Interpret a flubble as an interval on a DFS root-leaf path plus the subgraph
  induced/covered by bracket cycles between the two boundary black edges.
- `Fparent` contains `Fchild` when every real edge/vertex in the child's
  represented subgraph belongs to the parent's represented subgraph, or
  equivalently when the child's DFS interval is strictly nested inside the
  parent's interval. Lean should prove equivalence if both views are used.

Candidate obligations:

- `flubbleContainment_irrefl`
  - Strict containment is irreflexive.
- `flubbleContainment_trans`
  - Strict containment is transitive.
- `flubbleContainment_laminar`
  - Two flubbles are either disjoint, equal, or one contains the other. This is
    the key no-crossing property needed for a tree.
- `immediateParent_unique`
  - Every non-root flubble has a unique minimal strict container.
- `buildFlubbleTree_nodes_exact`
  - The tree nodes are exactly the real flubbles in a connected component plus
    any explicit dummy root node chosen by the implementation.
- `buildFlubbleTree_edges_are_immediateParents`
  - Tree edges encode immediate containment.
- `flubbleTree_is_tree`
  - For a connected component, the parent relation forms a rooted tree.
- `flubbleForest_components`
  - For a disconnected input, deconstruct outputs one tree per connected
    component.
- `flubbleTree_leaf_iff_minimal`
  - Leaves correspond to flubbles that contain no other flubble.
- `siblingFlubbles_are_consecutive_or_disjoint`
  - Optional theorem capturing the paper's statement about chains of minimal
    flubbles appearing as sibling leaves.

Dependencies:

- flubble boundary canonicalization;
- laminarity/no-overlap invariant;
- component decomposition;
- finite rooted tree definition.

Proof-risk notes:

- The paper invokes the Program Structure Tree analogy but does not define a
  full containment relation. Lean must define containment before it can prove a
  tree.
- The example tree has a dummy/root node `D`; architecture should decide
  whether `D` is part of `FlubbleTree` or only an output artifact.
- The paper states topological interpretation claims about depth/width and
  biological structure. Those are useful documentation, but not formal proof
  obligations unless a downstream task defines quantitative tree statistics.

## Count and Complexity Claims

Paper source: Abstract, Theorem 1, Corollary 1, Section 6.

### Theorem 1: At Most `|E|` Flubbles

Paper claim:

- In a biedged graph `G = (V, E)`, there are at most `|E|` flubbles.

Lean detector-surface theorem:

```lean
theorem PovuLean.Algorithms.Flubble.flubble_count_le_numEdges
    {g : Graph} {frame : TraversalFrame g}
    {classes : CycleClassAssignment g}
    (hInput : SupportedInput g frame)
    (hClasses : CycleClassAssignment.Correct frame classes) :
    (detectFlubbles frame classes).length <= g.edgeCount
```

This is intentionally stated over the existing canonical detector list:
`detectFlubbles` emits deduplicated canonical `IsFlubbleBoundary` values for a
supported traversal frame and certified cycle-class assignment.  It is not a
graph-only count of all pairwise cycle-equivalent edge pairs.

Likely proof dependencies:

- `backEdges_card_lt_treeEdges_card`
  - In the connected augmented DFS structure, the number of relevant back edges
    is bounded by the number of tree edges under the paper's assumptions.
- `flubbles_card_le_backEdges_card_add_one`
  - Each flubble is charged to one back edge, possibly with one extra root/top
    class.
- `treeEdges_card_le_edges_card`
  - Tree edges are a subset of graph edges.
- `backEdges_add_one_le_edges`
  - Handles the paper's "possibly +1" phrase.

Proof-risk notes:

- The exact theorem should probably be `<= |real black edges|` or `<= |E|`,
  depending on output semantics. The paper says `|E|`.
- The phrase "possibly +1" needs a precise case split. For a connected graph
  with at least one vertex, `#backEdges + 1 <= #edges` is true when there is at
  least one tree edge; degenerate one-vertex/no-edge components need an
  explicit convention.

### Corollary 1: Linear Flubble Tree Construction

Paper claim:

- The flubble tree can be computed in linear time in the size of `G`.

Candidate Lean theorems:

- `buildFlubbleForest_output_size_linear`
  - The number of output nodes/edges is `O(|E|)`.
- `buildFlubbleForest_linear_time`
  - Under a RAM/list-array cost model, the build procedure takes `O(|V| + |E|)`.
- `buildFlubbleForest_linear_space`
  - The auxiliary storage and output storage are `O(|V| + |E|)`.
- `deconstruct_linear_time`
  - Full workflow, including component decomposition, tip augmentation, DFS,
    cycle class assignment, flubble tree construction, and hairpin scan, is
    linear.

Suggested cost-model decomposition:

| Stage | Candidate theorem | Proof idea |
| --- | --- | --- |
| GFA-to-biedged graph already parsed | `biedgedGraph_size_eq_gfa_records` | Out of scope unless GFA parser is formalized. |
| Component decomposition | `connectedComponents_linear` | Each vertex/edge visited a constant number of times. |
| Tip detection and dummy root edges | `tipAugmentation_linear` | One incidence scan plus one edge per tip. |
| DFS spanning tree | `dfsTree_linear` | Standard DFS over finite adjacency arrays. |
| Bracket/class assignment | `cycleClassAssignment_linear` | Each tree/back edge contributes constant or amortized-constant updates. |
| Hairpin scan | `hairpinDetection_linear` | Traverse tree/class state once. |
| Flubble boundary extraction | `reverseDfs_flubbleExtraction_linear` | Reverse DFS visits each tree edge once. |
| Tree/forest parent construction | `flubbleTreeBuild_linear` | Requires laminar intervals or stack discipline. |

Proof-risk notes:

- A pure mathematical existence proof is not enough for linear time. Lean needs
  either a verified executable reference algorithm with a cost semantics or a
  clear separation where complexity is documented but not machine-checked.
- If using Lean as a conformance oracle rather than a high-performance
  implementation, downstream tasks may choose to prove count/output-size
  bounds first and leave detailed runtime proofs as a later phase.

## Out-of-Scope or Documentation-Only Paper Claims

These paper claims should be mentioned by downstream docs but do not need Lean
obligations in the first proof roadmap:

- Historical comparison with bubbles, superbubbles, ultrabubbles, mouths, and
  snarls, except for definitions needed to state differences.
- Biological claim that hairpins correspond to inversions, unless reverse-
  complement label semantics are in scope.
- Experimental runtime comparisons against `vg` and `BubbleGun`.
- Claim that no other tool detects hairpin inversions in pangenome graphs.
- Proposed future use of the flubble tree for progressive graph popping and
  hybrid graph/ED-string representations, unless a future decomposition task
  formalizes a language-preservation theorem.

Potential future obligations if decomposition is later formalized:

- `popMinimalFlubble_preserves_spelledLanguage`
- `postorderPopping_terminates`
- `postorderPopping_preserves_componentInterfaces`
- `flubbleTree_guides_valid_decomposition`

## Proof Dependency Order

Recommended dependency order for Lean tasks:

1. `FiniteGraph`, edge ids, subgraphs, deletion, connectedness, walks, cycles.
2. Biedged graph invariants, black/grey edge colors, valid alternating walks,
   labels, tips, and GFA-facing size definitions.
3. Cycle equivalence and disconnectability definitions, including equivalence
   relation proofs and black-edge cut variants.
4. DFS tree infrastructure: tree paths, ancestor order, non-tree/back edges,
   fundamental cycles, bracket sets.
5. Bracket characterization:
   `treeEdge_cycleEquivalent_iff_bracketSet_eq`.
6. Tip/dummy/virtual-edge augmentation preservation and filtering.
7. Hairpin definitions and detector soundness/completeness.
8. Flubble boundary definition, canonical open/close rule, reverse-DFS
   soundness/completeness.
9. Flubble containment, laminarity, parent uniqueness, tree/forest correctness.
10. Count bounds and then time/space bounds.

## Open Assumptions and Gaps for Downstream Tasks

1. Cycle semantics in biedged graphs:
   Decide whether cycle equivalence quantifies over all carrier cycles or only
   valid alternating biedged cycles. The algorithmic DFS may need carrier
   cycles, while genomic correctness likely needs valid cycles.

2. Edge identity:
   Use edge ids rather than endpoint-pair sets. This avoids ambiguity around
   grey links, labels, provenance, output boundary ids, and possible parallel
   links.

3. Black-edge disconnectability:
   Define it separately from ordinary edge disconnectability. State whether
   grey edges are never removable, whether boundary black edges are part of the
   cut, and how the disconnected subgraph is represented.

4. Compactness:
   Either make `CompactVariationGraph` a precondition for flubble theorems or
   prove that collapsing linear chains preserves all relevant outputs. The
   paper's "without loss of generality" is not enough for Lean.

5. DFS carrier:
   Specify the undirected graph on which DFS runs. If it is the biedged carrier
   with colored edges ignored for traversal, say so in every DFS theorem.

6. Back-edge laminarity/non-overlap:
   The paper relies on non-overlap language without a formal definition. This
   is the likely hardest mathematical gap. Downstream tasks need either:
   - a precise laminar interval theorem under povu-specific graph invariants; or
   - an implementation proof that does not require global laminarity.

7. Bracket-set equality:
   Prove `treeEdge_cycleEquivalent_iff_bracketSet_eq` before depending on
   cycle classes for flubbles. This is the central correctness bridge.

8. Virtual backedges:
   The paper says flubble-tree detection uses a modified `T'_s` with additional
   virtual backedges to distinguish hairpins and tips from flubbles. The exact
   construction must be recovered from the implementation or specified in the
   architecture document.

9. Hairpin detector completeness:
   Observation 1 is not a full detector proof. The Lean plan needs soundness
   and completeness statements for the actual reported boundary pairs.

10. Flubble boundary canonicalization:
    Define how a cycle-equivalence class maps to exactly two black tree edges,
    especially when a class contains more than two eligible edges or virtual
    artifacts.

11. Root and dummy flubble-tree node:
    Decide whether the root `D` shown in the paper is a real tree node in Lean
    or an output-only sentinel. This affects parent uniqueness and leaf/minimal
    theorems.

12. Degenerate components:
    State behavior for empty graphs, single-node components, all-tip
    components, self-loops, and components with no flubbles.

13. Complexity proof target:
    Decide whether the first Lean milestone proves only mathematical
    output-size bounds or also verifies executable linear-time algorithms under
    a cost model.

14. Paper typo/wording:
    Theorem 1 in the PDF says "at most a `|E|` flubbles"; read this as at most
    `|E|` flubbles. The Lean theorem should use the unambiguous bound.

## Suggested Lean Names Index

Foundational:

- `BiedgedGraph`
- `VariationGraph`
- `ValidBiedgedWalk`
- `Cycle`
- `SimpleCycle`
- `SpanningTree`
- `DFSSpanningTree`
- `Tip`
- `CompactVariationGraph`

Cycle/disconnectability:

- `EdgeCycleEquivalent`
- `NodeCycleEquivalent`
- `CEEP`
- `CEVP`
- `edgeCycleEquivalent_equivalence`
- `nodeCycleEquivalent_equivalence`
- `KBlackEdgeDisconnectable`
- `TwoBlackEdgeDisconnectable`

DFS/brackets:

- `TreeLt`
- `BackEdge`
- `AugmentedTree`
- `BracketOf`
- `BracketSet`
- `dfs_nonTreeEdge_comparable`
- `backEdge_creates_fundamentalCycle`
- `treeEdge_in_fundamentalCycle_iff_bracket`
- `treeEdge_cycleEquivalent_iff_bracketSet_eq`

Hairpins:

- `IsHairpin`
- `HairpinLoop`
- `HairpinStem`
- `HairpinBoundary`
- `hairpin_lifts_to_augmentedTree`
- `hairpin_has_unbracketed_treeEdge`
- `detectHairpins_sound`
- `detectHairpins_complete`

Flubbles:

- `IsFlubbleBoundary`
- `Flubble`
- `CanonicalClassBoundary`
- `flubbleBoundary_same_rootLeaf_path`
- `flubbleBoundary_has_multiple_paths`
- `flubbleBoundary_twoBlackEdgeDisconnectable`
- `flubble_boundary_is_ceep`
- `cycleClass_boundary_unique`
- `reverseDfs_flubbleBoundaries_sound`
- `reverseDfs_flubbleBoundaries_complete`
- `reverseDfs_no_duplicate_flubbles`

Flubble tree:

- `FlubbleContains`
- `StrictFlubbleContains`
- `ImmediateFlubbleParent`
- `FlubbleTree`
- `FlubbleForest`
- `MinimalFlubble`
- `flubbleContainment_laminar`
- `immediateParent_unique`
- `buildFlubbleTree_nodes_exact`
- `buildFlubbleTree_edges_are_immediateParents`
- `flubbleTree_is_tree`
- `flubbleForest_components`
- `flubbleTree_leaf_iff_minimal`

Bounds:

- `flubbles_card_le_backEdges_card_add_one`
- `backEdges_card_lt_treeEdges_card`
- `flubble_count_le_numEdges`
- `buildFlubbleForest_output_size_linear`
- `buildFlubbleForest_linear_time`
- `buildFlubbleForest_linear_space`
- `deconstruct_linear_time`
