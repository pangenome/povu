#set document(
  title: "What Lean Proves About Flubbles",
  author: "povu Lean proof notes",
)
#set page(paper: "us-letter", margin: (x: 0.75in, y: 0.7in))
#set text(font: "Libertinus Serif", size: 10pt)
#set par(justify: true, leading: 0.58em)
#set heading(numbering: "1.")

#show raw.where(block: true): it => block(
  fill: rgb("#f6f8fa"),
  inset: 8pt,
  radius: 2pt,
  width: 100%,
  it,
)

#align(center)[
  #text(size: 18pt, weight: "bold")[What Lean Proves About Flubbles]

  #text(size: 9pt)[Human-readable proof note, 2026-06-12]
]

#outline(title: [Contents], indent: auto)

= Short Answer

There are two natural questions.

#quote[
Are there only linearly many flubbles in a graph?
]

Yes, for the canonical flubbles used by the povu algorithm. Lean proves that,
for a supported traversal frame `F` over a graph `G` and a correct cycle-class
assignment `C`,

$ |"detectFlubbles"(F, C)| <= "edgeCount"(G). $

The formal theorem is
#link("https://github.com/pangenome/povu/blob/main/PovuLean/Algorithms/Flubble/Correctness.lean#L235")[#raw("PovuLean.Algorithms.Flubble.flubble_count_le_numEdges")].
It is a count theorem for the canonical detector output. It is not a theorem
about all pairwise cycle-equivalent edge pairs; that all-pairs interpretation
can be quadratic and is not the flubble definition used here.

#quote[
Is flubble finding linear time in graph size?
]

Yes for the indexed flubble-finding stage, under the ordinary word-RAM model
used in algorithms analysis. The indexed detector does one class-index lookup
and one class-index update per candidate, emits at most one canonical boundary
per candidate, and builds hierarchy metadata whose size is linear in the number
of emitted boundaries. Once the candidate stack is linearly bounded by the
graph-size measure, the whole indexed flubble stage is linear in graph size.

The formal stage theorem is
#link("https://github.com/pangenome/povu/blob/main/PovuLean/Complexity/Flubble.lean#L452")[#raw("PovuLean.Complexity.Flubble.indexed_flubble_stage_linear_in_graphSize")].
The candidate-stack version is
#link("https://github.com/pangenome/povu/blob/main/PovuLean/Complexity/Flubble.lean#L426")[#raw("indexed_flubble_stage_linear_in_candidateStack")].

The full decomposition theorem is broader and conditional: Lean composes the
flubble-stage result with named contracts for component decomposition,
augmentation, traversal-frame construction, cycle-class assignment, and
hairpin scanning. Those are standard linear graph-algorithm stages, but their
linearity is not all mechanized in Lean yet.

= The Flubble Definition Being Counted

The detector works over an ordered candidate stack

$ S(F) = [e_0, e_1, dots, e_(m-1)]. $

The stack contains the real black tree edges eligible to be flubble boundary
endpoints. A cycle-class assignment `C` labels candidate edges by the
cycle-equivalence class used in the flubble characterization.

For a candidate `e_i`, the canonical rule is:

$ j = min { k | i + 1 < k < m and C(e_i) = C(e_k) }. $

If such a `j` exists, the detector emits the ordered boundary

$ "Boundary.ordered"(e_i, e_j). $

The nonempty-gap condition `i + 1 < j` is part of the rule. The word
"canonical" matters: a candidate edge can start at most one emitted flubble,
namely the one closed by the first later same-class edge after the gap. This is
why the number of canonical flubbles is linear. Counting every same-class pair
would be a different relation and can be quadratic.

Lean names this semantic boundary predicate
#link("https://github.com/pangenome/povu/blob/main/PovuLean/Algorithms/Flubble/Spec.lean#L155")[#raw("IsFlubbleBoundary")].
The detector correctness theorem links that predicate to the executable list
returned by `detectFlubbles`:
#link("https://github.com/pangenome/povu/blob/main/PovuLean/Algorithms/Flubble/Correctness.lean#L210")[#raw("detectFlubbles_correct")].

= Why The Count Is Linear

The proof is simple once the canonical rule is fixed.

1. The raw detector visits each candidate as a possible opener.
2. For that opener it emits either zero boundaries or exactly one boundary.
3. Deduplication cannot increase the list length.

Therefore

$ |"detectFlubbles"(F, C)| <= |S(F)|. $

The stack `S(F)` is a filtered subsequence of the traversal frame's tree links.
For supported input, the tree-link list has no duplicates and every tree link
comes from the graph's link list. Since `edgeCount(G)` is the graph link count,

$ |S(F)| <= "edgeCount"(G). $

Combining the two inequalities gives

$ |"detectFlubbles"(F, C)| <= "edgeCount"(G). $

The underlying detector-size theorem is
#link("https://github.com/pangenome/povu/blob/main/PovuLean/Algorithms/Flubble/Detect.lean#L196")[#raw("detectFlubbles_length_le_supportedInput_graph_edgeCount")].
The paper-facing wrapper is
#link("https://github.com/pangenome/povu/blob/main/PovuLean/Algorithms/Flubble/Correctness.lean#L235")[#raw("flubble_count_le_numEdges")].

The theorem has a `CycleClassAssignment.Correct` hypothesis. That hypothesis is
not needed for the raw list-length inequality, but it is needed to read the
list as the list of semantic flubbles rather than merely same-class detector
outputs. The class-assignment contract is here:
#link("https://github.com/pangenome/povu/blob/main/PovuLean/Algorithms/Flubble/InputInvariant.lean#L93")[#raw("CycleClassAssignment.Correct")].

= Why Finding Them Is Linear

The direct specification says "for each opener, find the first later
same-class closer after a gap." A literal forward scan for every opener would
be quadratic. The algorithmic move is the indexed detector.

The indexed detector scans the candidate stack in reverse and maintains a
class-index table. At each candidate it performs:

+ one lookup: what is the next later candidate in this class?
+ one update: this candidate is now the latest seen candidate in this class.

Lean proves this exact operation count:
#link("https://github.com/pangenome/povu/blob/main/PovuLean/Algorithms/Flubble/Indexed.lean#L372")[#raw("detectStackIndexedRaw_one_lookup_and_update_per_candidate")].

Lean also proves the indexed detector returns the same boundaries as the direct
reference detector:
#link("https://github.com/pangenome/povu/blob/main/PovuLean/Algorithms/Flubble/Indexed.lean#L266")[#raw("detectFlubblesIndexed_eq_detectFlubbles")].
So the fast one-pass algorithm computes the same canonical flubble relation as
the obvious specification.

Under the usual word-RAM convention that class-index lookup and update are
constant-time operations, the indexed detector work is linear in `|S(F)|`.
Boundary emission is linear in the number of emitted boundaries, and the count
theorem above bounds that number by `|S(F)|` and by `edgeCount(G)`.

The hierarchy builder records one node per detected flubble and flat metadata
for parent/boundary references. Lean proves the relevant graph-edge-count
bounds:

+ #link("https://github.com/pangenome/povu/blob/main/PovuLean/Algorithms/FlubbleTree/Correctness.lean#L284")[#raw("buildHierarchy_nodes_length_le_supportedInput_graph_edgeCount")]
+ #link("https://github.com/pangenome/povu/blob/main/PovuLean/Algorithms/FlubbleTree/Correctness.lean#L309")[#raw("buildHierarchy_flatMetadataSlots_le_supportedInput_graph_edgeCount_twice")]

Putting these pieces together gives the flubble-stage theorem:

```lean
PovuLean.Complexity.Flubble.indexed_flubble_stage_linear_in_candidateStack
PovuLean.Complexity.Flubble.indexed_flubble_stage_linear_in_graphSize
```

The graph-size theorem takes the candidate-stack linear-size bound as a
hypothesis. In ordinary graph terms, that is the statement that the traversal
frame has not expanded the list of candidate edges superlinearly. The Lean
detector count theorem itself already proves the emitted flubble list is
bounded by `edgeCount(G)` under supported input.

= What Is Mechanized, Exactly?

Lean currently checks the following.

+ The semantic definition of canonical flubble boundaries over a traversal
  candidate stack.
+ Soundness and completeness of `detectFlubbles` for that boundary relation.
+ Duplicate freedom and canonicalization of the detector output.
+ The count theorem
  #link("https://github.com/pangenome/povu/blob/main/PovuLean/Algorithms/Flubble/Correctness.lean#L235")[#raw("flubble_count_le_numEdges")].
+ Equivalence of the indexed detector to the reference detector.
+ One lookup and one update per candidate in the indexed scan.
+ Linear output-size bounds for detected boundaries and flubble hierarchy
  metadata.
+ Linear-time flubble-stage composition under the explicit cost contracts.

The full semantic decomposition theorem is
#link("https://github.com/pangenome/povu/blob/main/PovuLean/Complexity/Decomposition.lean#L155")[#raw("conditional_semantic_decompose_linear")].
It proves the arithmetic composition: if the named upstream stages are linear,
and if the intermediate sizes stay linear in graph size, then the semantic
pipeline cost is linear in graph size. It also returns semantic VCF correctness
for the corresponding Lean witnesses.

= What Remains Outside This Proof

There are three boundaries worth keeping separate.

First, Lean has not introduced a graph-only set called something like
`flubbles G` and proved a cardinality theorem for that object. The current
count theorem is intentionally over the existing canonical detector list
`detectFlubbles frame classes`. That is the right theorem for the algorithm as
implemented in the Lean model.

Second, `CycleClassAssignment.Correct` is a genuine correctness obligation.
The count and runtime accounting assume a class assignment whose equal class
ids are sound and complete for cycle-equivalent candidate edges. Fixture-level
C++ audits compare exported classes with an independent cycle-equivalence
oracle, but that is conformance evidence, not a Lean proof for arbitrary C++
executions.

Third, the document is not a formal proof that the current C++ or Rust binaries
run in linear time. The C++ counters, Rust tests, and conformance harness are
implementation evidence connecting the real code to the Lean boundary. The
algorithm theorem itself is the standard word-RAM claim about the flubble
algorithm: canonical flubble count is linear, and the indexed flubble-finding
stage is linear once supplied with a linear-size candidate stack and constant
time class-index operations.

= Build And Render

The Lean names linked above are checked by:

```bash
lake build
```

The PDF is generated from this Typst source and should not be committed:

```bash
typst compile docs/lean4-proof/flubble_linear_proof.typ \
  /tmp/flubble_linear_proof.pdf
```
