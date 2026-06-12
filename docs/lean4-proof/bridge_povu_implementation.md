# Bridge From Current Povu To Lean Flubble Contracts

Date: 2026-06-12
Task: `bridge-povu-implementation`

This report maps the current C++ and Rust implementation surfaces to the Lean
linear-proof contracts. It is a bridge audit and conformance note, not a formal
runtime proof of current povu. The only formal claims are the Lean theorems in
`PovuLean`; C++ and Rust tests below are implementation evidence.

## Contract Boundary

`PovuLean/Complexity/Decomposition.lean` states the full semantic decomposition
composition theorem as a conditional theorem. Lines 12-30 name the external
stage obligations, and lines 32-36 state that implementation obligations are
limited to supplying these contracts for a concrete execution and relating that
execution to the semantic Lean witnesses. The theorem body composes costs only
after receiving `SemanticDecomposeStageContracts`,
`SemanticDecomposeIntermediateSizeBounds`, supported hierarchy input, and VCF
semantic-call witnesses.

The proof-owned semantic pieces are:

- `PovuLean/Algorithms/Flubble/InputInvariant.lean`: `candidateStack` is the
  filtered list of real black tree links, and `CycleClassAssignment.Correct`
  requires soundness and completeness against `CycleEquivalentEdges`.
- `PovuLean/Algorithms/Flubble/Spec.lean`: a base flubble boundary is a pair of
  real black tree edges, canonicalized by edge id, where the close edge is the
  first later same-class edge after a nonempty gap.
- `PovuLean/Algorithms/Flubble/Detect.lean`: `detectFlubbles` traverses the
  candidate stack, emits the canonical close-after-gap boundary, deduplicates
  exact pairs, and proves membership, duplicate-freedom, and size bounds.
- `PovuLean/Algorithms/FlubbleTree/Spec.lean`: hierarchy construction is defined
  over boundary spans in the candidate stack. Supported inputs require every
  boundary to have a span and all spans to be laminar.
- `PovuLean/Complexity/Flubble.lean` and
  `PovuLean/Complexity/Decomposition.lean`: extraction and hierarchy costs are
  proof-side cost records, not observed runtime counters from C++ or Rust.

## Implementation Stage Map

| Stage | Current implementation surface | Lean contract relation | Status |
| --- | --- | --- | --- |
| Component decomposition | C++ conformance CLI path configures and builds `povu`; componentization is visible in tests such as `tests/integration_tests/pvst_tests.cc`, where `bd::VG::componetize` feeds one component into `find_flubbles`. The Lean/Rust conformance harness exercises `povu gfa2vcf` over fixtures in `tests/lean4_conformance/src/main.rs`. | Must supply `componentDecompositionLinear` and a semantic relation from byte-level GFA components to `doc.toGraph`. | Conformance-tested at final VCF/structure-export level only. Not formally bridged. |
| Tip/dummy augmentation | C++ spanning-tree construction happens before `find_flubbles`; `src/povu/algorithms/flubbles.cpp` consumes `pst::Tree` and adds simplifying/capping backedges in `handle_vertex`. | Must supply `tipDummyAugmentationLinear`, `augmentedGraphLinear`, and preservation of supported semantic graph/frame inputs. | Open. The exact dummy/tip witness relation to Lean `TraversalFrame` is not exported. |
| Traversal-frame construction | C++ `compute_eq_class_stack` iterates spanning-tree vertices and collects black incoming tree edges into the equivalence-class stack. Lean `candidateStack` filters `frame.treeLinks` to real black edges. | Must show the C++ stack order equals the Lean `candidateStack frame` order for the same accepted graph/frame. | Conformance-tested indirectly by PVST and structure fixtures. Missing a machine-checkable frame/stack export. |
| Bracket and cycle-class assignment | C++ `handle_vertex` computes bracket lists, assigns new classes to bracket/backedge state, and labels incoming tree edges. `compute_eq_class_metadata` then records the next later stack index for each class. | Must supply `CycleClassAssignment.Correct`: equal classes imply cycle-equivalent boundary candidates and cycle-equivalent candidates receive equal classes. | Open for formal proof. C++ tests observe expected PVST labels/hierarchy, but do not prove soundness/completeness. |
| Flubble extraction | C++ `add_flubbles` emits a PVST vertex when `i + 1 < next_seen[i]`, matching the nonempty-gap condition against the next later same class. Rust `povu-rs/src/native_gfa.rs` exposes `detect_flubble_stack`, documented as a Rust port of Lean `detectStack`, with tests in `povu-rs/tests/native_gfa_vcf_tests.rs`. | Must match Lean `closeAfterGap?`, `firstSameClass?`, `Boundary.ordered`, and duplicate-free output. | Rust semantic detector is conformance-tested directly. C++ is conformance-tested indirectly through PVST/structure fixtures. |
| Hierarchy construction | C++ `add_flubbles` maintains a class stack and parent vertex pointer, adding PVST edges as it emits flubbles. Existing C++ `PVSTTest.VertexHierarchy` checks one nested hierarchy. The Lean/Rust conformance harness now derives stack spans from the exported C++ flubble debug stack and compares exported PVST parent edges against canonical nearest strict span-containment parents. Lean hierarchy uses span containment and rejects unsupported non-laminar boundaries. | Must show emitted parent edges equal canonical strict span-containment parents under supported laminar inputs. | Conformance-tested for selected exported fixtures, including nested and sibling hierarchy cases. Still not a proof for arbitrary supported inputs. |
| VCF/structure export | C++ conformance harness compares `povu gfa2vcf` VCF and `--structure-export` JSON to Lean references. Rust native VCF tests compare the Rust semantic extractor to Lean fixture expectations. | Lean `Pipeline.semanticGfaToVcf_correct` proves semantic emission when supplied semantic calls and hierarchy witnesses. | Conformance-tested. Not a proof that runtime C++ calls are the semantic witnesses. |

## New Conformance Evidence

This pass expands the Rust native detector tests in
`povu-rs/tests/native_gfa_vcf_tests.rs`:

- `flubble_stack_port_preserves_lean_first_later_same_class_rule` checks that a
  candidate closes at the first later same-class candidate, even when a farther
  same-class candidate can later open its own boundary.
- `flubble_stack_port_canonicalizes_and_deduplicates_like_lean_detector` checks
  ID-order canonicalization plus exact duplicate removal, matching the Lean
  `Boundary.ordered` and `uniqueBoundaries` behavior.

These tests are deliberately narrow. They support the Rust semantic port of the
Lean detector, but they do not prove C++ runtime conformance or linear runtime.

The existing Lean/Rust conformance harness already compares:

- final normalized VCF records from current `povu gfa2vcf` against Lean fixture
  references;
- canonical structure-export JSON against Lean structure references;
- expected implementation failures for unsupported/malformed GFA fixtures.

Those checks provide readable fixture-level diagnostics through
`tests/lean4_conformance/src/main.rs`. They remain conformance evidence only.

`bridge-check-cycle` adds a direct C++ class-id conformance check to the
Lean/Rust fixture harness. For the small supported structure fixtures, the
harness now carries an independent cycle-equivalence oracle over exported
candidate `tree_edge_id`s and checks both directions of the Lean
`CycleClassAssignment.Correct` obligation:

- soundness: any two exported candidate edges with the same C++ `class_id` must
  be grouped by the oracle;
- completeness: any two candidate edges grouped by the oracle must have the
  same exported C++ `class_id`.

Failures name the fixture, frame index, candidate edge ids, class ids, stack
orders, and the failed direction. Unit coverage includes the positive
`insertion-flubble` class assignment plus synthetic soundness and completeness
failures to keep diagnostics readable. This is still a finite fixture oracle,
not a general proof of C++ bracket-stack cycle equivalence for arbitrary inputs.

`bridge-check-hierarchy` adds a direct PVST parent-edge checker to the same
fixture harness. For each supported C++ structure export, the harness:

- reads the exported flubble debug stack entries, including stack order,
  boundary vertex id, orientation, tree-edge id, and class id;
- reconstructs the emitted base-flubble boundaries using the same first later
  same-class, nonempty-gap rule as the Lean `detectStack` semantics;
- canonicalizes boundary ids to match the PVST route labels;
- derives each boundary span as the interval between the opening and closing
  stack entries;
- computes the expected PVST parent as the nearest strict containing span, or
  the component dummy root for top-level boundaries;
- compares that expected parent with the actual exported `pvst_nodes[].parent`.

Hierarchy-check failures name the fixture, boundary id, node id, stack span,
expected parent, and actual parent. Unit coverage includes nested parenthood,
sibling top-level parenthood, reverse/reverse boundary canonicalization, and a
synthetic parent mismatch diagnostic. This remains finite conformance evidence
over exported fixture stacks, not a general theorem that C++ hierarchy
construction is correct for every supported laminar input.

## Formal, Conditional, Tested, And Open Items

Formally proved in Lean:

- semantic flubble detector membership, duplicate-freedom, canonical boundary
  properties, and output-size bounds;
- hierarchy specification properties for supported laminar boundary inputs;
- semantic VCF emission correctness under supplied semantic witnesses;
- arithmetic composition from named stage and size contracts to an input-linear
  semantic decomposition bound.

Conditionally proved in Lean:

- full semantic decomposition linearity, assuming the stage contracts and
  intermediate size bounds supplied to
  `conditional_semantic_decompose_linear`;
- flubble extraction/hierarchy linearity, assuming the explicit indexed-detector
  and hierarchy cost contracts.

Conformance-tested:

- C++ CLI VCF and structure-export behavior on the Lean conformance fixtures;
- exported C++ flubble candidate class ids against fixture-sized independent
  cycle-equivalence groups for selected small supported structure fixtures;
- exported C++ PVST parent edges against canonical stack-span containment for
  supported fixture stacks, including nested and sibling cases;
- C++ PVST nested hierarchy behavior on the integration graph in
  `tests/integration_tests/pvst_tests.cc`;
- Rust native GFA-to-VCF behavior on the Lean fixture corpus;
- Rust semantic `detect_flubble_stack` behavior for close-after-gap,
  first-later-same-class, canonical orientation, and duplicate removal.

Open implementation-proof bridge items:

- byte-level GFA parsing and componentization are not related to Lean
  `GFA.Document.Accepted` by a checked translator proof;
- C++ spanning-tree, tip/dummy augmentation, and bracket-stack state are not
  exported as a Lean `TraversalFrame` witness;
- C++ class ids are not proved against `CycleClassAssignment.Correct` for
  arbitrary supported inputs; the current checker covers only enumerated small
  fixture oracle groups;
- C++ PVST parent edges are checked against Lean stack-span parent semantics on
  exported fixtures, but not proved for arbitrary supported laminar inputs;
- no runtime cost instrumentation or machine-checked refinement connects
  C++/Rust loops, maps, stacks, and allocations to the Lean `StageCost` records;
- external library and container-operation assumptions are unnamed in code-level
  proof artifacts.

## Named External Assumptions Needed For A Formal Runtime Bridge

A future machine-checkable bridge must either prove or explicitly assume:

- parser refinement: accepted byte-level GFA records produce the same semantic
  graph as Lean `GFA.Document.toGraph`;
- component refinement: C++ component decomposition preserves segment/link/path
  semantics and has total output size linear in the input graph;
- frame refinement: the C++ spanning tree plus augmentation exports the same
  ordered tree-link list and real/black/provenance flags as the Lean
  `TraversalFrame`;
- class refinement: C++ bracket classes are sound and complete for Lean
  `CycleEquivalentEdges` on boundary candidates;
- extraction refinement: C++ `next_seen` and `add_flubbles` emit exactly Lean
  `detectFlubbles` boundaries for the exported stack/classes;
- hierarchy refinement: C++ PVST parent edges are exactly the Lean canonical
  laminar span parents;
- cost model assumptions: `std::vector`, `std::unordered_map`,
  `std::unordered_set`, `std::map`, and `std::list` operations used in the
  relevant loops satisfy the expected amortized or worst-case bounds under named
  preconditions;
- FFI/serialization preservation: Rust/C++ wrappers and JSON/VCF serialization
  preserve the semantic calls and hierarchy metadata.

## Proposed Follow-Up WG Tasks

- `bridge-export-flubble`: add a debug/conformance export for the C++
  candidate stack, tree-edge metadata, class ids, and `next_seen` table for each
  fixture.
- `bridge-check-cycle`: build a checker that validates exported C++
  classes against an independent cycle-equivalence oracle on small fixtures.
- `bridge-check-hierarchy`: compare C++ PVST parent edges with Lean/Rust
  stack-span containment on exported fixture stacks.
- `bridge-cost-instrumentation`: introduce named stage counters for component
  decomposition, augmentation, traversal-frame construction, class assignment,
  extraction, and hierarchy construction.
- `bridge-checked-translator`: define a machine-checkable translation artifact
  from accepted GFA bytes and exported C++ decomposition state into Lean
  semantic witnesses. This should consume the export/checker/instrumentation
  tasks, but WG max depth prevented encoding all of those dependency edges.

Until those tasks land, current povu has conformance evidence against the Lean
semantic contracts, but not a formal runtime proof or a formal implementation
refinement theorem.
