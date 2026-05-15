# Lean4 Formalization Architecture for povu

Task: `lean4-formalization-architecture`

Sources consumed:

- `docs/lean4-proof/paper_theorem_map.md`
- `docs/lean4-proof/povu_pipeline_inventory.md`
- `docs/lean4-proof/roadmap_quality_pass.md` for existing WG task/file
  ownership context

This document defines the Lean4 proof architecture for povu's GFA-to-VCF
pipeline. It is intentionally an architecture and task-boundary document only:
it does not scaffold a Lake project, add Lean modules, or change C++/Rust
implementation files.

## Architectural Position

The formalization should build a verified Lean reference pipeline over a
precise semantic input model, then use conformance tests to compare the current
povu implementation against that reference. The first proof target is not "the
existing C++ bytes are mechanically verified"; it is:

1. accepted GFA semantics produce well-formed finite pangenome graphs;
2. verified graph algorithms compute the intended flubbles, hairpins, and
   flubble forest for those graphs;
3. verified VCF semantics produce well-formed variant records faithful to the
   verified graph structures;
4. current povu behavior is tested for conformance against the verified
   semantic boundary.

This separates proof from implementation facts that are currently external:
`liteseq` byte parsing, filesystem IO, thread scheduling, queueing, stdout/file
formatting, and current C++ control flow.

## Consumed Assumptions and Disposition

The predecessor artifacts raised several open assumptions. The architecture
settles each one as follows, or tracks it as an explicit downstream obligation.

| Assumption or gap | Architectural decision |
| --- | --- |
| Cycle semantics in biedged graphs | Use two named cycle layers. `CarrierCycle` ranges over the finite undirected carrier used by DFS and bracket proofs. `ValidBiedgedCycle` adds black/grey alternation and sequence semantics. Algorithmic flubble correctness is first proved over carrier cycles plus well-formed biedged invariants; genomic/VCF semantics cite valid biedged walks where allele spelling matters. Theorems must state which layer they use. |
| Edge identity | Every graph representation uses stable edge ids. Edges are never modeled only as endpoint pairs. Parallel grey links, black sequence edges, virtual edges, labels, and output boundary ids all depend on identity. |
| Black-edge disconnectability | Define `KBlackEdgeDisconnectable` separately from ordinary edge disconnectability. Grey edges are not removable by black-edge cuts. Boundary convention: a flubble or hairpin boundary is represented by real black edge ids adjacent to the represented subgraph; proofs must state whether the boundary edges are included in the induced interior or only in its interface. |
| Compactness | First Lean milestones make `CompactVariationGraph` part of `PovuInputInvariant`. Collapsing linear chains and proving language/output preservation is a later optional strengthening, not a blocker for flubble correctness. |
| DFS carrier | DFS runs on the undirected carrier obtained from the biedged graph by forgetting color for traversal while retaining edge ids and provenance. This carrier choice appears explicitly in DFS, backedge, bracket, hairpin, and flubble theorems. |
| Backedge non-overlap or laminarity | Do not assume the paper phrase informally. The flubble module must introduce either `LaminarBackedgeIntervals` and prove it from povu-supported graph invariants, or prove the reverse-DFS state machine without relying on global laminarity. The proof obligation remains tracked in the flubble module. |
| Bracket-set equality | `treeEdge_cycleEquivalent_iff_bracketSet_eq` is a central theorem and must be completed before flubble boundary soundness/completeness depends on cycle classes. |
| Virtual backedges and dummy artifacts | Model edge provenance explicitly as `real`, `dummy`, or `virtual kind`. Algorithm outputs filter to real black boundaries. Virtual edges may be used in the reference algorithm but cannot appear as reported flubble, hairpin, or VCF boundary edges. |
| Hairpin detector completeness | The paper's observation is not enough. The hairpin module must prove soundness and completeness for the Lean reference detector, including the boundary-pair convention and near-miss cases. |
| Flubble boundary canonicalization | Cycle-equivalence classes may contain more than two eligible edges. The flubble module must define a canonical open/close rule over DFS order and prove uniqueness/no duplicates. |
| Root/dummy flubble-tree node | The Lean flubble tree includes an explicit dummy root in the tree data structure, matching PVST. Theorems distinguish real flubble nodes from the dummy root. Leaves/minimal-flubble theorems quantify over real nodes. |
| Degenerate components | Semantic behavior is specified instead of inherited silently from current C++. Empty components, components with fewer than three graph vertices, isolated vertices, no-flubble components, self-loops, and all-tip components must have named Lean cases or be rejected by `PovuInputInvariant`. |
| Complexity proof target | First prove output-size and structural bounds in Lean. Full linear-time proofs require an executable cost model and are a later phase unless a downstream task explicitly implements that model. External benchmarks and current C++ runtime claims remain validation evidence, not Lean theorems. |
| Paper typo in Theorem 1 | State the theorem unambiguously as `flubble_count_le_numEdges`. |
| GFA byte parsing by `liteseq` | Byte-level parsing remains outside the trusted Lean proof path initially. Lean owns a semantic GFA record model and optional pure parser spec; conformance tests compare `liteseq`/povu behavior against the semantic boundary. |
| Topology-only versus VCF-capable graph construction | Model two construction contracts: topology graph construction for `decompose`, and VCF-capable construction with labels/references for `call` and `gfa2vcf`. |
| PVST round-trip loses rich metadata | Treat `.pvst` as an external serialization format. The verified semantic PVST/flubble forest includes route-visible data needed by VCF. Current `.pvst` read/write behavior is tested by conformance, not trusted as proof. |
| Pressure valves and output order in current VCF generation | Limits such as `MAX_FLUBBLE_STEPS`, chunking, queueing, and map iteration order are implementation behavior. The Lean reference specifies deterministic semantic order; conformance either normalizes order where semantics permit or flags mismatches. |

## Trusted Computing Base

The trusted computing base (TCB) is concrete and deliberately narrow for the
Lean proof path.

Trusted for Lean proofs:

- the Lean kernel, elaborator, compiler, Lake build, and selected standard
  library/mathlib facts imported by the project;
- the definitions chosen for semantic GFA records, finite graph structures,
  valid walks, flubbles, hairpins, flubble trees, and VCF records;
- any theorem explicitly imported from external Lean libraries;
- small theorem-free wrappers used only to run executable Lean examples.

Outside Lean proof and validated by tests or conformance:

- `liteseq::gfa_new` and all byte-level GFA parsing behavior;
- C++ `mto::from_gfa`, `.pvst` read/write, VCF text writer, CLI parsing,
  filesystem IO, temp-directory creation/removal, threads, queues, and stdout
  or split-file routing;
- conversion scripts or harness code that normalize current povu output into
  the semantic comparison format;
- string formatting details for dates, floating-point `AF`, duplicated header
  lines, and textual VCF record rendering unless a later VCF task proves a pure
  Lean formatter;
- performance measurements of the production C++ implementation;
- any generated code or future extracted artifact until covered by its own
  conformance task.

Proven in Lean:

- well-formed semantic GFA construction implies the core graph invariants used
  by algorithms;
- graph decomposition, DFS/bracket infrastructure, cycle-equivalence
  characterization, flubble detection, hairpin detection, and flubble-tree
  construction are sound and complete for the supported invariant set;
- VCF semantic records derived from verified structures are well formed and
  faithful to reference paths, allele slices, genotypes, coordinates, and
  supported INFO fields;
- output-size bounds, including at most `|E|` flubbles under the final
  supported assumptions.

Tested externally:

- current povu accepts/rejects actual GFA files according to the semantic GFA
  subset;
- current povu's `.pvst` and VCF outputs match Lean semantic outputs after
  documented normalization;
- current C++ runtime, memory use, chunking, thread behavior, and pressure
  valves meet project expectations.

## Lean Module Layout

All modules live under namespace `PovuLean`. The scaffold task should create
only the build skeleton and import aggregators; substantive definitions belong
to the downstream owners listed here.

### Top-Level Aggregators

- `PovuLean.lean`
  - Root import aggregator.
- `PovuLean/Core.lean`
  - Imports stable core modules.
- `PovuLean/GFA.lean`
  - Imports GFA semantic modules.
- `PovuLean/Algorithms.lean`
  - Imports completed algorithm families.
- `PovuLean/VCF.lean`
  - Imports VCF semantic modules.
- `PovuLean/Pipeline.lean`
  - Added by synthesis, after all proof modules exist.

### Core: `PovuLean/Core/**`

Owned by `lean4-core-graph-model`.

Recommended modules:

- `Core/Ids.lean`
  - `SegmentId`, `Side`, `OrientedSegment`, `SideVertex`, `EdgeId`,
    `PathId`, `SampleId`, `HaplotypeId`, `ComponentId`.
- `Core/Finite.lean`
  - finite id sets, lookup invariants, no-duplicate facts.
- `Core/Graph.lean`
  - finite carrier graph with edge ids and endpoints.
- `Core/Biedged.lean`
  - black/grey edge colors, labels on black edges, side vertices, provenance,
    incident black uniqueness.
- `Core/Walk.lean`
  - walks, carrier cycles, simple cycles, valid black/grey alternating walks,
    sequence spelling hooks.
- `Core/Paths.lean`
  - reference paths, path steps, reverse traversal, PanSN/raw metadata hooks.
- `Core/Components.lean`
  - connected components, skipped/degenerate component policy.
- `Core/Cuts.lean`
  - edge cuts, black-edge cuts, disconnectability.
- `Core/DFS.lean`
  - generic DFS tree, ancestor order, tree paths, non-tree/backedge
    comparability on the carrier.
- `Core/CycleEquiv.lean`
  - edge and node cycle equivalence, CEEP/CEVP naming, equivalence proofs.
- `Core/Positions.lean`
  - zero-based half-open internal intervals and one-based VCF coordinate
    conversion facts.
- `Core/Complexity.lean`
  - size measures such as vertex count, real edge count, black edge count,
    output-size bounds; no runtime cost model unless later introduced.

### GFA: `PovuLean/GFA/**`

Owned by `lean4-gfa-spec`.

Recommended modules:

- `GFA/Syntax.lean`
  - semantic records for supported GFA 1.0 segment, link, and path/walk data.
- `GFA/WellFormed.lean`
  - unique numeric segment ids, link endpoint validity, supported orientation
    and overlap policy, reference/path consistency.
- `GFA/Semantics.lean`
  - mapping from semantic GFA records to `Core.Biedged` graphs.
- `GFA/ParserSpec.lean`
  - optional pure parser contract if downstream chooses to specify strings.
    If byte parsing stays external, this module states accepted semantic input
    assumptions rather than proving a byte parser.
- `GFA/ToGraph.lean`
  - construction theorems:
    `acceptedGfa_toGraph_wellFormed`,
    `topologyGraph_satisfies_decomposeInvariant`, and
    `vcfGraph_satisfies_callInvariant`.

### Flubble Algorithms: `PovuLean/Algorithms/Flubble/**`

Owned by `lean4-flubble-correctness`.

Recommended modules:

- `Algorithms/Flubble/InputInvariant.lean`
  - `PovuInputInvariant`, compactness, real/dummy/virtual provenance
    requirements, supported degenerate cases.
- `Algorithms/Flubble/AugmentedTree.lean`
  - dummy roots, tip augmentation, virtual backedges, real-edge filtering,
    `T'_s` representation.
- `Algorithms/Flubble/Brackets.lean`
  - bracket lists/sets, fundamental cycles, class state.
- `Algorithms/Flubble/Spec.lean`
  - `IsFlubbleBoundary`, `Flubble`, canonical open/close rule.
- `Algorithms/Flubble/Detect.lean`
  - pure executable reference detector over finite structures.
- `Algorithms/Flubble/Correctness.lean`
  - soundness, completeness, no duplicates, canonical output.
- `Algorithms/Flubble/Bounds.lean`
  - count and output-size bounds.

Shared intermediate definitions needed by hairpin or flubble-tree tasks should
be exported from this module family rather than copied into those tasks.

### Hairpin Algorithms: `PovuLean/Algorithms/Hairpin/**`

Owned by `lean4-hairpin-correctness`.

Recommended modules:

- `Algorithms/Hairpin/Spec.lean`
  - `IsHairpin`, loop/stem decomposition, supported inversion interpretation.
- `Algorithms/Hairpin/Detect.lean`
  - pure executable reference detector and boundary extraction.
- `Algorithms/Hairpin/Correctness.lean`
  - soundness, completeness, relationship to unbracketed tree edges, virtual
    edge separation from flubbles.
- `Algorithms/Hairpin/Examples.lean`
  - positive and near-miss hairpin examples.

### Flubble Tree: `PovuLean/Algorithms/FlubbleTree/**`

Owned by `lean4-flubble-tree-correctness`.

Recommended modules:

- `Algorithms/FlubbleTree/Spec.lean`
  - real flubble nodes, dummy root, containment, route-visible data.
- `Algorithms/FlubbleTree/Build.lean`
  - pure executable tree/forest builder from verified flubbles.
- `Algorithms/FlubbleTree/Correctness.lean`
  - laminarity, immediate parent uniqueness, tree/forest invariants,
    nodes-exact theorem, leaves/minimal theorem.
- `Algorithms/FlubbleTree/Ordering.lean`
  - deterministic child/order convention used by VCF and conformance.

### VCF: `PovuLean/VCF/**`

Owned by `lean4-vcf-semantics`.

Recommended modules:

- `VCF/Model.lean`
  - semantic VCF header, contig, sample, genotype, INFO, and record types.
- `VCF/Allele.lean`
  - haplotype slices, traversal strings, DNA extraction, REF/ALT construction.
- `VCF/Variant.lean`
  - supported variant kinds: `INS`, `DEL`, `SUB`, and `SUBR` when the upstream
    hairpin/SNE boundary is available.
- `VCF/Genotype.lean`
  - sample/ploidy metadata, PanSN grouping, raw-reference fallback semantics,
    missing allele handling.
- `VCF/Semantics.lean`
  - translation from verified flubble/hairpin/tree structures and selected
    reference paths to semantic records.
- `VCF/FormatSpec.lean`
  - optional pure formatting spec. If not implemented, state exactly which
    text rendering facts remain in the external TCB.
- `VCF/Correctness.lean`
  - record well-formedness and faithfulness theorems.

### Synthesis and Conformance

Owned by later tasks:

- `PovuLean/Pipeline.lean`
  - synthesis-owned integrated theorem statement after the component modules
    are completed.
- `docs/lean4-proof/proof_obligation_check.md`
  - synthesis report.
- `tests/lean4_conformance/**` or the convention chosen by the harness task
  - Rust/C++ versus Lean comparison harness.
- `docs/lean4-proof/conformance.md`
  - trusted-boundary and fixture documentation.

## Naming Conventions

- Namespace root: `PovuLean`.
- Structures and inductives use UpperCamelCase:
  `BiedgedGraph`, `GFA.Record`, `FlubbleTree`, `VcfRecord`.
- Predicates use either noun phrases or `Is...`:
  `WellFormed`, `CompactVariationGraph`, `IsHairpin`,
  `IsFlubbleBoundary`.
- Executable functions use lowerCamelCase:
  `buildAugmentedTree`, `detectFlubbles`, `buildFlubbleForest`,
  `emitVcfRecords`.
- Theorems use lowerCamelCase with a suffix that states proof role:
  `detectFlubbles_sound`, `detectFlubbles_complete`,
  `buildFlubbleForest_nodes_exact`.
- Use `Carrier` in names for graph-color-forgetting facts:
  `carrierCycle_of_validBiedgedCycle`.
- Use `Real`, `Dummy`, or `Virtual` in names where provenance matters:
  `reportedFlubbleBoundary_real`.
- Use `Semantic` or `Raw` to distinguish parsed/validated values from unchecked
  byte/string forms.
- Avoid ambiguous abbreviations in public theorem names. Short aliases such as
  `CEEP` and `CEVP` may exist, but the canonical theorem names should spell out
  cycle equivalence.

## Executable and Specification Split

Each major area should keep three layers separate.

### Parser and GFA Layer

Specification layer:

- semantic GFA records;
- accepted subset and rejected/undefined cases;
- construction of topology-only and VCF-capable graph views;
- theorem that accepted records satisfy core graph invariants.

Executable layer:

- optional pure Lean parser from `String` to semantic records;
- examples that evaluate semantic construction on small records.

External layer:

- actual byte input, `liteseq`, file paths, malformed-line diagnostics, and
  current C++ parser behavior.

If the pure Lean parser is not built in the first GFA task, the proof boundary
starts at semantic GFA records. The conformance harness then tests that povu's
external parser maps real files to equivalent semantic records or reports
documented unsupported behavior.

### Graph Algorithm Layer

Specification layer:

- mathematical definitions of carrier cycles, valid biedged walks,
  disconnectability, hairpins, flubbles, containment, and VCF-relevant route
  semantics.

Executable layer:

- pure Lean finite algorithms for component decomposition, tip augmentation,
  DFS, bracket/class computation, flubble detection, hairpin detection, and
  flubble tree construction.

Correctness layer:

- soundness/completeness theorems tying executable outputs to specifications;
- no-duplicate/canonical ordering theorems;
- real-edge filtering theorems for dummy and virtual artifacts.

The executable Lean algorithms do not need to mimic C++ data structures such as
thread pools, bracket-list mutability, or PVST serialization internals. They do
need enough deterministic structure for conformance diagnostics to map current
povu outputs back to semantic flubbles, hairpins, and records.

### VCF Layer

Specification layer:

- semantic VCF records independent of text serialization;
- coordinate conventions, allele sequences, traversal strings, genotypes, INFO
  fields, and selected reference/sample behavior.

Executable layer:

- pure Lean translation from verified graph variants and selected references to
  semantic VCF records;
- optional text renderer if downstream chooses to verify formatting.

External layer:

- current C++ VCF headers, `fileDate`, stream routing, duplicate `GT` header
  line, floating-point `AF` formatting, record ordering from map/chunk
  iteration, and file IO.

Conformance should compare semantic records first. Byte-for-byte VCF comparison
is a later, narrower target after the formatting spec is stabilized.

## Data Representation Strategy

### Graphs and Edges

Use a finite graph with explicit id sets:

- vertices are `SideVertex` values, representing a segment id and side;
- black edges connect the two sides of a segment and carry the segment label;
- grey edges represent GFA links between sides;
- every edge has an `EdgeId`, color, endpoint pair, and provenance;
- every side vertex has at most one incident real black edge;
- graph labels live on black edge ids, not on carrier vertices.

This matches the paper's biedged view while preserving the identity needed by
povu boundaries, parallel links, labels, and VCF traversal strings.

### Paths, Samples, and References

Represent reference paths independently from the graph:

- `ReferencePath` has a stable `PathId`, raw tag, optional PanSN sample,
  optional haplotype id, and ordered steps;
- each `PathStep` stores segment id, orientation, step index, and derived
  sequence interval/locus data when labels are available;
- `SampleMetadata` groups references into VCF genotype columns, including
  PanSN phased columns and raw-reference fallback behavior;
- the vertex-to-reference step matrix used by current povu is represented
  semantically as a relation or finite map with sorted-step invariants.

Theorems should distinguish "all graph paths" from "selected references for
VCF output", because the inventory shows that current headers and genotype
columns may include more graph metadata than the selected reference ids.

### Positions and Coordinates

Use zero-based half-open coordinates internally for reference slices and region
filters. Define a conversion to one-based VCF `POS`.

The VCF module should encode current povu's documented per-variant POS rules:

- `DEL` and `INS`: left-anchor locus minus one after moving to the next step;
- `SUB` and `SUBR`: locus of the first changed step;
- region filters are half-open internally and, in current povu, use route
  endpoint loci rather than all internal vertices.

The Lean semantic spec may choose a stricter or cleaner region-overlap model,
but any difference from current povu must be captured in conformance
documentation.

### Alleles and Traversals

Represent alleles as semantic values, not just strings:

- `Traversal` is an ordered list of oriented segment ids or black edge ids;
- `HaplotypeSlice` identifies a reference path, start step, and length;
- `AlleleSequence` is derived from graph labels plus traversal orientation;
- REF and ALT construction records which context endpoints were dropped or
  retained;
- equality/grouping of alternate alleles is by semantic sequence plus
  traversal policy chosen by the VCF module.

Reverse-complement semantics are required before claiming that hairpins are
biological inversions. Until then, hairpin correctness is graph-structural and
`SUBR` VCF semantics are a supported record convention rather than a proof of
biological inversion interpretation.

### PVST and Flubble Trees

Use one semantic tree type for the proof:

- one dummy root per component tree;
- real nodes for generic flubbles and supported subfamilies when modeled;
- parent/child relation over semantic containment;
- route-visible start/end ids and direction for VCF;
- explicit height/depth facts for `LV`.

Current `.pvst` text round-tripping is conformance scope. The Lean tree should
not depend on fields that the inventory says are lost after `.pvst` reading
unless those fields are reconstructed or proven unnecessary.

### VCF Records

Use a semantic `VcfRecord` with:

- `chrom`, `pos`, `id`, `ref`, `alts`, `qual`, `filter`;
- genotype column map;
- INFO fields `AC`, `AF`, `AN`, `NS`, `AT`, `VARTYPE`, `TANGLED`, and
  optional `ES` and `LV`;
- proof that REF and ALT strings are nonempty in the supported subset;
- proof that genotype allele numbers reference existing REF/ALT entries;
- deterministic ordering key for records, even if conformance must normalize
  current C++ order initially.

## Theorem Roadmap

The proof roadmap follows the dataflow from semantic input to semantic VCF.

### 1. Core Graph Invariants

Main obligations:

- finite id lookup soundness and no-duplicate ids;
- black-edge uniqueness at side vertices;
- valid alternating walk definitions;
- carrier cycle and valid biedged cycle relationship;
- connected components partition graph vertices/edges;
- black-edge cuts and disconnectability definitions.

Representative theorem names:

- `edgeCycleEquivalent_equivalence`
- `nodeCycleEquivalent_equivalence`
- `validBiedgedCycle_to_carrierCycle`
- `components_partition_vertices`
- `kBlackEdgeDisconnectable_cutSound`

### 2. GFA Semantics to Graph Construction

Main obligations:

- accepted GFA semantic records have unique numeric segment ids and valid link
  endpoints;
- segment records create real labelled black edges;
- link records create real grey edges;
- path records create valid reference steps;
- topology-only and VCF-capable construction satisfy their respective
  invariants.

Representative theorem names:

- `acceptedGfa_toBiedged_wellFormed`
- `gfaSegments_are_realBlackEdges`
- `gfaLinks_are_realGreyEdges`
- `gfaPaths_are_validReferencePaths`
- `vcfGraph_hasLabelsAndReferences`

### 3. DFS, Augmentation, Brackets, and Cycle Classes

Main obligations:

- tip/dummy augmentation connects components as specified;
- real cycles are preserved by augmentation;
- DFS non-tree edges are ancestor comparable in the carrier;
- fundamental cycles correspond to backedges;
- bracket sets characterize tree-edge cycle equivalence.

Representative theorem names:

- `dummyRoot_connected`
- `dummyRoot_preserves_realCycles`
- `dfs_nonTreeEdge_comparable`
- `backEdge_creates_fundamentalCycle`
- `treeEdge_in_fundamentalCycle_iff_bracket`
- `treeEdge_cycleEquivalent_iff_bracketSet_eq`

### 4. Flubble Detection Correctness

Main obligations:

- canonical open/close rule for cycle-equivalence classes;
- emitted boundaries are real black tree edges;
- soundness: every emitted pair is a flubble boundary;
- completeness: every supported flubble boundary is emitted once;
- no duplicate flubbles;
- count/output-size bounds.

Representative theorem names:

- `cycleClass_boundary_unique`
- `detectFlubbles_reportedBoundaries_real`
- `detectFlubbles_sound`
- `detectFlubbles_complete`
- `detectFlubbles_noDuplicates`
- `flubble_count_le_numEdges`

### 5. Hairpin Detection Correctness

Main obligations:

- formal loop/stem hairpin specification;
- paper observation as a supporting theorem;
- executable detector soundness/completeness;
- virtual backedges separate tips/hairpins from flubbles;
- boundary pairs match the chosen convention.

Representative theorem names:

- `hairpin_lifts_to_augmentedTree`
- `hairpin_has_unbracketed_treeEdge`
- `detectHairpins_sound`
- `detectHairpins_complete`
- `hairpin_virtualBackedges_separate_from_flubbles`

### 6. Flubble Tree Correctness

Main obligations:

- flubble containment relation;
- laminarity or equivalent no-crossing theorem;
- unique immediate parent;
- produced forest nodes exactly match real flubbles plus dummy roots;
- parent edges encode immediate containment;
- leaves correspond to minimal flubbles;
- deterministic child ordering.

Representative theorem names:

- `flubbleContainment_laminar`
- `immediateParent_unique`
- `buildFlubbleTree_nodes_exact`
- `buildFlubbleTree_edges_are_immediateParents`
- `buildFlubbleForest_components`
- `flubbleTree_leaf_iff_minimal`

### 7. VCF Semantic Correctness

Main obligations:

- selected references and genotype columns are well formed;
- RoV/route endpoints exist on selected references;
- allele slices and traversal strings match graph labels;
- REF/ALT construction follows the supported variant type rules;
- INFO and genotype fields are internally consistent;
- semantic records are faithful to verified flubble/hairpin/tree structures.

Representative theorem names:

- `selectedReferences_wellFormed`
- `hapSlice_sequence_matchesTraversal`
- `refAlt_nonempty`
- `genotypes_referenceExistingAlleles`
- `infoFields_consistent`
- `emitVcfRecords_sound`

### 8. Pipeline Synthesis

Main obligations:

- integrated module imports have no duplicated definitions;
- assumptions from GFA/core/algorithm/VCF modules align;
- final theorem states the verified semantic pipeline boundary;
- residual trusted/external parts are documented for conformance.

Representative theorem name:

- `semanticGfaToVcf_correct`

## Complexity Claims

Lean should prove structural bounds first:

- connected-component output partitions are bounded by input vertices/edges;
- tip/dummy/virtual augmentation adds bounded edges under stated rules;
- detected flubble count is at most `|E|` for the final graph size measure;
- flubble forest node/edge counts are linear in the number of real flubbles and
  therefore in input edges under the count theorem;
- semantic VCF record counts are bounded by the chosen RoV/reference/allele
  model where such a bound is meaningful.

Lean should not initially claim wall-clock linear time for current povu. A full
linear-time proof requires:

- executable algorithms whose operations are modeled at array/list/map cost;
- proofs that each traversal scans each vertex/edge a constant or amortized
  constant number of times;
- explicit assumptions for finite map/set implementations;
- a separate story for current C++ threads, queues, `liteseq`, and pressure
  valves.

Therefore:

- prove `flubble_count_le_numEdges` and flubble-forest output-size bounds in
  the first algorithm wave;
- track executable cost-model theorems such as
  `buildFlubbleForest_linear_time` as later optional obligations;
- validate current production performance externally with benchmarks and
  regression tests, not as Lean proof evidence.

## Extraction and Conformance Strategy

The near-term strategy is "Lean as executable oracle", not production
extraction.

1. Lean modules define pure semantic reference functions for graph
   construction from semantic GFA, flubble/hairpin/tree detection, and VCF
   semantic record emission.
2. The conformance harness runs current povu on small GFA fixtures and captures
   graph/PVST/VCF outputs.
3. Harness adapters normalize current outputs into semantic comparison forms:
   records, route endpoints, boundary ids, allele traversals, genotype maps,
   and INFO fields.
4. Comparisons prefer semantic equality over byte-for-byte text equality.
   Formatting-sensitive checks can be added after `VCF/FormatSpec.lean`
   stabilizes.
5. Mismatches are classified as:
   - implementation bug;
   - Lean spec bug;
   - unsupported input;
   - formatting-only difference;
   - intentionally external behavior.
6. Future language targets should either call extracted/ported verified
   algorithms or pass the same conformance suite. They should not be treated as
   independent rewrites.

The harness should include enough diagnostics to map a mismatch back to the
semantic object: GFA fixture, component id, boundary ids, reference path, VCF
record key, and differing semantic fields.

## Incremental Lake Milestones

Every implementation milestone should end with `lake build` from a clean
checkout. Trusted modules in the proof path should contain no `sorry` or
`admit`; draft modules with placeholders must be outside the trusted import
path and documented by the owning task.

Recommended sequence:

1. `lean4-scaffold`
   - Create Lake project, root namespace, import aggregators, and documented
     build command.
   - Validation: `lake build`.
2. `lean4-core-graph-model`
   - Add core ids, finite graph, biedged graph, paths, walks, components, cuts,
     DFS hooks, and cycle equivalence basics.
   - Validation: `lake build` plus small examples for path, cycle,
     branch/rejoin, and inversion-like traversal.
3. `lean4-gfa-spec`
   - Add semantic GFA records, accepted subset, and graph-construction
     contracts.
   - Validation: `lake build` plus accepted and rejected semantic examples.
4. `lean4-flubble-correctness`
   - Add input invariant, augmented tree, brackets, flubble spec, reference
     detector, correctness, and count bounds.
   - Validation: `lake build` plus simple, nested, overlapping/negative
     examples.
5. `lean4-hairpin-correctness`
   - Add hairpin spec, detector, and correctness using flubble intermediate
     interfaces.
   - Validation: `lake build` plus positive and near-miss hairpin examples.
6. `lean4-flubble-tree-correctness`
   - Add containment, tree builder, tree/forest correctness, route-visible
     data, and deterministic ordering.
   - Validation: `lake build` plus nested, sibling, single, and no-flubble
     examples.
7. `lean4-vcf-semantics`
   - Add semantic VCF model, allele/genotype/INFO semantics, and record
     correctness.
   - Validation: `lake build` plus SNP-like, insertion/deletion-like,
     complex/nested, hairpin/SUBR, no-variant, and ordering examples where
     supported.
8. `lean4-proof-synthesis-check`
   - Integrate modules into `PovuLean/Pipeline.lean`, reconcile assumptions,
     and publish `proof_obligation_check.md`.
   - Validation: full `lake build`; no placeholders in trusted import path.
9. `lean4-conformance-harness`
   - Build the current-povu versus Lean semantic comparison harness.
   - Validation: documented harness command, existing project tests, and
     `lake build`.
10. `lean4-e2e-validation-corpus`
    - Expand deterministic fixtures across parser, graph, flubble, hairpin,
      flubble-tree, and VCF behavior.
    - Validation: full conformance command plus fixture notes.

## WG Graph Refinement

The existing WG roadmap already has non-overlapping file ownership for the
main proof modules. One dependency refinement follows from the pipeline
inventory: the end-to-end corpus should wait for `fix-sne-subr`, because that
task addresses a known current-povu SUBR/SNE generation issue and the corpus is
expected to include hairpin inversion/SUBR cases. The architecture task should
record that dependency in WG metadata.

No additional proof-module fanout is needed here. The current graph already
separates core graph, GFA, flubble, hairpin, flubble-tree, VCF, synthesis,
conformance, corpus, and future language-targeting work with non-overlapping
file scopes.

## Acceptance Checklist Mapping

- Consumes predecessor artifacts:
  - The assumption table resolves or tracks the paper theorem-map gaps and the
    pipeline inventory's external parser, PVST, VCF, edge-case, and testing
    gaps.
- Independent module boundaries:
  - Core graph, GFA, flubble, hairpin, flubble-tree, VCF, synthesis, and
    conformance ownership are separated by module path and task.
- Concrete trusted base:
  - Parser/string/IO, C++ production behavior, VCF byte formatting, harness
    normalization, and performance claims are explicitly outside Lean proof
    unless later verified.
- Incremental `lake build` validation:
  - The milestone sequence gives a build gate after each proof stage.
