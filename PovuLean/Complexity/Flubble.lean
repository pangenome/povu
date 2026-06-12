import PovuLean.Algorithms.Flubble.Indexed
import PovuLean.Algorithms.FlubbleTree.Correctness
import PovuLean.Complexity.Basic

/-!
Proof vocabulary for future flubble linear-time obligations.

The structures below name the size measures and external stage contracts that a
later proof may assume explicitly.  They are not axioms, and this module does
not state that the current povu C++ or Rust implementation satisfies them.  Any
runtime theorem about those implementations must pass through a later concrete
cost model or machine-checkable conformance bridge.
-/

namespace PovuLean
namespace Complexity
namespace Flubble

open Core
open Algorithms

/-- Input-side graph size measures used by flubble cost statements. -/
structure GraphSizes where
  vertices : Nat
  edges : Nat
  deriving Repr, DecidableEq

namespace GraphSizes

def ofGraph (g : Graph) : GraphSizes :=
  { vertices := g.segments.length
    edges := g.edgeCount }

def total (sizes : GraphSizes) : Nat :=
  sizes.vertices + sizes.edges

@[simp] theorem ofGraph_vertices (g : Graph) :
    (ofGraph g).vertices = g.segments.length := rfl

@[simp] theorem ofGraph_edges (g : Graph) :
    (ofGraph g).edges = g.edgeCount := rfl

end GraphSizes

/--
Flubble decomposition size measures.

`candidateStackLength`, `detectedFlubbleCount`, and `hierarchyOutputSize` are
kept separate so later proofs can state whether a stage is linear in the input
graph, in a reused output list, or in a hierarchy representation.
-/
structure DecompositionSizes where
  graph : GraphSizes
  candidateStackLength : Nat
  detectedFlubbleCount : Nat
  hierarchyOutputSize : Nat
  deriving Repr, DecidableEq

namespace DecompositionSizes

def inputSize (sizes : DecompositionSizes) : Nat :=
  sizes.graph.total

def outputSize (sizes : DecompositionSizes) : Nat :=
  sizes.detectedFlubbleCount + sizes.hierarchyOutputSize

def ofFrame {g : Graph}
    (frame : TraversalFrame g)
    (classes : Algorithms.Flubble.CycleClassAssignment g)
    (hierarchy : Algorithms.FlubbleTree.Hierarchy) :
    DecompositionSizes :=
  { graph := GraphSizes.ofGraph g
    candidateStackLength := Algorithms.Flubble.candidateStack frame |>.length
    detectedFlubbleCount :=
      (Algorithms.Flubble.detectFlubbles frame classes).length
    hierarchyOutputSize := hierarchy.nodes.length }

@[simp] theorem ofFrame_candidateStackLength {g : Graph}
    (frame : TraversalFrame g)
    (classes : Algorithms.Flubble.CycleClassAssignment g)
    (hierarchy : Algorithms.FlubbleTree.Hierarchy) :
    (ofFrame frame classes hierarchy).candidateStackLength =
      (Algorithms.Flubble.candidateStack frame).length := rfl

@[simp] theorem ofFrame_detectedFlubbleCount {g : Graph}
    (frame : TraversalFrame g)
    (classes : Algorithms.Flubble.CycleClassAssignment g)
    (hierarchy : Algorithms.FlubbleTree.Hierarchy) :
    (ofFrame frame classes hierarchy).detectedFlubbleCount =
      (Algorithms.Flubble.detectFlubbles frame classes).length := rfl

@[simp] theorem ofFrame_hierarchyOutputSize {g : Graph}
    (frame : TraversalFrame g)
    (classes : Algorithms.Flubble.CycleClassAssignment g)
    (hierarchy : Algorithms.FlubbleTree.Hierarchy) :
    (ofFrame frame classes hierarchy).hierarchyOutputSize =
      hierarchy.nodes.length := rfl

end DecompositionSizes

/-- A named stage-cost slot for later detector and hierarchy proofs. -/
structure StageCost where
  cost : Nat
  deriving Repr, DecidableEq

/-- Concrete linear stage contract over one input size parameter. -/
def LinearStageWith (stage : StageCost) (n c k : Nat) : Prop :=
  LinearBoundWith stage.cost n c k

/-- Existential linear stage contract over one input size parameter. -/
def LinearStage (stage : StageCost) (n : Nat) : Prop :=
  LinearBound stage.cost n

/--
Explicit operation-cost assumptions for a class-index table.

Later indexed detector proofs should accept this structure as a theorem
parameter when they need constant-time lookup/update obligations.  Supplying a
value of this structure is outside this scaffolding module.
-/
structure ClassIndexContract where
  lookupCost : Nat
  updateCost : Nat
  lookupConstant : ∃ k, lookupCost ≤ k
  updateConstant : ∃ k, updateCost ≤ k

/--
Per-candidate operation contract for the indexed detector.

`visitedCandidates` is the candidate-stack length consumed by the scan.
`lookupCount` and `updateCount` are stated separately so a concrete map model can
plug in its lookup/update costs while reusing the detector's exact one-lookup
and one-update-per-candidate theorem.
-/
structure IndexedDetectorOperationContract where
  visitedCandidates : Nat
  lookupCount : Nat
  updateCount : Nat
  lookupCount_eq : lookupCount = visitedCandidates
  updateCount_eq : updateCount = visitedCandidates
  classIndex : ClassIndexContract

namespace IndexedDetectorOperationContract

def operationCost (contract : IndexedDetectorOperationContract) : Nat :=
  contract.lookupCount * contract.classIndex.lookupCost
    + contract.updateCount * contract.classIndex.updateCost

theorem operationCost_linear_in_candidates
    (contract : IndexedDetectorOperationContract) :
    LinearBound contract.operationCost contract.visitedCandidates := by
  rcases contract.classIndex.lookupConstant with ⟨lookupK, hLookup⟩
  rcases contract.classIndex.updateConstant with ⟨updateK, hUpdate⟩
  refine ⟨lookupK + updateK, 0, ?_⟩
  unfold LinearBoundWith operationCost
  rw [contract.lookupCount_eq, contract.updateCount_eq]
  calc
    contract.visitedCandidates * contract.classIndex.lookupCost +
        contract.visitedCandidates * contract.classIndex.updateCost
        ≤ contract.visitedCandidates * lookupK +
            contract.visitedCandidates * updateK :=
      Nat.add_le_add
        (Nat.mul_le_mul_left contract.visitedCandidates hLookup)
        (Nat.mul_le_mul_left contract.visitedCandidates hUpdate)
    _ = lookupK * contract.visitedCandidates +
          updateK * contract.visitedCandidates := by
      rw [Nat.mul_comm contract.visitedCandidates lookupK,
        Nat.mul_comm contract.visitedCandidates updateK]
    _ = (lookupK + updateK) * contract.visitedCandidates := by
      rw [Nat.add_mul]
    _ = (lookupK + updateK) * contract.visitedCandidates + 0 := by
      omega

end IndexedDetectorOperationContract

private theorem linearBoundWith_compose_size {cost m n c k cm km : Nat}
    (hCost : LinearBoundWith cost m c k)
    (hSize : LinearBoundWith m n cm km) :
    LinearBoundWith cost n (c * cm) (c * km + k) := by
  unfold LinearBoundWith at hCost hSize ⊢
  calc
    cost ≤ c * m + k := hCost
    _ ≤ c * (cm * n + km) + k :=
      Nat.add_le_add_right (Nat.mul_le_mul_left c hSize) k
    _ = (c * cm) * n + (c * km + k) := by
      rw [Nat.mul_add, Nat.mul_assoc]
      omega

private theorem linearBound_compose_size {cost m n : Nat}
    (hCost : LinearBound cost m)
    (hSize : LinearBound m n) :
    LinearBound cost n := by
  rcases hCost with ⟨c, k, hCost⟩
  rcases hSize with ⟨cm, km, hSize⟩
  exact ⟨c * cm, c * km + k, linearBoundWith_compose_size hCost hSize⟩

/--
Bounds connecting traversal-frame materialization to the input graph.

These fields are assumptions to be passed to future theorems, not global facts
about the current implementation.
-/
structure TraversalFrameSizeContract {g : Graph}
    (frame : TraversalFrame g) where
  treeLinksLinear :
    LinearBound frame.treeLinks.length (GraphSizes.ofGraph g).total
  candidateStackLinear :
    LinearBound (Algorithms.Flubble.candidateStack frame).length
      (GraphSizes.ofGraph g).total

/--
Output-size reuse contract for stages that consume already-produced flubble or
hierarchy data rather than walking the graph again.
-/
structure OutputSizeReuseContract (sizes : DecompositionSizes) where
  flubblesBoundedByOutput :
    sizes.detectedFlubbleCount ≤ sizes.outputSize
  hierarchyBoundedByOutput :
    sizes.hierarchyOutputSize ≤ sizes.outputSize
  outputLinearInInput :
    LinearBound sizes.outputSize sizes.inputSize

namespace OutputSizeReuseContract

theorem flubbles_linear_in_input {sizes : DecompositionSizes}
    (h : OutputSizeReuseContract sizes) :
    LinearBound sizes.detectedFlubbleCount sizes.inputSize :=
  LinearBound.reuse_output_size h.flubblesBoundedByOutput h.outputLinearInInput

theorem hierarchy_linear_in_input {sizes : DecompositionSizes}
    (h : OutputSizeReuseContract sizes) :
    LinearBound sizes.hierarchyOutputSize sizes.inputSize :=
  LinearBound.reuse_output_size h.hierarchyBoundedByOutput h.outputLinearInInput

end OutputSizeReuseContract

/-!
## Conditional indexed extraction costs

The next definitions state the strongest linear theorem currently supported by
the verified Lean detector and hierarchy proofs.  The theorem is intentionally
conditional: lookup/update costs for the class-index table, scan bookkeeping,
boundary emission, and hierarchy construction are explicit assumptions.  This
does not assert that the current povu C++ or Rust implementation satisfies
those assumptions; a later bridge must prove any implementation-specific cost
model.
-/

def boundaryOutputSize {g : Graph}
    (frame : TraversalFrame g)
    (classes : Algorithms.Flubble.CycleClassAssignment g) : Nat :=
  (Algorithms.Flubble.detectFlubblesIndexed frame classes).length

def hierarchyConstructionOutputSize {g : Graph}
    (frame : TraversalFrame g)
    (classes : Algorithms.Flubble.CycleClassAssignment g) : Nat :=
  (Algorithms.FlubbleTree.buildHierarchy frame classes).nodes.length
    + (Algorithms.FlubbleTree.buildHierarchy frame classes).flatMetadataSlots

/--
Named cost components for verified Lean flubble extraction and flat hierarchy
construction.  These are proof-side cost slots, not an executable timing model.
-/
structure ExtractionPipelineCosts where
  indexedOperations : Nat
  scanBookkeeping : Nat
  boundaryEmission : Nat
  hierarchyConstruction : Nat
  deriving Repr, DecidableEq

namespace ExtractionPipelineCosts

def total (costs : ExtractionPipelineCosts) : Nat :=
  costs.indexedOperations
    + costs.scanBookkeeping
    + costs.boundaryEmission
    + costs.hierarchyConstruction

end ExtractionPipelineCosts

/--
Explicit assumptions needed to turn the indexed detector's operation counts
into an end-to-end extraction cost theorem.

* `indexOperations` carries the O(1) class-index lookup/update assumptions.
* `indexedOperations_cost_eq` ties the operation-cost slot to that index model.
* `indexVisited_eq_candidateStack` ties the modeled scan length to the frame's
  candidate stack.
* `scanBookkeepingLinear` covers per-candidate detector work not represented as
  class-index lookup/update operations.
* `boundaryEmissionLinearInBoundaries` covers work proportional to the emitted
  boundary list.
* `hierarchyConstructionLinearInOutput` covers finite-map/list work needed to
  construct the flat hierarchy from the emitted boundary and reference slots.
-/
structure ConditionalExtractionCostContract {g : Graph}
    (frame : TraversalFrame g)
    (classes : Algorithms.Flubble.CycleClassAssignment g)
    (costs : ExtractionPipelineCosts) where
  indexOperations : IndexedDetectorOperationContract
  indexedOperations_cost_eq :
    costs.indexedOperations = indexOperations.operationCost
  indexVisited_eq_candidateStack :
    indexOperations.visitedCandidates =
      (Algorithms.Flubble.candidateStack frame).length
  scanBookkeepingLinear :
    LinearBound costs.scanBookkeeping
      (Algorithms.Flubble.candidateStack frame).length
  boundaryEmissionLinearInBoundaries :
    LinearBound costs.boundaryEmission (boundaryOutputSize frame classes)
  hierarchyConstructionLinearInOutput :
    LinearBound costs.hierarchyConstruction
      (hierarchyConstructionOutputSize frame classes)

theorem boundaryOutputSize_linear_in_candidateStack {g : Graph}
    (frame : TraversalFrame g)
    (classes : Algorithms.Flubble.CycleClassAssignment g) :
    LinearBoundWith (boundaryOutputSize frame classes)
      (Algorithms.Flubble.candidateStack frame).length 1 0 := by
  unfold boundaryOutputSize LinearBoundWith
  simpa using
    Algorithms.Flubble.detectFlubblesIndexed_length_le_candidateStack_length
      frame classes

theorem hierarchyConstructionOutputSize_linear_in_candidateStack {g : Graph}
    {frame : TraversalFrame g}
    {classes : Algorithms.Flubble.CycleClassAssignment g}
    (hSupported :
      Algorithms.FlubbleTree.SupportedHierarchyInput
        (Algorithms.Flubble.candidateStack frame)
        (Algorithms.Flubble.detectFlubbles frame classes)) :
    LinearBoundWith (hierarchyConstructionOutputSize frame classes)
      (Algorithms.Flubble.candidateStack frame).length 3 0 := by
  unfold hierarchyConstructionOutputSize LinearBoundWith
  rw [Algorithms.FlubbleTree.buildHierarchy_nodes_length,
    Algorithms.FlubbleTree.buildHierarchy_flatMetadataSlots_length hSupported]
  have hCount :
      (Algorithms.Flubble.detectFlubbles frame classes).length ≤
        (Algorithms.Flubble.candidateStack frame).length :=
    Algorithms.Flubble.detectFlubbles_length_le_candidateStack_length
      frame classes
  calc
    (Algorithms.Flubble.detectFlubbles frame classes).length +
        ((Algorithms.Flubble.detectFlubbles frame classes).length +
          (Algorithms.Flubble.detectFlubbles frame classes).length)
        ≤ (Algorithms.Flubble.candidateStack frame).length +
            ((Algorithms.Flubble.candidateStack frame).length +
              (Algorithms.Flubble.candidateStack frame).length) :=
      Nat.add_le_add hCount (Nat.add_le_add hCount hCount)
    _ = 3 * (Algorithms.Flubble.candidateStack frame).length + 0 := by
      omega

/--
Conditional linear-time theorem for verified Lean flubble boundary extraction
and flat hierarchy construction.

The conclusion is linear in the candidate-stack length and also exposes the
composed output bounds used by the proof: emitted indexed boundaries are at
most one per candidate, and the hierarchy node/reference output is at most
three slots per candidate under the supported hierarchy-input condition.  The
theorem is not a full graph-decomposition runtime theorem and makes no claim
about the current povu C++ or Rust implementation runtime.
-/
theorem conditional_extraction_linear_in_candidateStack {g : Graph}
    {frame : TraversalFrame g}
    {classes : Algorithms.Flubble.CycleClassAssignment g}
    {costs : ExtractionPipelineCosts}
    (hSupported :
      Algorithms.FlubbleTree.SupportedHierarchyInput
        (Algorithms.Flubble.candidateStack frame)
        (Algorithms.Flubble.detectFlubbles frame classes))
    (hCost : ConditionalExtractionCostContract frame classes costs) :
    LinearBound costs.total
      (Algorithms.Flubble.candidateStack frame).length
      ∧ LinearBoundWith (boundaryOutputSize frame classes)
        (Algorithms.Flubble.candidateStack frame).length 1 0
      ∧ LinearBoundWith (hierarchyConstructionOutputSize frame classes)
        (Algorithms.Flubble.candidateStack frame).length 3 0 := by
  have hIndex :
      LinearBound costs.indexedOperations
        (Algorithms.Flubble.candidateStack frame).length := by
    simpa [hCost.indexedOperations_cost_eq,
      hCost.indexVisited_eq_candidateStack] using
      IndexedDetectorOperationContract.operationCost_linear_in_candidates
        hCost.indexOperations
  have hBoundarySizeWith :=
    boundaryOutputSize_linear_in_candidateStack frame classes
  have hBoundarySize : LinearBound (boundaryOutputSize frame classes)
      (Algorithms.Flubble.candidateStack frame).length :=
    LinearBound.with_constants hBoundarySizeWith
  have hBoundaryCost :
      LinearBound costs.boundaryEmission
        (Algorithms.Flubble.candidateStack frame).length :=
    linearBound_compose_size
      hCost.boundaryEmissionLinearInBoundaries hBoundarySize
  have hHierarchySizeWith :=
    hierarchyConstructionOutputSize_linear_in_candidateStack
      (g := g) (frame := frame) (classes := classes) hSupported
  have hHierarchySize :
      LinearBound (hierarchyConstructionOutputSize frame classes)
        (Algorithms.Flubble.candidateStack frame).length :=
    LinearBound.with_constants hHierarchySizeWith
  have hHierarchyCost :
      LinearBound costs.hierarchyConstruction
        (Algorithms.Flubble.candidateStack frame).length :=
    linearBound_compose_size
      hCost.hierarchyConstructionLinearInOutput hHierarchySize
  have hTotal :
      LinearBound costs.total
        (Algorithms.Flubble.candidateStack frame).length := by
    unfold ExtractionPipelineCosts.total
    exact LinearBound.seq
      (LinearBound.seq
        (LinearBound.seq hIndex hCost.scanBookkeepingLinear)
        hBoundaryCost)
      hHierarchyCost
  exact ⟨hTotal, hBoundarySizeWith, hHierarchySizeWith⟩

/--
Algorithms-facing linearity statement for the indexed flubble stage.

Under the standard random-access class-index convention captured by
`ConditionalExtractionCostContract`, the indexed detector, boundary emission,
and flat hierarchy construction are linear in the candidate-stack length.
-/
theorem indexed_flubble_stage_linear_in_candidateStack {g : Graph}
    {frame : TraversalFrame g}
    {classes : Algorithms.Flubble.CycleClassAssignment g}
    {costs : ExtractionPipelineCosts}
    (hSupported :
      Algorithms.FlubbleTree.SupportedHierarchyInput
        (Algorithms.Flubble.candidateStack frame)
        (Algorithms.Flubble.detectFlubbles frame classes))
    (hCost : ConditionalExtractionCostContract frame classes costs) :
    LinearBound costs.total
      (Algorithms.Flubble.candidateStack frame).length :=
  (conditional_extraction_linear_in_candidateStack
    (g := g)
    (frame := frame)
    (classes := classes)
    (costs := costs)
    hSupported
    hCost).1

/--
Graph-size corollary for the indexed flubble stage.

Once the candidate-stack length is bounded linearly by
`segments.length + edgeCount`, the indexed flubble detector/extractor and flat
hierarchy construction are linear in that graph-size measure.
-/
theorem indexed_flubble_stage_linear_in_graphSize {g : Graph}
    {frame : TraversalFrame g}
    {classes : Algorithms.Flubble.CycleClassAssignment g}
    {costs : ExtractionPipelineCosts}
    (hSupported :
      Algorithms.FlubbleTree.SupportedHierarchyInput
        (Algorithms.Flubble.candidateStack frame)
        (Algorithms.Flubble.detectFlubbles frame classes))
    (hCandidateStack :
      LinearBound (Algorithms.Flubble.candidateStack frame).length
        (GraphSizes.ofGraph g).total)
    (hCost : ConditionalExtractionCostContract frame classes costs) :
    LinearBound costs.total (GraphSizes.ofGraph g).total :=
  LinearBound.compose_size
    (indexed_flubble_stage_linear_in_candidateStack
      (g := g)
      (frame := frame)
      (classes := classes)
      (costs := costs)
      hSupported
      hCost)
    hCandidateStack

/--
Named costs for the flubble decomposition stages that downstream proofs may
compose.  The fields are intentionally abstract; this record is only a theorem
interface until bridged to a concrete implementation model.
-/
structure DecompositionStageCosts where
  classIndex : StageCost
  traversal : StageCost
  detection : StageCost
  hierarchy : StageCost
  outputReuse : StageCost
  deriving Repr, DecidableEq

namespace DecompositionStageCosts

def detectorPrefix (costs : DecompositionStageCosts) : Nat :=
  costs.classIndex.cost + costs.traversal.cost + costs.detection.cost

def hierarchySuffix (costs : DecompositionStageCosts) : Nat :=
  costs.hierarchy.cost + costs.outputReuse.cost

def total (costs : DecompositionStageCosts) : Nat :=
  costs.detectorPrefix + costs.hierarchySuffix

theorem detectorPrefix_linear_with {costs : DecompositionStageCosts}
    {n cIndex kIndex cTraversal kTraversal cDetection kDetection : Nat}
    (hIndex : LinearStageWith costs.classIndex n cIndex kIndex)
    (hTraversal : LinearStageWith costs.traversal n cTraversal kTraversal)
    (hDetection : LinearStageWith costs.detection n cDetection kDetection) :
    LinearBoundWith costs.detectorPrefix n
      ((cIndex + cTraversal) + cDetection)
      ((kIndex + kTraversal) + kDetection) := by
  unfold detectorPrefix
  exact LinearBoundWith.seq (LinearBoundWith.seq hIndex hTraversal) hDetection

theorem hierarchySuffix_linear_with {costs : DecompositionStageCosts}
    {n cHierarchy kHierarchy cReuse kReuse : Nat}
    (hHierarchy : LinearStageWith costs.hierarchy n cHierarchy kHierarchy)
    (hReuse : LinearStageWith costs.outputReuse n cReuse kReuse) :
    LinearBoundWith costs.hierarchySuffix n
      (cHierarchy + cReuse) (kHierarchy + kReuse) := by
  unfold hierarchySuffix
  exact LinearBoundWith.seq hHierarchy hReuse

theorem total_linear_with {costs : DecompositionStageCosts}
    {n cPrefix kPrefix cSuffix kSuffix : Nat}
    (hPrefix : LinearBoundWith costs.detectorPrefix n cPrefix kPrefix)
    (hSuffix : LinearBoundWith costs.hierarchySuffix n cSuffix kSuffix) :
    LinearBoundWith costs.total n
      (cPrefix + cSuffix) (kPrefix + kSuffix) := by
  unfold total
  exact LinearBoundWith.seq hPrefix hSuffix

end DecompositionStageCosts

end Flubble
end Complexity
end PovuLean
