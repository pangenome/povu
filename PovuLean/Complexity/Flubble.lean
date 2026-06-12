import PovuLean.Algorithms.FlubbleTree.Spec
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
