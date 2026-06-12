import PovuLean.Complexity.Flubble
import PovuLean.Pipeline

/-!
# Conditional semantic decomposition linear-time composition

This module states and proves the full graph-decomposition composition theorem
for the semantic Lean pipeline only.  It does not claim that the current povu
C++ or Rust implementation runs in linear time; that requires a later
machine-checkable implementation bridge such as `bridge-povu-implementation`.

The theorem below is intentionally conditional.  The following contracts remain
external to Lean in this batch and appear as named theorem hypotheses:

* `componentDecompositionLinear`: component decomposition cost is linear in the
  input graph-size measure.
* `tipDummyAugmentationLinear`: tip/dummy augmentation cost is linear in the
  component-decomposition graph-size measure.
* `traversalFrameConstructionLinear`: DFS traversal-frame construction cost is
  linear in the augmented graph-size measure.
* `cycleClassAssignmentLinear`: cycle-class assignment cost is linear in the
  traversal-frame-size measure.
* `hairpinScanLinear`: hairpin scan cost is linear in the candidate-stack
  measure.
* `ConditionalExtractionCostContract`: indexed flubble extraction, boundary
  emission, and flubble hierarchy construction costs satisfy the explicit
  operation and output-size contracts from `PovuLean.Complexity.Flubble`.
* `SemanticDecomposeIntermediateSizeBounds`: intermediate component,
  augmentation, traversal-frame, candidate-stack, and hairpin-output sizes are
  each linear in the input graph size.

Proved in this batch: the arithmetic composition that turns those contracts and
the upstream flubble extraction theorem into one input-linear semantic
decomposition cost bound.  Implementation/conformance obligations are limited
to supplying these contracts for a concrete povu execution and relating that
execution to the semantic Lean witnesses used here.
-/

namespace PovuLean
namespace Complexity
namespace Decomposition

open Core
open Algorithms

/--
Intermediate size measures that must be bounded to avoid hidden superlinear
outputs in the full semantic decomposition composition theorem.
-/
structure SemanticDecomposeIntermediateSizes {g : Graph}
    (frame : TraversalFrame g) where
  inputGraphSize : Nat
  componentGraphSize : Nat
  augmentedGraphSize : Nat
  traversalFrameSize : Nat
  candidateStackSize : Nat
  hairpinBoundaryCount : Nat
  deriving Repr, DecidableEq

/--
Named intermediate-size bounds required by the composition theorem.

These are theorem hypotheses, not global facts about the executable
specification or current implementation.
-/
structure SemanticDecomposeIntermediateSizeBounds {g : Graph}
    {frame : TraversalFrame g}
    (sizes : SemanticDecomposeIntermediateSizes frame) where
  inputGraphSize_eq :
    sizes.inputGraphSize = (Flubble.GraphSizes.ofGraph g).total
  candidateStackSize_eq :
    sizes.candidateStackSize =
      (Algorithms.Flubble.candidateStack frame).length
  componentGraphLinear :
    LinearBound sizes.componentGraphSize sizes.inputGraphSize
  augmentedGraphLinear :
    LinearBound sizes.augmentedGraphSize sizes.inputGraphSize
  traversalFrameLinear :
    LinearBound sizes.traversalFrameSize sizes.inputGraphSize
  candidateStackLinear :
    LinearBound sizes.candidateStackSize sizes.inputGraphSize
  hairpinBoundaryCountLinear :
    LinearBound sizes.hairpinBoundaryCount sizes.inputGraphSize

/--
Named proof-side cost slots for the semantic graph-decomposition pipeline.

`flubbleExtractionAndHierarchy` reuses the upstream extraction cost model.  Its
fields explicitly include indexed detector operations, scan bookkeeping,
boundary emission, and flubble hierarchy construction.
-/
structure SemanticDecomposePipelineCosts where
  componentDecomposition : Flubble.StageCost
  tipDummyAugmentation : Flubble.StageCost
  traversalFrameConstruction : Flubble.StageCost
  cycleClassAssignment : Flubble.StageCost
  hairpinScan : Flubble.StageCost
  flubbleExtractionAndHierarchy : Flubble.ExtractionPipelineCosts
  deriving Repr, DecidableEq

namespace SemanticDecomposePipelineCosts

def externalPrefix (costs : SemanticDecomposePipelineCosts) : Nat :=
  costs.componentDecomposition.cost
    + costs.tipDummyAugmentation.cost
    + costs.traversalFrameConstruction.cost
    + costs.cycleClassAssignment.cost
    + costs.hairpinScan.cost

def total (costs : SemanticDecomposePipelineCosts) : Nat :=
  costs.externalPrefix + costs.flubbleExtractionAndHierarchy.total

end SemanticDecomposePipelineCosts

/--
Explicit stage contracts for the full semantic decomposition pipeline.

Each field names the size measure appropriate to that stage.  The theorem uses
`SemanticDecomposeIntermediateSizeBounds` to compose those measures back to the
input graph size.
-/
structure SemanticDecomposeStageContracts {g : Graph}
    {frame : TraversalFrame g}
    (classes : Algorithms.Flubble.CycleClassAssignment g)
    (sizes : SemanticDecomposeIntermediateSizes frame)
    (costs : SemanticDecomposePipelineCosts) where
  componentDecompositionLinear :
    Flubble.LinearStage costs.componentDecomposition sizes.inputGraphSize
  tipDummyAugmentationLinear :
    Flubble.LinearStage costs.tipDummyAugmentation sizes.componentGraphSize
  traversalFrameConstructionLinear :
    Flubble.LinearStage costs.traversalFrameConstruction sizes.augmentedGraphSize
  cycleClassAssignmentLinear :
    Flubble.LinearStage costs.cycleClassAssignment sizes.traversalFrameSize
  hairpinScanLinear :
    Flubble.LinearStage costs.hairpinScan sizes.candidateStackSize
  flubbleExtractionAndHierarchy :
    Flubble.ConditionalExtractionCostContract
      frame
      classes
      costs.flubbleExtractionAndHierarchy

/--
Conditional full graph-decomposition linear-time theorem for the semantic Lean
pipeline.

The conclusion is a proof-side cost statement over `costs.total`, linear in the
input graph-size measure `segments.length + edgeCount`.  The theorem uses the
upstream flubble extraction theorem for the flubble boundary and hierarchy
portion, while every unproved paper/implementation stage remains a named
contract hypothesis.  It also returns the semantic VCF correctness theorem to
make clear that this is the Lean semantic pipeline boundary, not a runtime claim
about current povu C++ or Rust code.
-/
theorem conditional_semantic_decompose_linear {doc : GFA.Document}
    {accepted : doc.Accepted}
    {frame : TraversalFrame doc.toGraph}
    {classes : Algorithms.Flubble.CycleClassAssignment doc.toGraph}
    {scan : Algorithms.Hairpin.HairpinScanAssignment doc.toGraph}
    {calls : List VCF.VariantCall}
    {sizes : SemanticDecomposeIntermediateSizes frame}
    {costs : SemanticDecomposePipelineCosts}
    (hInput : Algorithms.Flubble.SupportedGFAInput doc accepted frame)
    (hClasses : Algorithms.Flubble.CycleClassAssignment.Correct frame classes)
    (hScan : Algorithms.Hairpin.HairpinScanAssignment.Correct frame scan)
    (hHierarchySupported :
      Algorithms.FlubbleTree.SupportedHierarchyInput
        (Algorithms.Flubble.candidateStack frame)
        (Algorithms.Flubble.detectFlubbles frame classes))
    (hCalls :
      VCF.ReferenceCallSet
        frame
        scan
        (Algorithms.FlubbleTree.buildHierarchy frame classes)
        calls)
    (hSizeBounds : SemanticDecomposeIntermediateSizeBounds sizes)
    (hContracts : SemanticDecomposeStageContracts classes sizes costs) :
    LinearBound costs.total sizes.inputGraphSize
      ∧ VCF.EmissionCorrect
        frame
        classes
        scan
        (Algorithms.FlubbleTree.buildHierarchy frame classes)
        calls
        (VCF.emitRecords calls) := by
  have hTip :
      LinearBound costs.tipDummyAugmentation.cost sizes.inputGraphSize :=
    LinearBound.compose_size
      hContracts.tipDummyAugmentationLinear
      hSizeBounds.componentGraphLinear
  have hTraversal :
      LinearBound costs.traversalFrameConstruction.cost sizes.inputGraphSize :=
    LinearBound.compose_size
      hContracts.traversalFrameConstructionLinear
      hSizeBounds.augmentedGraphLinear
  have hCycle :
      LinearBound costs.cycleClassAssignment.cost sizes.inputGraphSize :=
    LinearBound.compose_size
      hContracts.cycleClassAssignmentLinear
      hSizeBounds.traversalFrameLinear
  have hHairpin :
      LinearBound costs.hairpinScan.cost sizes.inputGraphSize :=
    LinearBound.compose_size
      hContracts.hairpinScanLinear
      hSizeBounds.candidateStackLinear
  have hExtractionCandidate :
      LinearBound costs.flubbleExtractionAndHierarchy.total
        sizes.candidateStackSize :=
    by
      rw [hSizeBounds.candidateStackSize_eq]
      exact
        (Flubble.conditional_extraction_linear_in_candidateStack
          (g := doc.toGraph)
          (frame := frame)
          (classes := classes)
          (costs := costs.flubbleExtractionAndHierarchy)
          hHierarchySupported
          hContracts.flubbleExtractionAndHierarchy).1
  have hExtraction :
      LinearBound costs.flubbleExtractionAndHierarchy.total sizes.inputGraphSize :=
    LinearBound.compose_size
      hExtractionCandidate
      hSizeBounds.candidateStackLinear
  have hExternal :
      LinearBound costs.externalPrefix sizes.inputGraphSize := by
    unfold SemanticDecomposePipelineCosts.externalPrefix
    exact LinearBound.seq
      (LinearBound.seq
        (LinearBound.seq
          (LinearBound.seq
            hContracts.componentDecompositionLinear
            hTip)
          hTraversal)
        hCycle)
      hHairpin
  have hTotal :
      LinearBound costs.total sizes.inputGraphSize := by
    unfold SemanticDecomposePipelineCosts.total
    exact LinearBound.seq hExternal hExtraction
  have hSemantic :
      VCF.EmissionCorrect
        frame
        classes
        scan
        (Algorithms.FlubbleTree.buildHierarchy frame classes)
        calls
        (VCF.emitRecords calls) :=
    Pipeline.semanticGfaToVcf_correct
      hInput
      hClasses
      hScan
      hHierarchySupported
      hCalls
  exact ⟨hTotal, hSemantic⟩

end Decomposition
end Complexity
end PovuLean
