import PovuLean.Algorithms.Hairpin.Correctness
import PovuLean.Core.Examples
import PovuLean.GFA.Examples

/-!
Checked examples for the verified hairpin reference detector.

The examples exercise positive reporting, near misses that still contain a loop
but do not satisfy the reverse scan, and the reverse traversal convention used
for the stem.
-/

namespace PovuLean
namespace Algorithms
namespace Hairpin
namespace Examples

open Core
open Core.Examples

def blackEdge (id reverseId : EdgeId) (source target : OrientedSegment) : Link :=
  Link.between id reverseId source target EdgeColor.black EdgeProvenance.real

def greyEdge (id reverseId : EdgeId) (source target : OrientedSegment) : Link :=
  Link.between id reverseId source target EdgeColor.grey EdgeProvenance.real

def stemOuter : Link := blackEdge 0 10 (fwd 0) (fwd 1)
def stemInner : Link := blackEdge 1 11 (fwd 1) (fwd 2)
def terminalEdge : Link := blackEdge 2 12 (fwd 2) (fwd 4)
def loopForward : Link := greyEdge 3 13 (fwd 2) (fwd 3)
def loopBack : Link := greyEdge 4 14 (fwd 3) (fwd 2)

def hairpinGraph : Graph :=
  { segments :=
      [ seg 0 "O", seg 1 "S", seg 2 "L", seg 3 "C", seg 4 "T" ]
    links :=
      [ stemOuter
      , stemInner
      , terminalEdge
      , loopForward
      , loopBack
      , stemOuter.reverse
      , stemInner.reverse
      , terminalEdge.reverse
      , loopForward.reverse
      , loopBack.reverse ] }

example : hairpinGraph.WellFormed := by
  native_decide

example : hairpinGraph.BidirectedChecked := by
  native_decide

def hairpinEdgeSymmetric : hairpinGraph.EdgeSymmetric :=
  GFA.edgeSymmetric_of_bidirectedChecked
    (by native_decide : hairpinGraph.BidirectedChecked)

def positiveFrame : TraversalFrame hairpinGraph :=
  { roots := [fwd 0]
    discovered := [fwd 0, fwd 1, fwd 2, fwd 3, fwd 4]
    treeLinks := [stemOuter, stemInner]
    roots_valid := by
      intro v h
      simp at h
      cases h
      native_decide
    discovered_valid := by
      intro v h
      simp at h
      rcases h with h | h | h | h | h <;> cases h <;> native_decide
    treeLinks_supported := by
      intro e h
      simp at h
      rcases h with h | h <;> cases h <;> native_decide }

def terminalNearMissFrame : TraversalFrame hairpinGraph :=
  { roots := [fwd 0]
    discovered := [fwd 0, fwd 1, fwd 2, fwd 3, fwd 4]
    treeLinks := [stemOuter, terminalEdge, stemInner]
    roots_valid := by
      intro v h
      simp at h
      cases h
      native_decide
    discovered_valid := by
      intro v h
      simp at h
      rcases h with h | h | h | h | h <;> cases h <;> native_decide
    treeLinks_supported := by
      intro e h
      simp at h
      rcases h with h | h | h <;> cases h <;> native_decide }

def positiveScan : HairpinScanAssignment hairpinGraph :=
  { startsStem := fun edge => edge.id == stemInner.id
    extendsStem := fun edge => edge.id == stemOuter.id
    terminatesStem := fun _ => false }

def noStartNearMissScan : HairpinScanAssignment hairpinGraph :=
  { startsStem := fun _ => false
    extendsStem := fun edge => edge.id == stemOuter.id
    terminatesStem := fun _ => false }

def terminalNearMissScan : HairpinScanAssignment hairpinGraph :=
  { startsStem := fun edge => edge.id == stemInner.id
    extendsStem := fun edge => edge.id == stemOuter.id
    terminatesStem := fun edge => edge.id == terminalEdge.id }

def positiveBoundary : Boundary :=
  { outerEdge := stemOuter, innerEdge := stemInner }

/-- The detector scans the DFS candidate stack in reverse order. -/
example : reverseCandidateStack positiveFrame = [stemInner, stemOuter] := by
  native_decide

example :
    reverseCandidateStack terminalNearMissFrame =
      [stemInner, terminalEdge, stemOuter] := by
  native_decide

example : detectHairpins positiveFrame positiveScan = [positiveBoundary] := by
  native_decide

example : positiveBoundary ∈ detectHairpins positiveFrame positiveScan := by
  native_decide

/-- A loop alone is not enough: without a bracket-empty start, nothing reports. -/
example : detectHairpins positiveFrame noStartNearMissScan = [] := by
  native_decide

/-- A terminal event before the extension is a near miss, not a hairpin report. -/
example :
    detectHairpins terminalNearMissFrame terminalNearMissScan = [] := by
  native_decide

def positiveLoopCycle : Cycle hairpinGraph :=
  { root := fwd 2
    edges := [loopForward, loopBack]
    valid :=
      IsWalk.cons loopForward rfl rfl (by native_decide) (by native_decide)
        (IsWalk.cons loopBack rfl rfl (by native_decide) (by native_decide)
          (IsWalk.nil (by native_decide)))
    nonempty := by native_decide }

def positiveLoopCut : EdgeCut hairpinGraph :=
  { edges := [stemInner]
    supported := by
      intro e h
      simp at h
      cases h
      native_decide }

def positiveLoopWitness : HairpinLoop hairpinGraph positiveLoopCycle := by
  constructor
  · native_decide
  · refine ⟨positiveLoopCut, rfl, ?_, ?_⟩
    · intro e h
      simp [positiveLoopCut] at h
      cases h
      rfl
    · intro e h
      simp [positiveLoopCut] at h
      cases h
      rfl

def positiveStemInbound : Walk hairpinGraph :=
  { start := fwd 0
    finish := fwd 2
    edges := [stemOuter, stemInner]
    valid :=
      IsWalk.cons stemOuter rfl rfl (by native_decide) (by native_decide)
        (IsWalk.cons stemInner rfl rfl (by native_decide) (by native_decide)
          (IsWalk.nil (by native_decide))) }

def positiveStemOutbound : Walk hairpinGraph :=
  positiveStemInbound.reverse hairpinEdgeSymmetric

example :
    positiveStemOutbound.edges =
      IsWalk.reverseEdges positiveStemInbound.edges := rfl

example : positiveStemOutbound.start = positiveStemInbound.finish.reverse := rfl

example : positiveStemOutbound.finish = positiveStemInbound.start.reverse := rfl

def positiveStemWitness :
    HairpinStem positiveBoundary positiveStemInbound positiveStemOutbound := by
  constructor
  · native_decide
  constructor
  · native_decide
  constructor
  · exact ReverseStemConvention.of_walk_reverse
      hairpinEdgeSymmetric positiveStemInbound
  constructor
  · native_decide
  · native_decide

def positiveHairpinShape : HasHairpinShape hairpinGraph positiveBoundary := by
  refine ⟨positiveLoopCycle, positiveStemInbound, positiveStemOutbound, ?_, ?_, rfl⟩
  · exact positiveLoopWitness
  · exact positiveStemWitness

example : IsPaperHairpinBoundary positiveFrame positiveBoundary := by
  exact
    ⟨ by native_decide
    , by native_decide
    , by native_decide
    , by
        unfold Boundary.CanonicalIds positiveBoundary stemOuter stemInner
        exact Nat.zero_le 1
    , positiveHairpinShape ⟩

example :
    NoDuplicateBoundaries (detectHairpins positiveFrame positiveScan) := by
  exact detectHairpins_noDuplicates positiveFrame positiveScan

end Examples
end Hairpin
end Algorithms
end PovuLean
