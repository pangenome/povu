import PovuLean.Algorithms.Flubble.Correctness
import PovuLean.Core.Examples
import PovuLean.GFA.Examples

/-!
Checked examples for the verified flubble reference detector.

These examples exercise the executable canonical stack behavior independently
of the not-yet-implemented Rust/Lean conformance harness.
-/

namespace PovuLean
namespace Algorithms
namespace Flubble
namespace Examples

open Core
open Core.Examples

def blackEdge (id reverseId : EdgeId) (source target : OrientedSegment) : Link :=
  Link.between id reverseId source target EdgeColor.black EdgeProvenance.real

def simpleA : Link := blackEdge 0 8 (fwd 0) (fwd 1)
def simpleX : Link := blackEdge 1 9 (fwd 1) (fwd 2)
def simpleZ : Link := blackEdge 2 10 (fwd 2) (fwd 3)
def simpleArc : Link := Link.between 3 11 (fwd 3) (fwd 0)
def simpleArcRc : Link := simpleArc.reverse
def simpleARc : Link := simpleA.reverse
def simpleXRc : Link := simpleX.reverse
def simpleZRc : Link := simpleZ.reverse

def simpleFlubbleGraph : Graph :=
  { segments := [seg 0 "S", seg 1 "A", seg 2 "T", seg 3 "J"]
    links :=
      [ simpleA, simpleX, simpleZ, simpleArc
      , simpleARc, simpleXRc, simpleZRc, simpleArcRc ] }

example : simpleFlubbleGraph.WellFormed := by
  native_decide

example : simpleFlubbleGraph.BidirectedChecked := by
  native_decide

def simpleFrame : TraversalFrame simpleFlubbleGraph :=
  { roots := [fwd 0]
    discovered := [fwd 0, fwd 1, fwd 2, fwd 3]
    treeLinks := [simpleA, simpleX, simpleZ]
    roots_valid := by
      intro v h
      simp at h
      cases h
      native_decide
    discovered_valid := by
      intro v h
      simp at h
      rcases h with h | h | h | h <;> cases h <;> native_decide
    treeLinks_supported := by
      intro e h
      simp at h
      rcases h with h | h | h <;> cases h <;> native_decide }

def simpleClasses : CycleClassAssignment simpleFlubbleGraph :=
  { classOf := fun edge =>
      if edge.id = simpleA.id then 7
      else if edge.id = simpleZ.id then 7
      else edge.id + 100 }

def simpleBoundary : Boundary :=
  { openEdge := simpleA, closeEdge := simpleZ }

example : detectFlubbles simpleFrame simpleClasses = [simpleBoundary] := by
  native_decide

example : simpleBoundary ∈ detectFlubbles simpleFrame simpleClasses := by
  native_decide

def negativeClasses : CycleClassAssignment simpleFlubbleGraph :=
  { classOf := fun edge => edge.id }

example : detectFlubbles simpleFrame negativeClasses = [] := by
  native_decide

def nestedOuterOpen : Link := blackEdge 0 20 (fwd 0) (fwd 1)
def nestedInnerOpen : Link := blackEdge 1 21 (fwd 1) (fwd 2)
def nestedMiddle : Link := blackEdge 2 22 (fwd 2) (fwd 3)
def nestedInnerClose : Link := blackEdge 3 23 (fwd 3) (fwd 4)
def nestedOuterClose : Link := blackEdge 4 24 (fwd 4) (fwd 5)
def nestedReturn : Link := Link.between 5 25 (fwd 5) (fwd 0)

def nestedGraph : Graph :=
  { segments :=
      [ seg 0 "A", seg 1 "B", seg 2 "C", seg 3 "D", seg 4 "E", seg 5 "F" ]
    links :=
      [ nestedOuterOpen
      , nestedInnerOpen
      , nestedMiddle
      , nestedInnerClose
      , nestedOuterClose
      , nestedReturn
      , nestedOuterOpen.reverse
      , nestedInnerOpen.reverse
      , nestedMiddle.reverse
      , nestedInnerClose.reverse
      , nestedOuterClose.reverse
      , nestedReturn.reverse ] }

example : nestedGraph.WellFormed := by
  native_decide

example : nestedGraph.BidirectedChecked := by
  native_decide

def nestedFrame : TraversalFrame nestedGraph :=
  { roots := [fwd 0]
    discovered := [fwd 0, fwd 1, fwd 2, fwd 3, fwd 4, fwd 5]
    treeLinks :=
      [ nestedOuterOpen
      , nestedInnerOpen
      , nestedMiddle
      , nestedInnerClose
      , nestedOuterClose ]
    roots_valid := by
      intro v h
      simp at h
      cases h
      native_decide
    discovered_valid := by
      intro v h
      simp at h
      rcases h with h | h | h | h | h | h <;> cases h <;> native_decide
    treeLinks_supported := by
      intro e h
      simp at h
      rcases h with h | h | h | h | h <;> cases h <;> native_decide }

def nestedClasses : CycleClassAssignment nestedGraph :=
  { classOf := fun edge =>
      if edge.id = nestedOuterOpen.id then 10
      else if edge.id = nestedOuterClose.id then 10
      else if edge.id = nestedInnerOpen.id then 20
      else if edge.id = nestedInnerClose.id then 20
      else edge.id + 100 }

def nestedInnerBoundary : Boundary :=
  { openEdge := nestedInnerOpen, closeEdge := nestedInnerClose }

def nestedOuterBoundary : Boundary :=
  { openEdge := nestedOuterOpen, closeEdge := nestedOuterClose }

example :
    detectFlubbles nestedFrame nestedClasses =
      [nestedOuterBoundary, nestedInnerBoundary] := by
  native_decide

example : nestedOuterBoundary ∈ detectFlubbles nestedFrame nestedClasses := by
  native_decide

example : nestedInnerBoundary ∈ detectFlubbles nestedFrame nestedClasses := by
  native_decide

def overlappingClasses : CycleClassAssignment nestedGraph :=
  { classOf := fun edge =>
      if edge.id = nestedOuterOpen.id then 30
      else if edge.id = nestedInnerClose.id then 30
      else if edge.id = nestedInnerOpen.id then 40
      else if edge.id = nestedOuterClose.id then 40
      else edge.id + 100 }

def overlappingLeftBoundary : Boundary :=
  { openEdge := nestedOuterOpen, closeEdge := nestedInnerClose }

def overlappingRightBoundary : Boundary :=
  { openEdge := nestedInnerOpen, closeEdge := nestedOuterClose }

/-- Crossing class pairs exercise the overlapping-stack case. -/
example :
    detectFlubbles nestedFrame overlappingClasses =
      [overlappingLeftBoundary, overlappingRightBoundary] := by
  native_decide

def adjacentOpen : Link := blackEdge 0 30 (fwd 0) (fwd 1)
def adjacentClose : Link := blackEdge 1 31 (fwd 1) (fwd 2)

def adjacentGraph : Graph :=
  { segments := [seg 0 "A", seg 1 "B", seg 2 "C"]
    links :=
      [ adjacentOpen
      , adjacentClose
      , adjacentOpen.reverse
      , adjacentClose.reverse ] }

def adjacentFrame : TraversalFrame adjacentGraph :=
  { roots := [fwd 0]
    discovered := [fwd 0, fwd 1, fwd 2]
    treeLinks := [adjacentOpen, adjacentClose]
    roots_valid := by
      intro v h
      simp at h
      cases h
      native_decide
    discovered_valid := by
      intro v h
      simp at h
      rcases h with h | h | h <;> cases h <;> native_decide
    treeLinks_supported := by
      intro e h
      simp at h
      rcases h with h | h <;> cases h <;> native_decide }

def adjacentClasses : CycleClassAssignment adjacentGraph :=
  { classOf := fun _ => 1 }

/-- Adjacent same-class edges are not a base flubble boundary in this scope. -/
example : detectFlubbles adjacentFrame adjacentClasses = [] := by
  native_decide

example :
    NoDuplicateBoundaries (detectFlubbles nestedFrame nestedClasses) := by
  exact detectFlubbles_noDuplicates nestedFrame nestedClasses

end Examples
end Flubble
end Algorithms
end PovuLean
