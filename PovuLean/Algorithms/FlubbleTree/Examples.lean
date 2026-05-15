import PovuLean.Algorithms.FlubbleTree.Correctness
import PovuLean.Algorithms.Flubble.Examples

/-!
Checked examples for the flubble hierarchy builder.

The examples exercise hierarchy construction over already verified flubble
outputs: one flubble, nested flubbles, sibling flubbles, no-flubble output, and
an empty graph/frame.
-/

namespace PovuLean
namespace Algorithms
namespace FlubbleTree
namespace Examples

open Core
open Core.Examples


def singleHierarchy : Hierarchy :=
  buildHierarchy Flubble.Examples.simpleFrame Flubble.Examples.simpleClasses

example :
    singleHierarchy.nodes =
      [ { boundary := Flubble.Examples.simpleBoundary, parent? := none } ] := by
  native_decide

def nestedHierarchy : Hierarchy :=
  buildHierarchy Flubble.Examples.nestedFrame Flubble.Examples.nestedClasses

example :
    nestedHierarchy.nodes =
      [ { boundary := Flubble.Examples.nestedOuterBoundary, parent? := none }
      , { boundary := Flubble.Examples.nestedInnerBoundary,
          parent? := some Flubble.Examples.nestedOuterBoundary } ] := by
  native_decide

example :
    canonicalParent?
      (Flubble.candidateStack Flubble.Examples.nestedFrame)
      (Flubble.detectFlubbles Flubble.Examples.nestedFrame Flubble.Examples.nestedClasses)
      Flubble.Examples.nestedInnerBoundary =
        some Flubble.Examples.nestedOuterBoundary := by
  native_decide

def siblingOpenLeft : Link := Flubble.Examples.blackEdge 0 30 (fwd 0) (fwd 1)
def siblingMiddleLeft : Link := Flubble.Examples.blackEdge 1 31 (fwd 1) (fwd 2)
def siblingCloseLeft : Link := Flubble.Examples.blackEdge 2 32 (fwd 2) (fwd 3)
def siblingGap : Link := Flubble.Examples.blackEdge 3 33 (fwd 3) (fwd 4)
def siblingOpenRight : Link := Flubble.Examples.blackEdge 4 34 (fwd 4) (fwd 5)
def siblingMiddleRight : Link := Flubble.Examples.blackEdge 5 35 (fwd 5) (fwd 6)
def siblingCloseRight : Link := Flubble.Examples.blackEdge 6 36 (fwd 6) (fwd 7)
def siblingReturn : Link := Link.between 7 37 (fwd 7) (fwd 0)

def siblingGraph : Graph :=
  { segments :=
      [ seg 0 "A", seg 1 "B", seg 2 "C", seg 3 "D"
      , seg 4 "E", seg 5 "F", seg 6 "G", seg 7 "H" ]
    links :=
      [ siblingOpenLeft
      , siblingMiddleLeft
      , siblingCloseLeft
      , siblingGap
      , siblingOpenRight
      , siblingMiddleRight
      , siblingCloseRight
      , siblingReturn
      , siblingOpenLeft.reverse
      , siblingMiddleLeft.reverse
      , siblingCloseLeft.reverse
      , siblingGap.reverse
      , siblingOpenRight.reverse
      , siblingMiddleRight.reverse
      , siblingCloseRight.reverse
      , siblingReturn.reverse ] }

example : siblingGraph.WellFormed := by
  native_decide

example : siblingGraph.BidirectedChecked := by
  native_decide

def siblingFrame : TraversalFrame siblingGraph :=
  { roots := [fwd 0]
    discovered := [fwd 0, fwd 1, fwd 2, fwd 3, fwd 4, fwd 5, fwd 6, fwd 7]
    treeLinks :=
      [ siblingOpenLeft
      , siblingMiddleLeft
      , siblingCloseLeft
      , siblingGap
      , siblingOpenRight
      , siblingMiddleRight
      , siblingCloseRight ]
    roots_valid := by
      intro v h
      simp at h
      cases h
      native_decide
    discovered_valid := by
      intro v h
      simp at h
      rcases h with h | h | h | h | h | h | h | h <;>
        cases h <;> native_decide
    treeLinks_supported := by
      intro e h
      simp at h
      rcases h with h | h | h | h | h | h | h <;>
        cases h <;> native_decide }

def siblingClasses : Flubble.CycleClassAssignment siblingGraph :=
  { classOf := fun edge =>
      if edge.id = siblingOpenLeft.id then 10
      else if edge.id = siblingCloseLeft.id then 10
      else if edge.id = siblingOpenRight.id then 20
      else if edge.id = siblingCloseRight.id then 20
      else edge.id + 100 }

def siblingLeftBoundary : Boundary :=
  { openEdge := siblingOpenLeft, closeEdge := siblingCloseLeft }

def siblingRightBoundary : Boundary :=
  { openEdge := siblingOpenRight, closeEdge := siblingCloseRight }

example :
    Flubble.detectFlubbles siblingFrame siblingClasses =
      [siblingLeftBoundary, siblingRightBoundary] := by
  native_decide

def siblingHierarchy : Hierarchy :=
  buildHierarchy siblingFrame siblingClasses

example :
    siblingHierarchy.nodes =
      [ { boundary := siblingLeftBoundary, parent? := none }
      , { boundary := siblingRightBoundary, parent? := none } ] := by
  native_decide

def noFlubbleHierarchy : Hierarchy :=
  buildHierarchy Flubble.Examples.adjacentFrame Flubble.Examples.adjacentClasses

example : noFlubbleHierarchy.nodes = [] := by
  native_decide

def emptyGraph : Graph :=
  { segments := []
    links := [] }

def emptyFrame : TraversalFrame emptyGraph :=
  { roots := []
    discovered := []
    treeLinks := []
    roots_valid := by
      intro v h
      cases h
    discovered_valid := by
      intro v h
      cases h
    treeLinks_supported := by
      intro e h
      cases h }

def emptyClasses : Flubble.CycleClassAssignment emptyGraph :=
  { classOf := fun edge => edge.id }

def emptyHierarchy : Hierarchy :=
  buildHierarchy emptyFrame emptyClasses

example : Flubble.detectFlubbles emptyFrame emptyClasses = [] := by
  native_decide

example : emptyHierarchy.nodes = [] := by
  native_decide

end Examples
end FlubbleTree
end Algorithms
end PovuLean
