import PovuLean.Algorithms.FlubbleTree.Spec

/-!
Executable hierarchy builder for verified flubble boundaries.

The builder is intentionally a pure transformation over already verified
flubble outputs.  It does not discover flubbles; it only assigns each reported
boundary the nearest containing boundary in the candidate-stack span order.
-/

namespace PovuLean
namespace Algorithms
namespace FlubbleTree

open Core

def nestedBool (stack : List Link) (parent child : Boundary) : Bool :=
  match boundarySpan? stack parent, boundarySpan? stack child with
  | some parentSpan, some childSpan =>
      parentSpan.containsBool childSpan
  | _, _ => false

theorem nestedBool_eq_true {stack : List Link} {parent child : Boundary} :
    nestedBool stack parent child = true ↔ NestedIn stack parent child := by
  unfold nestedBool NestedIn HasSpan
  cases hParent : boundarySpan? stack parent with
  | none =>
      simp
  | some parentSpan =>
      cases hChild : boundarySpan? stack child with
      | none =>
          simp
      | some childSpan =>
          simp [Span.containsBool_eq_true]

def boundaryKeyLt (left right : Boundary) : Bool :=
  if left.openEdge.id < right.openEdge.id then
    true
  else if right.openEdge.id < left.openEdge.id then
    false
  else
    left.closeEdge.id < right.closeEdge.id

/--
Deterministic tie-breaker between two already-known parent candidates.  Smaller
span wins; equal spans fall back to the flubble boundary key.
-/
def betterParent (stack : List Link) (left right : Boundary) : Boundary :=
  match boundarySpan? stack left, boundarySpan? stack right with
  | some leftSpan, some rightSpan =>
      if leftSpan.width < rightSpan.width then
        left
      else if rightSpan.width < leftSpan.width then
        right
      else if boundaryKeyLt left right then
        left
      else
        right
  | some _, none => left
  | none, some _ => right
  | none, none =>
      if boundaryKeyLt left right then left else right

theorem betterParent_eq_left_or_right {stack : List Link}
    {left right : Boundary} :
    betterParent stack left right = left ∨
      betterParent stack left right = right := by
  unfold betterParent
  cases boundarySpan? stack left <;>
    cases boundarySpan? stack right <;>
    repeat
      first
      | split
      | simp

/-- Select the nearest containing parent from a deterministic boundary list. -/
def chooseNearestParent? (stack : List Link) (child : Boundary) :
    List Boundary → Option Boundary
  | [] => none
  | candidate :: rest =>
      let tailParent? := chooseNearestParent? stack child rest
      if nestedBool stack candidate child then
        match tailParent? with
        | none => some candidate
        | some tailParent =>
            some (betterParent stack candidate tailParent)
      else
        tailParent?

def canonicalParent? (stack : List Link) (boundaries : List Boundary)
    (child : Boundary) : Option Boundary :=
  chooseNearestParent? stack child boundaries

def nodeFor (stack : List Link) (boundaries : List Boundary)
    (boundary : Boundary) : Node :=
  { boundary := boundary
    parent? := canonicalParent? stack boundaries boundary }

/-- Build the parent-map representation from a stack and a boundary list. -/
def buildHierarchyFrom (stack : List Link) (boundaries : List Boundary) :
    Hierarchy :=
  { nodes := boundaries.map (nodeFor stack boundaries) }

/-- Build the hierarchy directly from the verified flubble detector output. -/
def buildHierarchy {g : Graph} (frame : TraversalFrame g)
    (classes : Flubble.CycleClassAssignment g) : Hierarchy :=
  buildHierarchyFrom
    (Flubble.candidateStack frame)
    (Flubble.detectFlubbles frame classes)

end FlubbleTree
end Algorithms
end PovuLean
