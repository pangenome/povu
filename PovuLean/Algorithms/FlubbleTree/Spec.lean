import PovuLean.Algorithms.Flubble.Correctness

/-!
Specification layer for the flubble hierarchy.

The base flubble proof reports canonical boundary pairs.  This module gives the
next layer a stack-span semantics: each boundary occupies the interval between
its two endpoints in the certified candidate stack, and parenthood means strict
span containment.  Non-laminar pairs are recorded as unsupported hierarchy
inputs rather than folded into an arbitrary tree.
-/

namespace PovuLean
namespace Algorithms
namespace FlubbleTree

open Core

abbrev Boundary := Flubble.Boundary

/-- Inclusive endpoint positions of one flubble boundary in a candidate stack. -/
structure Span where
  start : Nat
  finish : Nat
  deriving Repr, DecidableEq

namespace Span

def WellFormed (span : Span) : Prop :=
  span.start < span.finish

def Contains (outer inner : Span) : Prop :=
  outer.start < inner.start ∧ inner.finish < outer.finish

def Disjoint (left right : Span) : Prop :=
  left.finish < right.start ∨ right.finish < left.start

def width (span : Span) : Nat :=
  span.finish - span.start

def containsBool (outer inner : Span) : Bool :=
  decide (outer.start < inner.start) &&
    decide (inner.finish < outer.finish)

theorem containsBool_eq_true {outer inner : Span} :
    outer.containsBool inner = true ↔ Span.Contains outer inner := by
  unfold containsBool
  rw [Bool.and_eq_true]
  constructor
  · intro h
    exact
      ⟨(decide_eq_true_eq (p := outer.start < inner.start)).mp h.1,
        (decide_eq_true_eq (p := inner.finish < outer.finish)).mp h.2⟩
  · intro h
    exact
      ⟨(decide_eq_true_eq (p := outer.start < inner.start)).mpr h.1,
        (decide_eq_true_eq (p := inner.finish < outer.finish)).mpr h.2⟩

end Span

/-- First position of `needle` in a list, if present. -/
def indexOf? [DecidableEq α] (needle : α) : List α → Option Nat
  | [] => none
  | head :: tail =>
      if needle = head then
        some 0
      else
        Option.map Nat.succ (indexOf? needle tail)

/--
Span occupied by a boundary in a stack.  The endpoint order is normalized by
position; equal or missing endpoints produce `none`.
-/
def boundarySpan? (stack : List Link) (boundary : Boundary) : Option Span :=
  match indexOf? boundary.openEdge stack, indexOf? boundary.closeEdge stack with
  | some openIndex, some closeIndex =>
      if openIndex < closeIndex then
        some { start := openIndex, finish := closeIndex }
      else if closeIndex < openIndex then
        some { start := closeIndex, finish := openIndex }
      else
        none
  | _, _ => none

def HasSpan (stack : List Link) (boundary : Boundary) (span : Span) : Prop :=
  boundarySpan? stack boundary = some span

def BoundariesHaveSpans (stack : List Link) (boundaries : List Boundary) : Prop :=
  ∀ boundary, boundary ∈ boundaries → ∃ span, HasSpan stack boundary span

def NestedIn (stack : List Link) (parent child : Boundary) : Prop :=
  ∃ parentSpan childSpan,
    HasSpan stack parent parentSpan
      ∧ HasSpan stack child childSpan
      ∧ Span.Contains parentSpan childSpan

def DisjointIn (stack : List Link) (left right : Boundary) : Prop :=
  ∃ leftSpan rightSpan,
    HasSpan stack left leftSpan
      ∧ HasSpan stack right rightSpan
      ∧ Span.Disjoint leftSpan rightSpan

def LaminarPair (stack : List Link) (left right : Boundary) : Prop :=
  NestedIn stack left right ∨ NestedIn stack right left ∨ DisjointIn stack left right

def LaminarBoundaries (stack : List Link) (boundaries : List Boundary) : Prop :=
  ∀ left, left ∈ boundaries →
    ∀ right, right ∈ boundaries →
      left ≠ right → LaminarPair stack left right

/--
Supported hierarchy inputs are exactly the base flubble outputs whose endpoints
all appear in the stack and whose spans are laminar.
-/
def SupportedHierarchyInput (stack : List Link) (boundaries : List Boundary) : Prop :=
  BoundariesHaveSpans stack boundaries ∧ LaminarBoundaries stack boundaries

/--
Unsupported graph/output case for this task: at least one pair cannot be placed
in a laminar nesting relation.  This covers crossing intervals and malformed
boundary endpoints from an external flubble provider.
-/
def UnsupportedNonLaminarBoundaries
    (stack : List Link) (boundaries : List Boundary) : Prop :=
  ∃ left, left ∈ boundaries ∧
    ∃ right, right ∈ boundaries ∧
      left ≠ right ∧ ¬ LaminarPair stack left right

def UnsupportedMissingSpan (stack : List Link) (boundaries : List Boundary) : Prop :=
  ∃ boundary, boundary ∈ boundaries ∧ boundarySpan? stack boundary = none

def UnsupportedHierarchyInput (stack : List Link) (boundaries : List Boundary) : Prop :=
  UnsupportedMissingSpan stack boundaries
    ∨ UnsupportedNonLaminarBoundaries stack boundaries

structure Node where
  boundary : Boundary
  parent? : Option Boundary
  deriving Repr, DecidableEq

structure Hierarchy where
  nodes : List Node
  deriving Repr, DecidableEq

namespace Hierarchy

def boundaries (hierarchy : Hierarchy) : List Boundary :=
  hierarchy.nodes.map Node.boundary

end Hierarchy

/-- Parent interpretation for one hierarchy node. -/
def NodeParentCorrect
    (stack : List Link) (boundaries : List Boundary) (node : Node) : Prop :=
  node.boundary ∈ boundaries ∧
    match node.parent? with
    | none =>
        ¬ ∃ parent, parent ∈ boundaries ∧ NestedIn stack parent node.boundary
    | some parent =>
        parent ∈ boundaries ∧ NestedIn stack parent node.boundary

/--
Forest invariant for the flat parent-map representation: nodes correspond
exactly to the supplied boundaries, boundaries are duplicate-free, the input is
laminar, and every stored parent is the canonical nesting parent selected by the
builder.
-/
def IsForestHierarchy
    (stack : List Link) (boundaries : List Boundary)
    (hierarchy : Hierarchy) : Prop :=
  hierarchy.boundaries = boundaries
    ∧ hierarchy.boundaries.Nodup
    ∧ SupportedHierarchyInput stack boundaries
    ∧ ∀ node, node ∈ hierarchy.nodes →
        NodeParentCorrect stack boundaries node

def CorrespondsToFlubbles {g : Graph}
    (frame : TraversalFrame g) (classes : Flubble.CycleClassAssignment g)
    (hierarchy : Hierarchy) : Prop :=
  ∀ boundary,
    boundary ∈ hierarchy.boundaries ↔
      Flubble.IsFlubbleBoundary frame classes boundary

/--
Output contract needed by downstream emitters: the hierarchy preserves the base
detector's deterministic boundary order, keeps duplicate freedom, and every
boundary has the flubble-level canonical endpoint orientation.
-/
def OutputOrderCanonical {g : Graph}
    (frame : TraversalFrame g) (classes : Flubble.CycleClassAssignment g)
    (hierarchy : Hierarchy) : Prop :=
  hierarchy.boundaries = Flubble.detectFlubbles frame classes
    ∧ hierarchy.boundaries.Nodup
    ∧ ∀ boundary, boundary ∈ hierarchy.boundaries →
        boundary.CanonicalIds

def IsCorrectHierarchy {g : Graph}
    (frame : TraversalFrame g) (classes : Flubble.CycleClassAssignment g)
    (hierarchy : Hierarchy) : Prop :=
  IsForestHierarchy
      (Flubble.candidateStack frame)
      (Flubble.detectFlubbles frame classes)
      hierarchy
    ∧ CorrespondsToFlubbles frame classes hierarchy
    ∧ OutputOrderCanonical frame classes hierarchy

theorem supportedHierarchyInput_not_unsupported {stack : List Link}
    {boundaries : List Boundary}
    (hSupported : SupportedHierarchyInput stack boundaries) :
    ¬ UnsupportedHierarchyInput stack boundaries := by
  intro hUnsupported
  cases hUnsupported with
  | inl hMissing =>
      rcases hMissing with ⟨boundary, hMem, hNone⟩
      rcases hSupported.1 boundary hMem with ⟨span, hSpan⟩
      unfold HasSpan at hSpan
      rw [hSpan] at hNone
      cases hNone
  | inr hNonLaminar =>
      rcases hNonLaminar with
        ⟨left, hLeft, right, hRight, hNe, hNotLaminar⟩
      exact hNotLaminar (hSupported.2 left hLeft right hRight hNe)

end FlubbleTree
end Algorithms
end PovuLean
