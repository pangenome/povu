import PovuLean.Algorithms.Flubble.InputInvariant

/-!
Paper-facing flubble boundary specification.

Within the supported scope, a reported flubble is represented by the two real
black tree edges that open and close one cycle-equivalence class in the
canonical stack order.  Hairpin detection, hierarchy construction, PVST records,
and VCF route data are intentionally outside this module.
-/

namespace PovuLean
namespace Algorithms
namespace Flubble

open Core

/-- Boundary pair reported by the base flubble detector. -/
structure Boundary where
  openEdge : Link
  closeEdge : Link
  deriving Repr, DecidableEq

namespace Boundary

instance : BEq Boundary where
  beq left right := decide (left = right)

instance : LawfulBEq Boundary where
  rfl := by
    intro boundary
    exact decide_eq_true rfl
  eq_of_beq := by
    intro left right h
    exact of_decide_eq_true h

/-- Stable key used by downstream deterministic emitters. -/
def key (boundary : Boundary) : Nat × Nat :=
  (boundary.openEdge.id, boundary.closeEdge.id)

/-- The detector always emits boundaries in this normalized orientation. -/
def CanonicalIds (boundary : Boundary) : Prop :=
  boundary.openEdge.id ≤ boundary.closeEdge.id

/-- Canonicalize a pair by edge id for deterministic output. -/
def ordered (left right : Link) : Boundary :=
  if left.id < right.id then
    { openEdge := left, closeEdge := right }
  else
    { openEdge := right, closeEdge := left }

theorem ordered_canonicalIds {left right : Link} :
    (ordered left right).CanonicalIds := by
  unfold ordered CanonicalIds
  by_cases hLt : left.id < right.id
  · simp [hLt, Nat.le_of_lt hLt]
  · simp [hLt]
    exact Nat.le_of_not_gt hLt

theorem ordered_distinct {left right : Link}
    (hNe : left ≠ right) :
    (ordered left right).openEdge ≠ (ordered left right).closeEdge := by
  unfold ordered
  by_cases hLt : left.id < right.id
  · simpa [hLt] using hNe
  · simp [hLt]
    exact hNe.symm

end Boundary

namespace CycleClassAssignment

/-- The executable class comparison used by the reference detector. -/
def sameClassBool {g : Graph} (classes : CycleClassAssignment g)
    (left right : Link) : Bool :=
  decide (classes.sameClass left right)

theorem sameClassBool_eq_true {g : Graph}
    {classes : CycleClassAssignment g} {left right : Link} :
    classes.sameClassBool left right = true ↔ classes.sameClass left right := by
  unfold sameClassBool
  rw [decide_eq_true_eq (p := classes.sameClass left right)]

end CycleClassAssignment

/--
Finds the first later edge in the same class.  The surrounding detector only
uses this after skipping at least one candidate, matching the non-adjacent
open/close convention used by povu's equivalence-class stack.
-/
def firstSameClass? {g : Graph} (classes : CycleClassAssignment g)
    (left : Link) : List Link → Option Link
  | [] => none
  | candidate :: rest =>
      if classes.sameClassBool left candidate then
        some candidate
      else
        firstSameClass? classes left rest

/--
The next closing edge for `left`, requiring at least one intervening candidate.
If the immediately following candidate is already in the same class, `left`
does not open a flubble in this canonical base detector.
-/
def closeAfterGap? {g : Graph} (classes : CycleClassAssignment g)
    (left : Link) : List Link → Option Link
  | [] => none
  | immediate :: rest =>
      if classes.sameClassBool left immediate then
        none
      else
        firstSameClass? classes left rest

/-- `right` is the first later same-class edge in `rest`. -/
def FirstSameClassInRest {g : Graph} (classes : CycleClassAssignment g)
    (left : Link) (rest : List Link) (right : Link) : Prop :=
  firstSameClass? classes left rest = some right

/--
`right` closes `left` after at least one intervening candidate in the canonical
stack order.
-/
def ClosesAfterGap {g : Graph} (classes : CycleClassAssignment g)
    (left : Link) (rest : List Link) (right : Link) : Prop :=
  closeAfterGap? classes left rest = some right

/--
Canonical class-boundary relation in a candidate stack.  This is the open/close
rule that resolves classes with more than two eligible edges into deterministic
base flubble boundaries.
-/
inductive CanonicalClassBoundary {g : Graph}
    (classes : CycleClassAssignment g) :
    List Link → Boundary → Prop where
  | here {left right : Link} {rest : List Link}
      (hClose : closeAfterGap? classes left rest = some right) :
      CanonicalClassBoundary classes (left :: rest) (Boundary.ordered left right)
  | there {left : Link} {rest : List Link} {boundary : Boundary}
      (tail : CanonicalClassBoundary classes rest boundary) :
      CanonicalClassBoundary classes (left :: rest) boundary

/-- Paper-level real/cycle-equivalent boundary facts, independent of output. -/
def IsPaperBoundary {g : Graph} (frame : TraversalFrame g)
    (boundary : Boundary) : Prop :=
  IsBoundaryCandidate frame boundary.openEdge
    ∧ IsBoundaryCandidate frame boundary.closeEdge
    ∧ boundary.openEdge ≠ boundary.closeEdge
    ∧ boundary.CanonicalIds
    ∧ CycleEquivalentEdges g boundary.openEdge boundary.closeEdge

/--
Supported flubble boundary: a paper-level pair plus the canonical open/close
rule over the certified class stack.
-/
def IsFlubbleBoundary {g : Graph} (frame : TraversalFrame g)
    (classes : CycleClassAssignment g) (boundary : Boundary) : Prop :=
  IsPaperBoundary frame boundary
    ∧ CanonicalClassBoundary classes (candidateStack frame) boundary

/-- Downstream deterministic emitters need list-level duplicate freedom. -/
def NoDuplicateBoundaries (boundaries : List Boundary) : Prop :=
  boundaries.Nodup

end Flubble
end Algorithms
end PovuLean
