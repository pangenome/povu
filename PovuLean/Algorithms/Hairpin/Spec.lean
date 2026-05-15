import PovuLean.Algorithms.Hairpin.InputInvariant

/-!
Paper-facing hairpin inversion boundary specification.

The paper describes a hairpin as a loop preceded and followed by the same stem,
with the latter traversed in opposite directions.  In this supported Lean scope
we make that convention explicit with `ReverseStemConvention`: an outbound stem
walk must be the reverse-complement traversal of the inbound stem walk.

The executable detector below does not rebuild the C++ bracket list.  Instead it
scans the reverse real-black candidate stack using a certified
`HairpinScanAssignment`.  The certificate is the boundary between the eventual
DFS/bracket conformance proof and the reference algorithm proved here.
-/

namespace PovuLean
namespace Algorithms
namespace Hairpin

open Core

/-- Boundary pair reported for a supported hairpin stem. -/
structure Boundary where
  outerEdge : Link
  innerEdge : Link
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

/-- Stable key used by deterministic downstream emitters. -/
def key (boundary : Boundary) : Nat × Nat :=
  (boundary.outerEdge.id, boundary.innerEdge.id)

/-- Deterministic orientation used by this reference detector. -/
def CanonicalIds (boundary : Boundary) : Prop :=
  boundary.outerEdge.id ≤ boundary.innerEdge.id

/-- Canonicalize a pair by edge id while keeping both real boundary records. -/
def ordered (left right : Link) : Boundary :=
  if left.id < right.id then
    { outerEdge := left, innerEdge := right }
  else
    { outerEdge := right, innerEdge := left }

theorem ordered_canonicalIds {left right : Link} :
    (ordered left right).CanonicalIds := by
  unfold ordered CanonicalIds
  by_cases hLt : left.id < right.id
  · simp [hLt, Nat.le_of_lt hLt]
  · simp [hLt]
    exact Nat.le_of_not_gt hLt

theorem ordered_distinct {left right : Link}
    (hNe : left ≠ right) :
    (ordered left right).outerEdge ≠ (ordered left right).innerEdge := by
  unfold ordered
  by_cases hLt : left.id < right.id
  · simpa [hLt] using hNe
  · simp [hLt]
    exact hNe.symm

end Boundary

/-- A supported one-black-edge cut witness for the loop side of a hairpin. -/
def OneBlackEdgeCut (g : Graph) : Prop :=
  ∃ cut : EdgeCut g,
    cut.edges.length = 1
      ∧ cut.AllBlack
      ∧ cut.AllReal

/-- The loop component of a supported hairpin: an edge-simple cycle plus a one-black-edge cut. -/
def HairpinLoop (g : Graph) (cycle : CarrierCycle g) : Prop :=
  cycle.edges.Nodup ∧ OneBlackEdgeCut g

/--
The reverse-traversal convention for the stem.  It is separated from biological
reverse-complement sequence claims; this theorem layer is graph-structural.
-/
def ReverseStemConvention {g : Graph}
    (inbound outbound : Walk g) : Prop :=
  outbound.edges = IsWalk.reverseEdges inbound.edges
    ∧ outbound.start = inbound.finish.reverse
    ∧ outbound.finish = inbound.start.reverse

namespace ReverseStemConvention

theorem of_walk_reverse {g : Graph}
    (hSym : g.EdgeSymmetric) (inbound : Walk g) :
    ReverseStemConvention inbound (inbound.reverse hSym) := by
  exact ⟨rfl, rfl, rfl⟩

end ReverseStemConvention

/-- The cycle-free stem plus the two reported boundary edges. -/
def HairpinStem {g : Graph} (boundary : Boundary)
    (inbound outbound : Walk g) : Prop :=
  inbound.edges ≠ []
    ∧ inbound.orientedVertices.Nodup
    ∧ ReverseStemConvention inbound outbound
    ∧ boundary.outerEdge ∈ inbound.edges
    ∧ boundary.innerEdge ∈ inbound.edges

/--
Structural hairpin witness in the supported core graph model: a loop attached
to an inbound stem whose outbound traversal is the reverse-complement of that
same stem.
-/
def HasHairpinShape (g : Graph) (boundary : Boundary) : Prop :=
  ∃ loop : CarrierCycle g,
  ∃ inbound outbound : Walk g,
    HairpinLoop g loop
      ∧ HairpinStem boundary inbound outbound
      ∧ inbound.finish = loop.root

/-- Paper-level real hairpin boundary facts, independent of output membership. -/
def IsPaperHairpinBoundary {g : Graph} (frame : TraversalFrame g)
    (boundary : Boundary) : Prop :=
  IsBoundaryCandidate frame boundary.outerEdge
    ∧ IsBoundaryCandidate frame boundary.innerEdge
    ∧ boundary.outerEdge ≠ boundary.innerEdge
    ∧ boundary.CanonicalIds
    ∧ HasHairpinShape g boundary

/--
Find the first stem-extension edge before a terminal event.  A terminal event
models the C++ scan reaching a leaf/root before the hairpin end was extended,
which makes the candidate a near miss rather than a report.
-/
def firstExtensionBeforeTerminal? {g : Graph}
    (scan : HairpinScanAssignment g) : List Link → Option Link
  | [] => none
  | candidate :: rest =>
      if scan.terminatesStem candidate then
        none
      else if scan.extendsStem candidate then
        some candidate
      else
        firstExtensionBeforeTerminal? scan rest

/-- `right` is the first extension found before termination in `rest`. -/
def FirstExtensionBeforeTerminal {g : Graph}
    (scan : HairpinScanAssignment g) (rest : List Link) (right : Link) : Prop :=
  firstExtensionBeforeTerminal? scan rest = some right

/--
Canonical reverse-scan relation.  The head edge starts a stem, and the first
later extension before termination closes it.
-/
inductive CanonicalScanBoundary {g : Graph}
    (scan : HairpinScanAssignment g) :
    List Link → Boundary → Prop where
  | here {left right : Link} {rest : List Link}
      (hStart : scan.StartsStem left)
      (hClose : firstExtensionBeforeTerminal? scan rest = some right) :
      CanonicalScanBoundary scan (left :: rest) (Boundary.ordered left right)
  | there {left : Link} {rest : List Link} {boundary : Boundary}
      (tail : CanonicalScanBoundary scan rest boundary) :
      CanonicalScanBoundary scan (left :: rest) boundary

/--
Supported hairpin boundary: a structural hairpin plus the canonical reverse
scan rule over the certified stack.
-/
def IsHairpinBoundary {g : Graph} (frame : TraversalFrame g)
    (scan : HairpinScanAssignment g) (boundary : Boundary) : Prop :=
  IsPaperHairpinBoundary frame boundary
    ∧ CanonicalScanBoundary scan (reverseCandidateStack frame) boundary

namespace HairpinScanAssignment

/--
Correctness contract for the scan assignment.  Later bracket-stack conformance
work should prove that povu's concrete reverse DFS state supplies this contract.
-/
structure Correct {g : Graph} (frame : TraversalFrame g)
    (scan : HairpinScanAssignment g) : Prop where
  sound :
    ∀ {boundary : Boundary},
      CanonicalScanBoundary scan (reverseCandidateStack frame) boundary →
      HasHairpinShape g boundary
  complete :
    ∀ {boundary : Boundary},
      IsPaperHairpinBoundary frame boundary →
      CanonicalScanBoundary scan (reverseCandidateStack frame) boundary

end HairpinScanAssignment

/-- Downstream deterministic emitters need list-level duplicate freedom. -/
def NoDuplicateBoundaries (boundaries : List Boundary) : Prop :=
  boundaries.Nodup

end Hairpin
end Algorithms
end PovuLean
