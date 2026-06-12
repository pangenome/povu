import PovuLean.Algorithms.Flubble.Detect

/-!
Correctness theorem for the Lean flubble reference detector.

The theorem is intentionally about the Lean reference algorithm, not the current
C++ implementation.  The remaining conformance gap for `lean4-conformance-harness`
is to show that povu's bracket/class stack supplies a `CycleClassAssignment`
that satisfies `CycleClassAssignment.Correct` and uses the same canonical stack
order as `candidateStack`.
-/

namespace PovuLean
namespace Algorithms
namespace Flubble

open Core

theorem firstSameClass?_sameClass {g : Graph}
    {classes : CycleClassAssignment g} {left right : Link} {rest : List Link}
    (h : firstSameClass? classes left rest = some right) :
    classes.sameClass left right := by
  induction rest with
  | nil =>
      simp [firstSameClass?] at h
  | cons candidate tail ih =>
      unfold firstSameClass? at h
      by_cases hSame : classes.sameClassBool left candidate = true
      · simp [hSame] at h
        cases h
        exact (CycleClassAssignment.sameClassBool_eq_true).mp hSame
      · simp [hSame] at h
        exact ih h

theorem closeAfterGap?_sameClass {g : Graph}
    {classes : CycleClassAssignment g} {left right : Link} {rest : List Link}
    (h : closeAfterGap? classes left rest = some right) :
    classes.sameClass left right := by
  cases rest with
  | nil =>
      simp [closeAfterGap?] at h
  | cons immediate tail =>
      unfold closeAfterGap? at h
      by_cases hSame : classes.sameClassBool left immediate = true
      · simp [hSame] at h
      · simp [hSame] at h
        exact firstSameClass?_sameClass h

theorem firstSameClass?_mem {g : Graph}
    {classes : CycleClassAssignment g} {left right : Link} {rest : List Link}
    (h : firstSameClass? classes left rest = some right) :
    right ∈ rest := by
  induction rest with
  | nil =>
      simp [firstSameClass?] at h
  | cons candidate tail ih =>
      unfold firstSameClass? at h
      by_cases hSame : classes.sameClassBool left candidate = true
      · simp [hSame] at h
        cases h
        exact List.Mem.head _
      · simp [hSame] at h
        exact List.Mem.tail _ (ih h)

theorem closeAfterGap?_mem {g : Graph}
    {classes : CycleClassAssignment g} {left right : Link} {rest : List Link}
    (h : closeAfterGap? classes left rest = some right) :
    right ∈ rest := by
  cases rest with
  | nil =>
      simp [closeAfterGap?] at h
  | cons immediate tail =>
      unfold closeAfterGap? at h
      by_cases hSame : classes.sameClassBool left immediate = true
      · simp [hSame] at h
      · simp [hSame] at h
        exact List.Mem.tail _ (firstSameClass?_mem h)

theorem ordered_edges_mem_of_mem_stack {left right : Link}
    {stack : List Link}
    (hLeft : left ∈ stack) (hRight : right ∈ stack) :
    (Boundary.ordered left right).openEdge ∈ stack
      ∧ (Boundary.ordered left right).closeEdge ∈ stack := by
  unfold Boundary.ordered
  by_cases hLt : left.id < right.id
  · simp [hLt, hLeft, hRight]
  · simp [hLt, hLeft, hRight]

theorem ordered_candidates_of_candidates {g : Graph}
    {frame : TraversalFrame g} {left right : Link}
    (hLeft : IsBoundaryCandidate frame left)
    (hRight : IsBoundaryCandidate frame right) :
    IsBoundaryCandidate frame (Boundary.ordered left right).openEdge
      ∧ IsBoundaryCandidate frame (Boundary.ordered left right).closeEdge := by
  unfold Boundary.ordered
  by_cases hLt : left.id < right.id
  · simp [hLt, hLeft, hRight]
  · simp [hLt, hLeft, hRight]

theorem ordered_sameClass {g : Graph}
    {classes : CycleClassAssignment g} {left right : Link}
    (hSame : classes.sameClass left right) :
    classes.sameClass (Boundary.ordered left right).openEdge
      (Boundary.ordered left right).closeEdge := by
  unfold Boundary.ordered
  by_cases hLt : left.id < right.id
  · simpa [hLt] using hSame
  · simp [hLt, CycleClassAssignment.sameClass]
    exact hSame.symm

theorem closeAfterGap?_right_ne_left {g : Graph}
    {classes : CycleClassAssignment g} {left right : Link} {rest : List Link}
    (hStackNodup : (left :: rest).Nodup)
    (hClose : closeAfterGap? classes left rest = some right) :
    left ≠ right := by
  intro hEq
  have hRightMem : right ∈ rest := closeAfterGap?_mem hClose
  subst hEq
  exact (List.nodup_cons.mp hStackNodup).1 hRightMem

theorem canonicalClassBoundary_edges_mem {g : Graph}
    {classes : CycleClassAssignment g} {stack : List Link}
    {boundary : Boundary}
    (hBoundary : CanonicalClassBoundary classes stack boundary) :
    boundary.openEdge ∈ stack ∧ boundary.closeEdge ∈ stack := by
  induction hBoundary with
  | here hClose =>
      have hRightMem := closeAfterGap?_mem hClose
      exact ordered_edges_mem_of_mem_stack (List.Mem.head _) (List.Mem.tail _ hRightMem)
  | there tail ih =>
      exact ⟨List.Mem.tail _ ih.1, List.Mem.tail _ ih.2⟩

theorem canonicalClassBoundary_sameClass {g : Graph}
    {classes : CycleClassAssignment g} {stack : List Link}
    {boundary : Boundary}
    (hBoundary : CanonicalClassBoundary classes stack boundary) :
    classes.sameClass boundary.openEdge boundary.closeEdge := by
  induction hBoundary with
  | here hClose =>
      exact ordered_sameClass (closeAfterGap?_sameClass hClose)
  | there tail ih =>
      exact ih

theorem canonicalClassBoundary_distinct {g : Graph}
    {classes : CycleClassAssignment g} {stack : List Link}
    {boundary : Boundary}
    (hNodup : stack.Nodup)
    (hBoundary : CanonicalClassBoundary classes stack boundary) :
    boundary.openEdge ≠ boundary.closeEdge := by
  induction hBoundary with
  | here hClose =>
      have hNe := closeAfterGap?_right_ne_left hNodup hClose
      exact Boundary.ordered_distinct hNe
  | there tail ih =>
      exact ih (List.nodup_cons.mp hNodup).2

theorem canonicalClassBoundary_canonicalIds {g : Graph}
    {classes : CycleClassAssignment g} {stack : List Link}
    {boundary : Boundary}
    (hBoundary : CanonicalClassBoundary classes stack boundary) :
    boundary.CanonicalIds := by
  induction hBoundary with
  | here hClose =>
      exact Boundary.ordered_canonicalIds
  | there tail ih =>
      exact ih

/--
Soundness: every reported boundary is a real black tree-edge pair, is canonical,
and its two edges are cycle equivalent.
-/
theorem detectFlubbles_sound {g : Graph} {frame : TraversalFrame g}
    {classes : CycleClassAssignment g}
    (hInput : SupportedInput g frame)
    (hClasses : CycleClassAssignment.Correct frame classes)
    {boundary : Boundary}
    (hMem : boundary ∈ detectFlubbles frame classes) :
    IsFlubbleBoundary frame classes boundary := by
  unfold detectFlubbles at hMem
  rw [mem_detectStack] at hMem
  have hEdgeMem := canonicalClassBoundary_edges_mem hMem
  have hOpenCandidate := mem_candidateStack.mp hEdgeMem.1
  have hCloseCandidate := mem_candidateStack.mp hEdgeMem.2
  have hSame := canonicalClassBoundary_sameClass hMem
  have hCycle := hClasses.sound hOpenCandidate hCloseCandidate hSame
  have hDistinct :=
    canonicalClassBoundary_distinct hInput.candidateStack_nodup hMem
  have hCanonicalIds := canonicalClassBoundary_canonicalIds hMem
  exact
    ⟨ ⟨hOpenCandidate, hCloseCandidate, hDistinct, hCanonicalIds, hCycle⟩
    , hMem ⟩

/-- Completeness: every supported canonical flubble boundary is reported. -/
theorem detectFlubbles_complete {g : Graph} {frame : TraversalFrame g}
    {classes : CycleClassAssignment g}
    {boundary : Boundary}
    (hBoundary : IsFlubbleBoundary frame classes boundary) :
    boundary ∈ detectFlubbles frame classes := by
  unfold detectFlubbles
  rw [mem_detectStack]
  exact hBoundary.2

/-- Duplicate/canonical output guarantee required before deterministic VCF emission. -/
theorem detectFlubbles_canonical_noDuplicates {g : Graph}
    (frame : TraversalFrame g) (classes : CycleClassAssignment g) :
    NoDuplicateBoundaries (detectFlubbles frame classes) :=
  detectFlubbles_noDuplicates frame classes

/-- Main bundled correctness theorem for supported core graph inputs. -/
theorem detectFlubbles_correct {g : Graph} {frame : TraversalFrame g}
    {classes : CycleClassAssignment g}
    (hInput : SupportedInput g frame)
    (hClasses : CycleClassAssignment.Correct frame classes) :
    (∀ boundary,
      boundary ∈ detectFlubbles frame classes →
        IsFlubbleBoundary frame classes boundary)
      ∧
    (∀ boundary,
      IsFlubbleBoundary frame classes boundary →
        boundary ∈ detectFlubbles frame classes)
      ∧
    NoDuplicateBoundaries (detectFlubbles frame classes) :=
  ⟨ fun _ hMem => detectFlubbles_sound hInput hClasses hMem
  , fun _ hBoundary => detectFlubbles_complete hBoundary
  , detectFlubbles_canonical_noDuplicates frame classes ⟩

/--
Canonical flubble count bound for supported core graph inputs.

This is the paper-facing linear count theorem for the current Lean detector
surface: the deduplicated list of canonical `IsFlubbleBoundary` witnesses
emitted by `detectFlubbles` has length at most the graph edge count.  It does
not count all pairwise cycle-equivalent edge pairs.
-/
theorem flubble_count_le_numEdges {g : Graph} {frame : TraversalFrame g}
    {classes : CycleClassAssignment g}
    (hInput : SupportedInput g frame)
    (_hClasses : CycleClassAssignment.Correct frame classes) :
    (detectFlubbles frame classes).length ≤ g.edgeCount :=
  detectFlubbles_length_le_supportedInput_graph_edgeCount classes hInput

/-- GFA semantic-input corollary of the core correctness theorem. -/
theorem detectFlubbles_correct_for_gfa {doc : GFA.Document}
    {accepted : doc.Accepted} {frame : TraversalFrame doc.toGraph}
    {classes : CycleClassAssignment doc.toGraph}
    (hInput : SupportedGFAInput doc accepted frame)
    (hClasses : CycleClassAssignment.Correct frame classes) :
    (∀ boundary,
      boundary ∈ detectFlubbles frame classes →
        IsFlubbleBoundary frame classes boundary)
      ∧
    (∀ boundary,
      IsFlubbleBoundary frame classes boundary →
        boundary ∈ detectFlubbles frame classes)
      ∧
    NoDuplicateBoundaries (detectFlubbles frame classes) :=
  detectFlubbles_correct hInput.toSupportedInput hClasses

end Flubble
end Algorithms
end PovuLean
