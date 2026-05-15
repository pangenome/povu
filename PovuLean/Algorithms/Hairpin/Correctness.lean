import PovuLean.Algorithms.Hairpin.Detect

/-!
Correctness theorem for the Lean hairpin reference detector.

The theorem is about the verified Lean reference algorithm.  The remaining
conformance obligation for a later task is to prove that povu's concrete
reverse DFS bracket state supplies a `HairpinScanAssignment.Correct`.
-/

namespace PovuLean
namespace Algorithms
namespace Hairpin

open Core

theorem firstExtensionBeforeTerminal?_extends {g : Graph}
    {scan : HairpinScanAssignment g} {right : Link} {rest : List Link}
    (h : firstExtensionBeforeTerminal? scan rest = some right) :
    scan.ExtendsStem right := by
  induction rest with
  | nil =>
      simp [firstExtensionBeforeTerminal?] at h
  | cons candidate tail ih =>
      unfold firstExtensionBeforeTerminal? at h
      by_cases hTerm : scan.terminatesStem candidate = true
      · simp [hTerm] at h
      · by_cases hExtend : scan.extendsStem candidate = true
        · simp [hTerm, hExtend] at h
          cases h
          exact hExtend
        · simp [hTerm, hExtend] at h
          exact ih h

theorem firstExtensionBeforeTerminal?_mem {g : Graph}
    {scan : HairpinScanAssignment g} {right : Link} {rest : List Link}
    (h : firstExtensionBeforeTerminal? scan rest = some right) :
    right ∈ rest := by
  induction rest with
  | nil =>
      simp [firstExtensionBeforeTerminal?] at h
  | cons candidate tail ih =>
      unfold firstExtensionBeforeTerminal? at h
      by_cases hTerm : scan.terminatesStem candidate = true
      · simp [hTerm] at h
      · by_cases hExtend : scan.extendsStem candidate = true
        · simp [hTerm, hExtend] at h
          cases h
          exact List.Mem.head _
        · simp [hTerm, hExtend] at h
          exact List.Mem.tail _ (ih h)

theorem ordered_edges_mem_of_mem_stack {left right : Link}
    {stack : List Link}
    (hLeft : left ∈ stack) (hRight : right ∈ stack) :
    (Boundary.ordered left right).outerEdge ∈ stack
      ∧ (Boundary.ordered left right).innerEdge ∈ stack := by
  unfold Boundary.ordered
  by_cases hLt : left.id < right.id
  · simp [hLt, hLeft, hRight]
  · simp [hLt, hLeft, hRight]

theorem firstExtensionBeforeTerminal?_right_ne_left {g : Graph}
    {scan : HairpinScanAssignment g} {left right : Link} {rest : List Link}
    (hStackNodup : (left :: rest).Nodup)
    (hClose : firstExtensionBeforeTerminal? scan rest = some right) :
    left ≠ right := by
  intro hEq
  have hRightMem : right ∈ rest := firstExtensionBeforeTerminal?_mem hClose
  subst hEq
  exact (List.nodup_cons.mp hStackNodup).1 hRightMem

theorem canonicalScanBoundary_edges_mem {g : Graph}
    {scan : HairpinScanAssignment g} {stack : List Link}
    {boundary : Boundary}
    (hBoundary : CanonicalScanBoundary scan stack boundary) :
    boundary.outerEdge ∈ stack ∧ boundary.innerEdge ∈ stack := by
  induction hBoundary with
  | here hStart hClose =>
      have hRightMem := firstExtensionBeforeTerminal?_mem hClose
      exact ordered_edges_mem_of_mem_stack (List.Mem.head _) (List.Mem.tail _ hRightMem)
  | there tail ih =>
      exact ⟨List.Mem.tail _ ih.1, List.Mem.tail _ ih.2⟩

theorem canonicalScanBoundary_distinct {g : Graph}
    {scan : HairpinScanAssignment g} {stack : List Link}
    {boundary : Boundary}
    (hNodup : stack.Nodup)
    (hBoundary : CanonicalScanBoundary scan stack boundary) :
    boundary.outerEdge ≠ boundary.innerEdge := by
  induction hBoundary with
  | here hStart hClose =>
      have hNe := firstExtensionBeforeTerminal?_right_ne_left hNodup hClose
      exact Boundary.ordered_distinct hNe
  | there tail ih =>
      exact ih (List.nodup_cons.mp hNodup).2

theorem canonicalScanBoundary_canonicalIds {g : Graph}
    {scan : HairpinScanAssignment g} {stack : List Link}
    {boundary : Boundary}
    (hBoundary : CanonicalScanBoundary scan stack boundary) :
    boundary.CanonicalIds := by
  induction hBoundary with
  | here hStart hClose =>
      exact Boundary.ordered_canonicalIds
  | there tail ih =>
      exact ih

/--
Soundness: every reported boundary is a real black tree-edge pair, is canonical,
and has a certified structural hairpin shape.
-/
theorem detectHairpins_sound {g : Graph} {frame : TraversalFrame g}
    {scan : HairpinScanAssignment g}
    (hInput : SupportedInput g frame)
    (hScan : HairpinScanAssignment.Correct frame scan)
    {boundary : Boundary}
    (hMem : boundary ∈ detectHairpins frame scan) :
    IsPaperHairpinBoundary frame boundary := by
  unfold detectHairpins at hMem
  rw [mem_detectStack] at hMem
  have hEdgeMem := canonicalScanBoundary_edges_mem hMem
  have hOuterCandidate := mem_reverseCandidateStack.mp hEdgeMem.1
  have hInnerCandidate := mem_reverseCandidateStack.mp hEdgeMem.2
  have hDistinct :=
    canonicalScanBoundary_distinct (reverseCandidateStack_nodup hInput) hMem
  have hCanonicalIds := canonicalScanBoundary_canonicalIds hMem
  have hShape := hScan.sound hMem
  exact ⟨hOuterCandidate, hInnerCandidate, hDistinct, hCanonicalIds, hShape⟩

/-- Completeness: every supported paper hairpin boundary is reported. -/
theorem detectHairpins_complete {g : Graph} {frame : TraversalFrame g}
    {scan : HairpinScanAssignment g}
    (hScan : HairpinScanAssignment.Correct frame scan)
    {boundary : Boundary}
    (hBoundary : IsPaperHairpinBoundary frame boundary) :
    boundary ∈ detectHairpins frame scan := by
  unfold detectHairpins
  rw [mem_detectStack]
  exact hScan.complete hBoundary

/-- Duplicate-free output guarantee required before deterministic downstream emission. -/
theorem detectHairpins_canonical_noDuplicates {g : Graph}
    (frame : TraversalFrame g) (scan : HairpinScanAssignment g) :
    NoDuplicateBoundaries (detectHairpins frame scan) :=
  detectHairpins_noDuplicates frame scan

/-- Main bundled correctness theorem for supported core graph inputs. -/
theorem detectHairpins_correct {g : Graph} {frame : TraversalFrame g}
    {scan : HairpinScanAssignment g}
    (hInput : SupportedInput g frame)
    (hScan : HairpinScanAssignment.Correct frame scan) :
    (∀ boundary,
      boundary ∈ detectHairpins frame scan →
        IsPaperHairpinBoundary frame boundary)
      ∧
    (∀ boundary,
      IsPaperHairpinBoundary frame boundary →
        boundary ∈ detectHairpins frame scan)
      ∧
    NoDuplicateBoundaries (detectHairpins frame scan) :=
  ⟨ fun _ hMem => detectHairpins_sound hInput hScan hMem
  , fun _ hBoundary => detectHairpins_complete hScan hBoundary
  , detectHairpins_canonical_noDuplicates frame scan ⟩

/-- GFA semantic-input corollary of the core correctness theorem. -/
theorem detectHairpins_correct_for_gfa {doc : GFA.Document}
    {accepted : doc.Accepted} {frame : TraversalFrame doc.toGraph}
    {scan : HairpinScanAssignment doc.toGraph}
    (hInput : SupportedGFAInput doc accepted frame)
    (hScan : HairpinScanAssignment.Correct frame scan) :
    (∀ boundary,
      boundary ∈ detectHairpins frame scan →
        IsPaperHairpinBoundary frame boundary)
      ∧
    (∀ boundary,
      IsPaperHairpinBoundary frame boundary →
        boundary ∈ detectHairpins frame scan)
      ∧
    NoDuplicateBoundaries (detectHairpins frame scan) :=
  detectHairpins_correct (Flubble.SupportedGFAInput.toSupportedInput hInput) hScan

end Hairpin
end Algorithms
end PovuLean
