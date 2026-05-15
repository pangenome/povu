import PovuLean.Algorithms.Hairpin.Spec

/-!
Executable Lean reference detector for supported hairpin boundaries.

The detector traverses the reverse real-black candidate stack.  A start event
opens a possible stem; the first later extension event before a terminal event
emits the canonical boundary pair.  Exact duplicate boundaries are removed.
-/

namespace PovuLean
namespace Algorithms
namespace Hairpin

open Core

/-- Insert a boundary only if it is not already present. -/
def insertUniqueBoundary (boundary : Boundary) (boundaries : List Boundary) :
    List Boundary :=
  if boundary ∈ boundaries then
    boundaries
  else
    boundary :: boundaries

theorem mem_insertUniqueBoundary {boundary existing : Boundary}
    {boundaries : List Boundary} :
    existing ∈ insertUniqueBoundary boundary boundaries ↔
      existing = boundary ∨ existing ∈ boundaries := by
  unfold insertUniqueBoundary
  by_cases hMem : boundary ∈ boundaries
  · simp [hMem]
    intro hEq
    exact hEq ▸ hMem
  · simp [hMem]

theorem nodup_insertUniqueBoundary {boundary : Boundary}
    {boundaries : List Boundary}
    (h : boundaries.Nodup) :
    (insertUniqueBoundary boundary boundaries).Nodup := by
  unfold insertUniqueBoundary
  by_cases hMem : boundary ∈ boundaries
  · simpa [hMem] using h
  · simp [hMem, h]

/-- Deduplicate while preserving membership. -/
def uniqueBoundaries : List Boundary → List Boundary
  | [] => []
  | boundary :: rest =>
      insertUniqueBoundary boundary (uniqueBoundaries rest)

theorem mem_uniqueBoundaries {existing : Boundary}
    {boundaries : List Boundary} :
    existing ∈ uniqueBoundaries boundaries ↔ existing ∈ boundaries := by
  induction boundaries with
  | nil =>
      simp [uniqueBoundaries]
  | cons head tail ih =>
      simp [uniqueBoundaries, mem_insertUniqueBoundary, ih]

theorem nodup_uniqueBoundaries (boundaries : List Boundary) :
    (uniqueBoundaries boundaries).Nodup := by
  induction boundaries with
  | nil =>
      simp [uniqueBoundaries]
  | cons head tail ih =>
      exact nodup_insertUniqueBoundary ih

/-- Raw recursive detector over a reverse DFS candidate stack. -/
def detectStackRaw {g : Graph} (scan : HairpinScanAssignment g) :
    List Link → List Boundary
  | [] => []
  | left :: rest =>
      let tail := detectStackRaw scan rest
      if scan.startsStem left then
        match firstExtensionBeforeTerminal? scan rest with
        | some right => Boundary.ordered left right :: tail
        | none => tail
      else
        tail

/-- Deduplicated reference detector over a supplied reverse DFS stack. -/
def detectStack {g : Graph} (scan : HairpinScanAssignment g)
    (stack : List Link) : List Boundary :=
  uniqueBoundaries (detectStackRaw scan stack)

/-- Reference detector over the reverse real-black stack of a traversal frame. -/
def detectHairpins {g : Graph} (frame : TraversalFrame g)
    (scan : HairpinScanAssignment g) : List Boundary :=
  detectStack scan (reverseCandidateStack frame)

theorem mem_detectStackRaw {g : Graph}
    {scan : HairpinScanAssignment g} {stack : List Link}
    {boundary : Boundary} :
    boundary ∈ detectStackRaw scan stack ↔
      CanonicalScanBoundary scan stack boundary := by
  induction stack with
  | nil =>
      constructor
      · intro h
        cases h
      · intro h
        cases h
  | cons left rest ih =>
      simp only [detectStackRaw]
      by_cases hStart : scan.startsStem left = true
      · simp [hStart]
        cases hClose : firstExtensionBeforeTerminal? scan rest with
        | none =>
            simp only [ih]
            constructor
            · intro h
              exact CanonicalScanBoundary.there h
            · intro h
              cases h with
              | here hHeadStart hHeadClose =>
                  rw [hClose] at hHeadClose
                  cases hHeadClose
              | there hTail =>
                  exact hTail
        | some right =>
            simp only [ih, List.mem_cons]
            constructor
            · intro h
              cases h with
              | inl hHead =>
                  cases hHead
                  exact CanonicalScanBoundary.here hStart hClose
              | inr hTail =>
                  exact CanonicalScanBoundary.there hTail
            · intro h
              cases h with
              | here hHeadStart hHeadClose =>
                  rw [hClose] at hHeadClose
                  cases hHeadClose
                  exact Or.inl rfl
              | there hTail =>
                  exact Or.inr hTail
      · simp [hStart, ih]
        constructor
        · intro h
          exact CanonicalScanBoundary.there h
        · intro h
          cases h with
          | here hHeadStart hHeadClose =>
              exfalso
              exact hStart hHeadStart
          | there hTail =>
              exact hTail

theorem mem_detectStack {g : Graph}
    {scan : HairpinScanAssignment g} {stack : List Link}
    {boundary : Boundary} :
    boundary ∈ detectStack scan stack ↔
      CanonicalScanBoundary scan stack boundary := by
  unfold detectStack
  rw [mem_uniqueBoundaries (existing := boundary)]
  exact mem_detectStackRaw (scan := scan) (stack := stack) (boundary := boundary)

theorem detectStack_noDuplicates {g : Graph}
    (scan : HairpinScanAssignment g) (stack : List Link) :
    NoDuplicateBoundaries (detectStack scan stack) :=
  nodup_uniqueBoundaries _

theorem detectHairpins_noDuplicates {g : Graph}
    (frame : TraversalFrame g) (scan : HairpinScanAssignment g) :
    NoDuplicateBoundaries (detectHairpins frame scan) :=
  detectStack_noDuplicates scan (reverseCandidateStack frame)

end Hairpin
end Algorithms
end PovuLean
