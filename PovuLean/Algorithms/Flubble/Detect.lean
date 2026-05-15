import PovuLean.Algorithms.Flubble.Spec

/-!
Executable Lean reference detector for base flubble boundaries.

The detector is intentionally small: it traverses the canonical candidate stack,
emits the next same-class edge after a nonempty gap, normalizes orientation by
edge id, and removes exact duplicates.  The correctness module proves that this
list is exactly the canonical flubble-boundary specification under the supplied
class-assignment certificate.
-/

namespace PovuLean
namespace Algorithms
namespace Flubble

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

/-- Raw recursive detector over a candidate stack. -/
def detectStackRaw {g : Graph} (classes : CycleClassAssignment g) :
    List Link → List Boundary
  | [] => []
  | left :: rest =>
      let tail := detectStackRaw classes rest
      match closeAfterGap? classes left rest with
      | some right => Boundary.ordered left right :: tail
      | none => tail

/-- Deduplicated, deterministic reference detector over a candidate stack. -/
def detectStack {g : Graph} (classes : CycleClassAssignment g)
    (stack : List Link) : List Boundary :=
  uniqueBoundaries (detectStackRaw classes stack)

/-- Reference detector over the canonical real-black stack of a frame. -/
def detectFlubbles {g : Graph} (frame : TraversalFrame g)
    (classes : CycleClassAssignment g) : List Boundary :=
  detectStack classes (candidateStack frame)

theorem mem_detectStackRaw {g : Graph}
    {classes : CycleClassAssignment g} {stack : List Link}
    {boundary : Boundary} :
    boundary ∈ detectStackRaw classes stack ↔
      CanonicalClassBoundary classes stack boundary := by
  induction stack with
  | nil =>
      constructor
      · intro h
        cases h
      · intro h
        cases h
  | cons left rest ih =>
      simp only [detectStackRaw]
      cases hClose : closeAfterGap? classes left rest with
      | none =>
          simp only [ih]
          constructor
          · intro h
            exact CanonicalClassBoundary.there h
          · intro h
            cases h with
            | here hHead =>
                rw [hClose] at hHead
                cases hHead
            | there hTail =>
                exact hTail
      | some right =>
          simp only [ih, List.mem_cons]
          constructor
          · intro h
            cases h with
            | inl hHead =>
                cases hHead
                exact CanonicalClassBoundary.here hClose
            | inr hTail =>
                exact CanonicalClassBoundary.there hTail
          · intro h
            cases h with
            | here hHead =>
                rw [hClose] at hHead
                cases hHead
                exact Or.inl rfl
            | there hTail =>
                exact Or.inr hTail

theorem mem_detectStack {g : Graph}
    {classes : CycleClassAssignment g} {stack : List Link}
    {boundary : Boundary} :
    boundary ∈ detectStack classes stack ↔
      CanonicalClassBoundary classes stack boundary := by
  unfold detectStack
  rw [mem_uniqueBoundaries (existing := boundary)]
  exact mem_detectStackRaw (classes := classes) (stack := stack) (boundary := boundary)

theorem detectStack_noDuplicates {g : Graph}
    (classes : CycleClassAssignment g) (stack : List Link) :
    NoDuplicateBoundaries (detectStack classes stack) :=
  nodup_uniqueBoundaries _

theorem detectFlubbles_noDuplicates {g : Graph}
    (frame : TraversalFrame g) (classes : CycleClassAssignment g) :
    NoDuplicateBoundaries (detectFlubbles frame classes) :=
  detectStack_noDuplicates classes (candidateStack frame)

end Flubble
end Algorithms
end PovuLean
