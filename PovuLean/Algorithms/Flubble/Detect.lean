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

theorem insertUniqueBoundary_length_le_succ (boundary : Boundary)
    (boundaries : List Boundary) :
    (insertUniqueBoundary boundary boundaries).length ≤ boundaries.length + 1 := by
  unfold insertUniqueBoundary
  by_cases hMem : boundary ∈ boundaries
  · simp [hMem]
  · simp [hMem]

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

theorem uniqueBoundaries_length_le (boundaries : List Boundary) :
    (uniqueBoundaries boundaries).length ≤ boundaries.length := by
  induction boundaries with
  | nil =>
      simp [uniqueBoundaries]
  | cons head tail ih =>
      exact Nat.le_trans
        (insertUniqueBoundary_length_le_succ head (uniqueBoundaries tail))
        (Nat.succ_le_succ ih)

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

theorem detectStackRaw_length_le_stack_length {g : Graph}
    (classes : CycleClassAssignment g) (stack : List Link) :
    (detectStackRaw classes stack).length ≤ stack.length := by
  induction stack with
  | nil =>
      simp [detectStackRaw]
  | cons left rest ih =>
      simp only [detectStackRaw]
      cases closeAfterGap? classes left rest with
      | none =>
          exact Nat.le_trans ih (Nat.le_succ rest.length)
      | some right =>
          exact Nat.succ_le_succ ih

theorem detectStack_length_le_stack_length {g : Graph}
    (classes : CycleClassAssignment g) (stack : List Link) :
    (detectStack classes stack).length ≤ stack.length := by
  exact Nat.le_trans
    (uniqueBoundaries_length_le (detectStackRaw classes stack))
    (detectStackRaw_length_le_stack_length classes stack)

theorem detectFlubbles_length_le_candidateStack_length {g : Graph}
    (frame : TraversalFrame g) (classes : CycleClassAssignment g) :
    (detectFlubbles frame classes).length ≤ (candidateStack frame).length := by
  exact detectStack_length_le_stack_length classes (candidateStack frame)

theorem candidateStack_length_le_treeLinks_length {g : Graph}
    (frame : TraversalFrame g) :
    (candidateStack frame).length ≤ frame.treeLinks.length := by
  unfold candidateStack
  exact List.length_filter_le _ _

theorem detectFlubbles_length_le_treeLinks_length {g : Graph}
    (frame : TraversalFrame g) (classes : CycleClassAssignment g) :
    (detectFlubbles frame classes).length ≤ frame.treeLinks.length := by
  exact Nat.le_trans
    (detectFlubbles_length_le_candidateStack_length frame classes)
    (candidateStack_length_le_treeLinks_length frame)

theorem nodup_length_le_of_subset {α : Type} [BEq α] [LawfulBEq α]
    {xs ys : List α} (hNodup : xs.Nodup) (hSubset : xs ⊆ ys) :
    xs.length ≤ ys.length := by
  induction xs generalizing ys with
  | nil =>
      simp
  | cons x rest ih =>
      have hRestNodup : rest.Nodup := (List.nodup_cons.mp hNodup).2
      have hNotMem : ¬x ∈ rest := (List.nodup_cons.mp hNodup).1
      have hMemYs : x ∈ ys := hSubset (List.Mem.head rest)
      have hRestSubsetErase : rest ⊆ ys.erase x := by
        intro y hy
        have hyYs : y ∈ ys := hSubset (List.Mem.tail x hy)
        have hyNe : y ≠ x := by
          intro hEq
          exact hNotMem (hEq ▸ hy)
        exact (List.mem_erase_of_ne hyNe).2 hyYs
      have hRestLength : rest.length ≤ (ys.erase x).length :=
        ih hRestNodup hRestSubsetErase
      have hEraseLength : (ys.erase x).length = ys.length - 1 :=
        List.length_erase_of_mem hMemYs
      have hEraseSucc : (ys.erase x).length + 1 = ys.length := by
        rw [hEraseLength]
        exact Nat.succ_pred_eq_of_pos (List.length_pos_of_mem hMemYs)
      simpa [Nat.succ_eq_add_one, hEraseSucc] using
        Nat.succ_le_succ hRestLength

theorem treeLinks_length_le_graph_links_length {g : Graph}
    (frame : TraversalFrame g) (hTreeLinksNodup : frame.treeLinks.Nodup) :
    frame.treeLinks.length ≤ g.links.length := by
  exact nodup_length_le_of_subset hTreeLinksNodup (by
    intro edge hEdge
    exact List.contains_iff_mem.mp (frame.treeLinks_supported edge hEdge))

theorem detectFlubbles_length_le_graph_links_length {g : Graph}
    (frame : TraversalFrame g) (classes : CycleClassAssignment g)
    (hTreeLinksNodup : frame.treeLinks.Nodup) :
    (detectFlubbles frame classes).length ≤ g.links.length := by
  exact Nat.le_trans
    (detectFlubbles_length_le_treeLinks_length frame classes)
    (treeLinks_length_le_graph_links_length frame hTreeLinksNodup)

theorem detectFlubbles_length_le_graph_edgeCount {g : Graph}
    (frame : TraversalFrame g) (classes : CycleClassAssignment g)
    (hTreeLinksNodup : frame.treeLinks.Nodup) :
    (detectFlubbles frame classes).length ≤ g.edgeCount := by
  unfold Graph.edgeCount
  exact detectFlubbles_length_le_graph_links_length frame classes hTreeLinksNodup

theorem detectFlubbles_length_le_supportedInput_graph_edgeCount {g : Graph}
    {frame : TraversalFrame g} (classes : CycleClassAssignment g)
    (hInput : SupportedInput g frame) :
    (detectFlubbles frame classes).length ≤ g.edgeCount :=
  detectFlubbles_length_le_graph_edgeCount frame classes hInput.treeLinks_nodup

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
