import PovuLean.Algorithms.Flubble.Detect

/-!
Indexed one-pass flubble detector interface.

This module keeps the canonical `detectFlubbles` reference detector intact and
adds a detector shaped like the intended linear implementation.  The scan walks
the candidate stack from right to left, carrying an explicit class-index state
whose lookup returns the nearest later occurrence of a cycle class.  Each
candidate performs one class lookup and one state update; replacing the concrete
association-list state by a constant-time finite map is isolated behind the
same lookup/update behavior proved below.
-/

namespace PovuLean
namespace Algorithms
namespace Flubble

open Core

/-- A candidate occurrence annotated with its stack index. -/
structure IndexedOccurrence where
  index : Nat
  edge : Link
  deriving Repr, DecidableEq

/--
Concrete class-index state used for the proof.

Entries are ordered by stack position.  While scanning from right to left, the
state for the already-processed suffix stores the nearest later occurrence of a
class before any farther occurrence of that same class.
-/
structure ClassIndexState where
  entries : List (Nat × IndexedOccurrence)
  deriving Repr, DecidableEq

namespace ClassIndexState

/-- Empty class-index state before any suffix has been processed. -/
def empty : ClassIndexState :=
  { entries := [] }

/-- Lookup the nearest indexed occurrence currently stored for a class id. -/
def lookup (classId : Nat) : ClassIndexState → Option IndexedOccurrence
  | { entries := [] } => none
  | { entries := (entryClass, occurrence) :: rest } =>
      if entryClass = classId then
        some occurrence
      else
        lookup classId { entries := rest }

/-- Update the index so this occurrence shadows farther occurrences. -/
def update (classId : Nat) (occurrence : IndexedOccurrence)
    (state : ClassIndexState) : ClassIndexState :=
  { entries := (classId, occurrence) :: state.entries }

end ClassIndexState

/-- Entries contributed by a stack suffix starting at `startIndex`. -/
def indexedEntriesFrom {g : Graph} (classes : CycleClassAssignment g)
    (startIndex : Nat) : List Link → List (Nat × IndexedOccurrence)
  | [] => []
  | edge :: rest =>
      (classes.classOf edge, { index := startIndex, edge := edge }) ::
        indexedEntriesFrom classes (startIndex + 1) rest

/-- Result of the indexed detector scan, including explicit operation counts. -/
structure IndexedScanResult where
  state : ClassIndexState
  boundaries : List Boundary
  visited : Nat
  lookups : Nat
  updates : Nat
  deriving Repr, DecidableEq

namespace IndexedScanResult

/-- A scan result with no processed candidates. -/
def emptyWithState (state : ClassIndexState) : IndexedScanResult :=
  { state := state, boundaries := [], visited := 0, lookups := 0, updates := 0 }

end IndexedScanResult

/--
One indexed detector step.

The caller supplies the state for the already-processed suffix.  The step uses
one lookup for `left`'s class and suppresses immediately adjacent same-class
occurrences, matching `closeAfterGap?`.
-/
def indexedBoundaryAt? {g : Graph} (classes : CycleClassAssignment g)
    (left : Link) (rest : List Link) (state : ClassIndexState) :
    Option Boundary :=
  match rest with
  | [] => none
  | immediate :: _ =>
      if classes.sameClassBool left immediate then
        none
      else
        match state.lookup (classes.classOf left) with
        | some occurrence => some (Boundary.ordered left occurrence.edge)
        | none => none

/--
Indexed suffix scan.  It processes `stack` from right to left; `startIndex`
is the absolute index of the first element of `stack`.
-/
def scanStackIndexedAux {g : Graph} (classes : CycleClassAssignment g)
    (startIndex : Nat) : List Link → ClassIndexState → IndexedScanResult
  | [], state => IndexedScanResult.emptyWithState state
  | left :: rest, state =>
      let suffix := scanStackIndexedAux classes (startIndex + 1) rest state
      let classId := classes.classOf left
      let occurrence : IndexedOccurrence := { index := startIndex, edge := left }
      let state' := suffix.state.update classId occurrence
      let emitted :=
        match indexedBoundaryAt? classes left rest suffix.state with
        | some boundary => boundary :: suffix.boundaries
        | none => suffix.boundaries
      { state := state'
        boundaries := emitted
        visited := suffix.visited + 1
        lookups := suffix.lookups + 1
        updates := suffix.updates + 1 }

/-- Indexed raw detector over a stack. -/
def detectStackIndexedRaw {g : Graph} (classes : CycleClassAssignment g)
    (stack : List Link) : List Boundary :=
  (scanStackIndexedAux classes 0 stack ClassIndexState.empty).boundaries

/-- Deduplicated indexed detector over a stack. -/
def detectStackIndexed {g : Graph} (classes : CycleClassAssignment g)
    (stack : List Link) : List Boundary :=
  uniqueBoundaries (detectStackIndexedRaw classes stack)

/-- Indexed detector over the canonical real-black stack of a frame. -/
def detectFlubblesIndexed {g : Graph} (frame : TraversalFrame g)
    (classes : CycleClassAssignment g) : List Boundary :=
  detectStackIndexed classes (candidateStack frame)

theorem scanStackIndexedAux_entries {g : Graph}
    (classes : CycleClassAssignment g) (stack : List Link)
    (startIndex : Nat) (state : ClassIndexState) :
    (scanStackIndexedAux classes startIndex stack state).state.entries =
      indexedEntriesFrom classes startIndex stack ++ state.entries := by
  induction stack generalizing startIndex state with
  | nil =>
      simp [scanStackIndexedAux, IndexedScanResult.emptyWithState,
        indexedEntriesFrom]
  | cons left rest ih =>
      simp [scanStackIndexedAux, indexedEntriesFrom, ClassIndexState.update,
        ih]

theorem lookup_indexedEntriesFrom_eq_firstSameClass? {g : Graph}
    (classes : CycleClassAssignment g) (left : Link)
    (stack : List Link) (startIndex : Nat) :
    Option.map IndexedOccurrence.edge
        (ClassIndexState.lookup (classes.classOf left)
          { entries := indexedEntriesFrom classes startIndex stack }) =
      firstSameClass? classes left stack := by
  induction stack generalizing startIndex with
  | nil =>
      simp [indexedEntriesFrom, ClassIndexState.lookup, firstSameClass?]
  | cons candidate rest ih =>
      by_cases hSame : classes.sameClassBool left candidate = true
      · have hClass : classes.classOf candidate = classes.classOf left := by
          exact (CycleClassAssignment.sameClassBool_eq_true.mp hSame).symm
        simp [indexedEntriesFrom, ClassIndexState.lookup, firstSameClass?,
          hSame, hClass]
      · have hClass : classes.classOf candidate ≠ classes.classOf left := by
          intro hEq
          have hSameClass : classes.sameClass left candidate := hEq.symm
          exact hSame
            (CycleClassAssignment.sameClassBool_eq_true.mpr hSameClass)
        simp [indexedEntriesFrom, ClassIndexState.lookup, firstSameClass?,
          hSame, hClass, ih]

theorem lookup_indexedEntriesFrom_boundary_eq_firstSameClass? {g : Graph}
    (classes : CycleClassAssignment g) (left : Link)
    (stack : List Link) (startIndex : Nat) :
    Option.map (fun occurrence => Boundary.ordered left occurrence.edge)
        (ClassIndexState.lookup (classes.classOf left)
          { entries := indexedEntriesFrom classes startIndex stack }) =
      Option.map (Boundary.ordered left) (firstSameClass? classes left stack) := by
  have hLookup :=
    lookup_indexedEntriesFrom_eq_firstSameClass? classes left stack startIndex
  cases h :
      ClassIndexState.lookup (classes.classOf left)
        { entries := indexedEntriesFrom classes startIndex stack } <;>
    simp [h] at hLookup ⊢
  · rw [← hLookup]
    rfl
  · rw [← hLookup]
    rfl

theorem indexedBoundaryAt?_eq_closeAfterGap? {g : Graph}
    (classes : CycleClassAssignment g) (left : Link) (rest : List Link)
    (startIndex : Nat) :
    indexedBoundaryAt? classes left rest
        { entries := indexedEntriesFrom classes startIndex rest } =
      Option.map (Boundary.ordered left) (closeAfterGap? classes left rest) := by
  cases rest with
  | nil =>
      simp [indexedBoundaryAt?, closeAfterGap?]
  | cons immediate tail =>
      by_cases hSame : classes.sameClassBool left immediate = true
      · simp [indexedBoundaryAt?, closeAfterGap?, hSame]
      · have hLookup :=
          lookup_indexedEntriesFrom_boundary_eq_firstSameClass?
            classes left (immediate :: tail) startIndex
        cases hOcc :
            ClassIndexState.lookup (classes.classOf left)
              { entries := indexedEntriesFrom classes startIndex
                  (immediate :: tail) } with
        | none =>
            simp [hOcc, firstSameClass?, hSame] at hLookup
            unfold indexedBoundaryAt? closeAfterGap?
            simp [hSame, hOcc]
            exact hLookup
        | some occurrence =>
            simp [hOcc, firstSameClass?, hSame] at hLookup
            unfold indexedBoundaryAt? closeAfterGap?
            simp [hSame, hOcc]
            exact hLookup

theorem scanStackIndexedAux_boundaries_eq_detectStackRaw {g : Graph}
    (classes : CycleClassAssignment g) (stack : List Link)
    (startIndex : Nat) :
    (scanStackIndexedAux classes startIndex stack ClassIndexState.empty).boundaries =
      detectStackRaw classes stack := by
  induction stack generalizing startIndex with
  | nil =>
      simp [scanStackIndexedAux, IndexedScanResult.emptyWithState,
        detectStackRaw]
  | cons left rest ih =>
      have hEntries :
          (scanStackIndexedAux classes (startIndex + 1) rest
            ClassIndexState.empty).state =
            { entries := indexedEntriesFrom classes (startIndex + 1) rest } := by
        have hEntries :=
          scanStackIndexedAux_entries classes rest (startIndex + 1)
            ClassIndexState.empty
        cases hState :
            (scanStackIndexedAux classes (startIndex + 1) rest
              ClassIndexState.empty).state with
        | mk entries =>
            rw [hState] at hEntries
            simp [ClassIndexState.empty] at hEntries
            cases hEntries
            rfl
      simp [scanStackIndexedAux, detectStackRaw, ih, hEntries,
        indexedBoundaryAt?_eq_closeAfterGap?]
      cases closeAfterGap? classes left rest <;> simp

theorem detectStackIndexedRaw_eq_detectStackRaw {g : Graph}
    (classes : CycleClassAssignment g) (stack : List Link) :
    detectStackIndexedRaw classes stack = detectStackRaw classes stack := by
  exact scanStackIndexedAux_boundaries_eq_detectStackRaw classes stack 0

theorem detectStackIndexed_eq_detectStack {g : Graph}
    (classes : CycleClassAssignment g) (stack : List Link) :
    detectStackIndexed classes stack = detectStack classes stack := by
  simp [detectStackIndexed, detectStack, detectStackIndexedRaw_eq_detectStackRaw]

theorem detectFlubblesIndexed_eq_detectFlubbles {g : Graph}
    (frame : TraversalFrame g) (classes : CycleClassAssignment g) :
    detectFlubblesIndexed frame classes = detectFlubbles frame classes := by
  simp [detectFlubblesIndexed, detectFlubbles, detectStackIndexed_eq_detectStack]

theorem mem_detectStackIndexed {g : Graph}
    {classes : CycleClassAssignment g} {stack : List Link}
    {boundary : Boundary} :
    boundary ∈ detectStackIndexed classes stack ↔
      CanonicalClassBoundary classes stack boundary := by
  rw [detectStackIndexed_eq_detectStack]
  exact mem_detectStack (classes := classes) (stack := stack)

theorem mem_detectFlubblesIndexed {g : Graph}
    {frame : TraversalFrame g} {classes : CycleClassAssignment g}
    {boundary : Boundary} :
    boundary ∈ detectFlubblesIndexed frame classes ↔
      CanonicalClassBoundary classes (candidateStack frame) boundary := by
  unfold detectFlubblesIndexed
  exact mem_detectStackIndexed (classes := classes)

theorem detectFlubblesIndexed_iff_isFlubbleBoundary {g : Graph}
    {frame : TraversalFrame g} {classes : CycleClassAssignment g}
    {boundary : Boundary} (hPaper : IsPaperBoundary frame boundary) :
    boundary ∈ detectFlubblesIndexed frame classes ↔
      IsFlubbleBoundary frame classes boundary := by
  rw [mem_detectFlubblesIndexed]
  constructor
  · intro hCanonical
    exact ⟨hPaper, hCanonical⟩
  · intro hBoundary
    exact hBoundary.2

theorem detectFlubblesIndexed_noDuplicates {g : Graph}
    (frame : TraversalFrame g) (classes : CycleClassAssignment g) :
    NoDuplicateBoundaries (detectFlubblesIndexed frame classes) := by
  rw [detectFlubblesIndexed_eq_detectFlubbles]
  exact detectFlubbles_noDuplicates frame classes

theorem detectFlubblesIndexed_length_le_candidateStack_length {g : Graph}
    (frame : TraversalFrame g) (classes : CycleClassAssignment g) :
    (detectFlubblesIndexed frame classes).length ≤ (candidateStack frame).length := by
  rw [detectFlubblesIndexed_eq_detectFlubbles]
  exact detectFlubbles_length_le_candidateStack_length frame classes

theorem detectFlubblesIndexed_length_le_treeLinks_length {g : Graph}
    (frame : TraversalFrame g) (classes : CycleClassAssignment g) :
    (detectFlubblesIndexed frame classes).length ≤ frame.treeLinks.length := by
  rw [detectFlubblesIndexed_eq_detectFlubbles]
  exact detectFlubbles_length_le_treeLinks_length frame classes

theorem detectFlubblesIndexed_length_le_graph_links_length {g : Graph}
    (frame : TraversalFrame g) (classes : CycleClassAssignment g)
    (hTreeLinksNodup : frame.treeLinks.Nodup) :
    (detectFlubblesIndexed frame classes).length ≤ g.links.length := by
  rw [detectFlubblesIndexed_eq_detectFlubbles]
  exact detectFlubbles_length_le_graph_links_length frame classes hTreeLinksNodup

theorem detectFlubblesIndexed_length_le_graph_edgeCount {g : Graph}
    (frame : TraversalFrame g) (classes : CycleClassAssignment g)
    (hTreeLinksNodup : frame.treeLinks.Nodup) :
    (detectFlubblesIndexed frame classes).length ≤ g.edgeCount := by
  rw [detectFlubblesIndexed_eq_detectFlubbles]
  exact detectFlubbles_length_le_graph_edgeCount frame classes hTreeLinksNodup

theorem detectFlubblesIndexed_length_le_supportedInput_graph_edgeCount {g : Graph}
    {frame : TraversalFrame g} (classes : CycleClassAssignment g)
    (hInput : SupportedInput g frame) :
    (detectFlubblesIndexed frame classes).length ≤ g.edgeCount := by
  exact detectFlubblesIndexed_length_le_graph_edgeCount
    frame classes hInput.treeLinks_nodup

theorem scanStackIndexedAux_visited {g : Graph}
    (classes : CycleClassAssignment g) (stack : List Link)
    (startIndex : Nat) (state : ClassIndexState) :
    (scanStackIndexedAux classes startIndex stack state).visited = stack.length := by
  induction stack generalizing startIndex state with
  | nil =>
      simp [scanStackIndexedAux, IndexedScanResult.emptyWithState]
  | cons left rest ih =>
      simp [scanStackIndexedAux, ih]

theorem scanStackIndexedAux_lookups {g : Graph}
    (classes : CycleClassAssignment g) (stack : List Link)
    (startIndex : Nat) (state : ClassIndexState) :
    (scanStackIndexedAux classes startIndex stack state).lookups = stack.length := by
  induction stack generalizing startIndex state with
  | nil =>
      simp [scanStackIndexedAux, IndexedScanResult.emptyWithState]
  | cons left rest ih =>
      simp [scanStackIndexedAux, ih]

theorem scanStackIndexedAux_updates {g : Graph}
    (classes : CycleClassAssignment g) (stack : List Link)
    (startIndex : Nat) (state : ClassIndexState) :
    (scanStackIndexedAux classes startIndex stack state).updates = stack.length := by
  induction stack generalizing startIndex state with
  | nil =>
      simp [scanStackIndexedAux, IndexedScanResult.emptyWithState]
  | cons left rest ih =>
      simp [scanStackIndexedAux, ih]

/--
Per-candidate shape of the indexed scan: over a stack of length `n`, the
detector performs exactly `n` state lookups and `n` state updates.
-/
theorem detectStackIndexedRaw_one_lookup_and_update_per_candidate {g : Graph}
    (classes : CycleClassAssignment g) (stack : List Link) :
    let result := scanStackIndexedAux classes 0 stack ClassIndexState.empty
    result.visited = stack.length
      ∧ result.lookups = stack.length
      ∧ result.updates = stack.length := by
  simp [scanStackIndexedAux_visited, scanStackIndexedAux_lookups,
    scanStackIndexedAux_updates]

end Flubble
end Algorithms
end PovuLean
