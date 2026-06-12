import PovuLean.Algorithms.FlubbleTree.Build

/-!
Correctness theorem for the flubble hierarchy builder.

The proof composes the base flubble detector theorem with the hierarchy
builder.  Under the named laminarity/support condition, every hierarchy node is
one discovered flubble, every discovered flubble appears exactly once, parent
links are strict nesting links, and output order/canonical boundary orientation
are preserved for downstream formats.
-/

namespace PovuLean
namespace Algorithms
namespace FlubbleTree

open Core

theorem chooseNearestParent?_some_sound {stack : List Link}
    {child parent : Boundary} {boundaries : List Boundary}
    (hParent : chooseNearestParent? stack child boundaries = some parent) :
    parent ∈ boundaries ∧ NestedIn stack parent child := by
  revert parent
  induction boundaries with
  | nil =>
      intro parent hParent
      simp [chooseNearestParent?] at hParent
  | cons candidate rest ih =>
      intro parent hParent
      unfold chooseNearestParent? at hParent
      by_cases hNested : nestedBool stack candidate child = true
      · cases hTail : chooseNearestParent? stack child rest with
        | none =>
            simp [hNested, hTail] at hParent
            cases hParent
            exact
              ⟨List.Mem.head _,
                (nestedBool_eq_true (stack := stack)
                  (parent := candidate) (child := child)).mp hNested⟩
        | some tailParent =>
            simp [hNested, hTail] at hParent
            cases hParent
            have hTailSound := ih hTail
            cases betterParent_eq_left_or_right
                (stack := stack) (left := candidate) (right := tailParent) with
            | inl hBetter =>
                constructor
                · rw [hBetter]
                  exact List.Mem.head _
                · rw [hBetter]
                  exact
                    (nestedBool_eq_true (stack := stack)
                      (parent := candidate) (child := child)).mp hNested
            | inr hBetter =>
                constructor
                · rw [hBetter]
                  exact List.Mem.tail _ hTailSound.1
                · rw [hBetter]
                  exact hTailSound.2
      · simp [hNested] at hParent
        have hTailSound := ih hParent
        exact ⟨List.Mem.tail _ hTailSound.1, hTailSound.2⟩

theorem chooseNearestParent?_none_no_parent {stack : List Link}
    {child : Boundary} {boundaries : List Boundary}
    (hNone : chooseNearestParent? stack child boundaries = none) :
    ¬ ∃ parent, parent ∈ boundaries ∧ NestedIn stack parent child := by
  induction boundaries with
  | nil =>
      intro hExists
      rcases hExists with ⟨parent, hMem, _hNested⟩
      cases hMem
  | cons candidate rest ih =>
      unfold chooseNearestParent? at hNone
      by_cases hNestedCandidate :
          nestedBool stack candidate child = true
      · cases hTail : chooseNearestParent? stack child rest <;>
          simp [hNestedCandidate, hTail] at hNone
      · simp [hNestedCandidate] at hNone
        intro hExists
        rcases hExists with ⟨parent, hMem, hNested⟩
        cases hMem with
        | head =>
            have hNestedTrue :
                nestedBool stack candidate child = true :=
              (nestedBool_eq_true (stack := stack)
                (parent := candidate) (child := child)).mpr hNested
            exact hNestedCandidate hNestedTrue
        | tail _ hMemRest =>
            exact ih hNone ⟨parent, hMemRest, hNested⟩

theorem map_nodeFor_boundaries (stack : List Link)
    (allBoundaries boundaries : List Boundary) :
    (boundaries.map (nodeFor stack allBoundaries)).map Node.boundary =
      boundaries := by
  induction boundaries with
  | nil =>
      simp
  | cons boundary rest ih =>
      simp [nodeFor, ih]

theorem buildHierarchyFrom_boundaries (stack : List Link)
    (boundaries : List Boundary) :
    (buildHierarchyFrom stack boundaries).boundaries = boundaries := by
  simp [buildHierarchyFrom, Hierarchy.boundaries,
    map_nodeFor_boundaries stack boundaries boundaries]

theorem buildHierarchyFrom_nodes_length (stack : List Link)
    (boundaries : List Boundary) :
    (buildHierarchyFrom stack boundaries).nodes.length = boundaries.length := by
  simp [buildHierarchyFrom]

theorem buildHierarchyFrom_nodes_length_le (stack : List Link)
    (boundaries : List Boundary) :
    (buildHierarchyFrom stack boundaries).nodes.length ≤ boundaries.length := by
  rw [buildHierarchyFrom_nodes_length]
  exact Nat.le_refl _

theorem buildHierarchyFrom_parentReferenceSlots_length
    (stack : List Link) (boundaries : List Boundary) :
    (buildHierarchyFrom stack boundaries).parentReferenceSlots =
      boundaries.length := by
  simp [Hierarchy.parentReferenceSlots, buildHierarchyFrom_nodes_length]

theorem buildHierarchyFrom_boundaryReferenceSlots_length
    (stack : List Link) (boundaries : List Boundary) :
    (buildHierarchyFrom stack boundaries).boundaryReferenceSlots =
      boundaries.length := by
  simp [Hierarchy.boundaryReferenceSlots, buildHierarchyFrom_nodes_length]

/--
For supported hierarchy inputs, parent and boundary-reference metadata remain
flat: exactly one optional parent slot and one boundary-reference slot are
stored per node.  No additional list-size invariant is needed because `Node`
has no list-valued parent/reference metadata field.
-/
theorem buildHierarchyFrom_flatMetadataSlots_length
    {stack : List Link} {boundaries : List Boundary}
    (_hSupported : SupportedHierarchyInput stack boundaries) :
    (buildHierarchyFrom stack boundaries).flatMetadataSlots =
      boundaries.length + boundaries.length := by
  simp [Hierarchy.flatMetadataSlots,
    buildHierarchyFrom_parentReferenceSlots_length,
    buildHierarchyFrom_boundaryReferenceSlots_length]

theorem buildHierarchyFrom_flatMetadataSlots_le_two_mul
    {stack : List Link} {boundaries : List Boundary}
    (hSupported : SupportedHierarchyInput stack boundaries) :
    (buildHierarchyFrom stack boundaries).flatMetadataSlots ≤
      boundaries.length + boundaries.length := by
  rw [buildHierarchyFrom_flatMetadataSlots_length hSupported]
  exact Nat.le_refl _

theorem buildHierarchyFrom_parent_correct {stack : List Link}
    {boundaries : List Boundary} {node : Node}
    (hNode : node ∈ (buildHierarchyFrom stack boundaries).nodes) :
    NodeParentCorrect stack boundaries node := by
  unfold buildHierarchyFrom at hNode
  rcases List.mem_map.mp hNode with ⟨boundary, hBoundaryMem, hNodeEq⟩
  cases hNodeEq
  unfold NodeParentCorrect nodeFor canonicalParent?
  constructor
  · exact hBoundaryMem
  · cases hParent : chooseNearestParent? stack boundary boundaries with
    | none =>
        exact chooseNearestParent?_none_no_parent hParent
    | some parent =>
        exact chooseNearestParent?_some_sound hParent

theorem buildHierarchyFrom_correct {stack : List Link}
    {boundaries : List Boundary}
    (hNodup : boundaries.Nodup)
    (hSupported : SupportedHierarchyInput stack boundaries) :
    IsForestHierarchy stack boundaries
      (buildHierarchyFrom stack boundaries) := by
  constructor
  · exact buildHierarchyFrom_boundaries stack boundaries
  constructor
  · rw [buildHierarchyFrom_boundaries]
    exact hNodup
  constructor
  · exact hSupported
  · intro node hNode
    exact buildHierarchyFrom_parent_correct hNode

theorem buildHierarchy_correct {g : Graph} {frame : TraversalFrame g}
    {classes : Flubble.CycleClassAssignment g}
    (hInput : Flubble.SupportedInput g frame)
    (hClasses : Flubble.CycleClassAssignment.Correct frame classes)
    (hSupported :
      SupportedHierarchyInput
        (Flubble.candidateStack frame)
        (Flubble.detectFlubbles frame classes)) :
    IsCorrectHierarchy frame classes (buildHierarchy frame classes) := by
  have hFlubbleCorrect := Flubble.detectFlubbles_correct hInput hClasses
  unfold IsCorrectHierarchy
  constructor
  · unfold buildHierarchy
    exact buildHierarchyFrom_correct hFlubbleCorrect.2.2 hSupported
  constructor
  · intro boundary
    unfold buildHierarchy
    rw [buildHierarchyFrom_boundaries]
    constructor
    · intro hMem
      exact hFlubbleCorrect.1 boundary hMem
    · intro hBoundary
      exact hFlubbleCorrect.2.1 boundary hBoundary
  · unfold OutputOrderCanonical buildHierarchy
    constructor
    · exact
        buildHierarchyFrom_boundaries
          (Flubble.candidateStack frame)
          (Flubble.detectFlubbles frame classes)
    constructor
    · rw [buildHierarchyFrom_boundaries]
      exact hFlubbleCorrect.2.2
    · intro boundary hMem
      rw [buildHierarchyFrom_boundaries] at hMem
      have hBoundary := hFlubbleCorrect.1 boundary hMem
      exact hBoundary.1.2.2.2.1

theorem buildHierarchy_nodes_length {g : Graph}
    (frame : TraversalFrame g)
    (classes : Flubble.CycleClassAssignment g) :
    (buildHierarchy frame classes).nodes.length =
      (Flubble.detectFlubbles frame classes).length := by
  unfold buildHierarchy
  exact buildHierarchyFrom_nodes_length
    (Flubble.candidateStack frame)
    (Flubble.detectFlubbles frame classes)

theorem buildHierarchy_nodes_length_le_detectFlubbles_length {g : Graph}
    (frame : TraversalFrame g)
    (classes : Flubble.CycleClassAssignment g) :
    (buildHierarchy frame classes).nodes.length ≤
      (Flubble.detectFlubbles frame classes).length := by
  rw [buildHierarchy_nodes_length]
  exact Nat.le_refl _

theorem buildHierarchy_parentReferenceSlots_length {g : Graph}
    (frame : TraversalFrame g)
    (classes : Flubble.CycleClassAssignment g) :
    (buildHierarchy frame classes).parentReferenceSlots =
      (Flubble.detectFlubbles frame classes).length := by
  unfold buildHierarchy
  exact buildHierarchyFrom_parentReferenceSlots_length
    (Flubble.candidateStack frame)
    (Flubble.detectFlubbles frame classes)

theorem buildHierarchy_boundaryReferenceSlots_length {g : Graph}
    (frame : TraversalFrame g)
    (classes : Flubble.CycleClassAssignment g) :
    (buildHierarchy frame classes).boundaryReferenceSlots =
      (Flubble.detectFlubbles frame classes).length := by
  unfold buildHierarchy
  exact buildHierarchyFrom_boundaryReferenceSlots_length
    (Flubble.candidateStack frame)
    (Flubble.detectFlubbles frame classes)

theorem buildHierarchy_flatMetadataSlots_length {g : Graph}
    {frame : TraversalFrame g}
    {classes : Flubble.CycleClassAssignment g}
    (hSupported :
      SupportedHierarchyInput
        (Flubble.candidateStack frame)
        (Flubble.detectFlubbles frame classes)) :
    (buildHierarchy frame classes).flatMetadataSlots =
      (Flubble.detectFlubbles frame classes).length +
        (Flubble.detectFlubbles frame classes).length := by
  unfold buildHierarchy
  exact buildHierarchyFrom_flatMetadataSlots_length hSupported

theorem buildHierarchy_nodes_length_le_graph_edgeCount {g : Graph}
    (frame : TraversalFrame g)
    (classes : Flubble.CycleClassAssignment g)
    (hTreeLinksNodup : frame.treeLinks.Nodup) :
    (buildHierarchy frame classes).nodes.length ≤ g.edgeCount := by
  exact Nat.le_trans
    (buildHierarchy_nodes_length_le_detectFlubbles_length frame classes)
    (Flubble.detectFlubbles_length_le_graph_edgeCount
      frame classes hTreeLinksNodup)

theorem buildHierarchy_nodes_length_le_supportedInput_graph_edgeCount
    {g : Graph} {frame : TraversalFrame g}
    (classes : Flubble.CycleClassAssignment g)
    (hInput : Flubble.SupportedInput g frame) :
    (buildHierarchy frame classes).nodes.length ≤ g.edgeCount :=
  buildHierarchy_nodes_length_le_graph_edgeCount
    frame classes hInput.treeLinks_nodup

theorem buildHierarchy_flatMetadataSlots_le_graph_edgeCount_twice
    {g : Graph} {frame : TraversalFrame g}
    {classes : Flubble.CycleClassAssignment g}
    (hTreeLinksNodup : frame.treeLinks.Nodup)
    (hSupported :
      SupportedHierarchyInput
        (Flubble.candidateStack frame)
        (Flubble.detectFlubbles frame classes)) :
    (buildHierarchy frame classes).flatMetadataSlots ≤
      g.edgeCount + g.edgeCount := by
  rw [buildHierarchy_flatMetadataSlots_length hSupported]
  exact Nat.add_le_add
    (Flubble.detectFlubbles_length_le_graph_edgeCount
      frame classes hTreeLinksNodup)
    (Flubble.detectFlubbles_length_le_graph_edgeCount
      frame classes hTreeLinksNodup)

theorem buildHierarchy_flatMetadataSlots_le_supportedInput_graph_edgeCount_twice
    {g : Graph} {frame : TraversalFrame g}
    {classes : Flubble.CycleClassAssignment g}
    (hInput : Flubble.SupportedInput g frame)
    (hSupported :
      SupportedHierarchyInput
        (Flubble.candidateStack frame)
        (Flubble.detectFlubbles frame classes)) :
    (buildHierarchy frame classes).flatMetadataSlots ≤
      g.edgeCount + g.edgeCount :=
  buildHierarchy_flatMetadataSlots_le_graph_edgeCount_twice
    hInput.treeLinks_nodup hSupported

theorem buildHierarchy_rejects_unsupported {stack : List Link}
    {boundaries : List Boundary}
    (hSupported : SupportedHierarchyInput stack boundaries) :
    ¬ UnsupportedHierarchyInput stack boundaries :=
  supportedHierarchyInput_not_unsupported hSupported

end FlubbleTree
end Algorithms
end PovuLean
