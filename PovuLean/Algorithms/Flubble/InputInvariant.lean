import PovuLean.GFA.Basic

/-!
Input-side contracts for the verified flubble reference detector.

The current core graph scaffold exposes traversal frames but does not yet build
the exact povu DFS tree/bracket stack.  The flubble detector therefore accepts a
finite, certified class assignment over a supplied traversal frame.  The
certificate is the formal boundary between the paper's cycle-equivalence
characterization and the executable reference enumeration proved in this module
family.
-/

namespace PovuLean
namespace Algorithms
namespace Flubble

open Core

/-- A real black tree edge can serve as a flubble boundary candidate. -/
def IsBoundaryCandidate {g : Graph} (frame : TraversalFrame g) (edge : Link) : Prop :=
  edge ∈ frame.treeLinks
    ∧ edge.color = EdgeColor.black
    ∧ edge.provenance = EdgeProvenance.real

instance decidableIsBoundaryCandidate {g : Graph}
    (frame : TraversalFrame g) (edge : Link) :
    Decidable (IsBoundaryCandidate frame edge) :=
  inferInstanceAs
    (Decidable
      ( edge ∈ frame.treeLinks
        ∧ edge.color = EdgeColor.black
        ∧ edge.provenance = EdgeProvenance.real))

/--
Executable candidate filter.  Its order is the canonical stack order used by
the reference detector; downstream tree construction may impose a stronger DFS
interpretation on `frame.treeLinks`.
-/
def candidateStack {g : Graph} (frame : TraversalFrame g) : List Link :=
  frame.treeLinks.filter (fun edge => decide (IsBoundaryCandidate frame edge))

theorem mem_candidateStack {g : Graph} {frame : TraversalFrame g} {edge : Link} :
    edge ∈ candidateStack frame ↔ IsBoundaryCandidate frame edge := by
  unfold candidateStack
  rw [List.mem_filter]
  constructor
  · intro h
    exact (decide_eq_true_eq (p := IsBoundaryCandidate frame edge)).mp h.2
  · intro h
    exact ⟨h.1, (decide_eq_true_eq (p := IsBoundaryCandidate frame edge)).mpr h⟩

/-- Flubble-level input invariant over an already constructed core graph. -/
structure SupportedInput (g : Graph) (frame : TraversalFrame g) : Prop where
  graph : g.BidirectedWellFormed
  treeLinks_nodup : frame.treeLinks.Nodup

namespace SupportedInput

theorem candidateStack_nodup {g : Graph} {frame : TraversalFrame g}
    (h : SupportedInput g frame) :
    (candidateStack frame).Nodup := by
  unfold candidateStack
  exact h.treeLinks_nodup.filter _

end SupportedInput

/--
A cycle-equivalence class assignment for tree edges.  This is the pure Lean
counterpart of the equivalence classes assigned by povu's bracket pass.
-/
structure CycleClassAssignment (g : Graph) where
  classOf : Link → Nat

namespace CycleClassAssignment

def sameClass {g : Graph} (classes : CycleClassAssignment g)
    (left right : Link) : Prop :=
  classes.classOf left = classes.classOf right

instance decidableSameClass {g : Graph} (classes : CycleClassAssignment g)
    (left right : Link) :
    Decidable (classes.sameClass left right) :=
  inferInstanceAs (Decidable (classes.classOf left = classes.classOf right))

end CycleClassAssignment

/--
Correctness contract for a class assignment on the supported candidate stack.
Soundness connects equal class ids to the paper's cycle-equivalence relation;
completeness records the reverse direction for later bracket/DFS refinements.
-/
structure CycleClassAssignment.Correct {g : Graph} (frame : TraversalFrame g)
    (classes : CycleClassAssignment g) : Prop where
  sound :
    ∀ {left right : Link},
      IsBoundaryCandidate frame left →
      IsBoundaryCandidate frame right →
      classes.sameClass left right →
      CycleEquivalentEdges g left right
  complete :
    ∀ {left right : Link},
      IsBoundaryCandidate frame left →
      IsBoundaryCandidate frame right →
      CycleEquivalentEdges g left right →
      classes.sameClass left right

/--
GFA-facing specialization: accepted semantic GFA documents already construct a
bidirected core graph.  A flubble detector invocation over such a document still
needs a certified traversal frame and cycle-class assignment, because the byte
parser and C++ bracket-stack conformance harness are intentionally outside this
Lean task.
-/
structure SupportedGFAInput (doc : GFA.Document)
    (accepted : doc.Accepted)
    (frame : TraversalFrame doc.toGraph) : Prop where
  treeLinks_nodup : frame.treeLinks.Nodup

namespace SupportedGFAInput

theorem toSupportedInput {doc : GFA.Document} {accepted : doc.Accepted}
    {frame : TraversalFrame doc.toGraph}
    (h : SupportedGFAInput doc accepted frame) :
    SupportedInput doc.toGraph frame :=
  { graph := GFA.Document.accepted_toGraph_bidirectedWellFormed accepted
    treeLinks_nodup := h.treeLinks_nodup }

end SupportedGFAInput

end Flubble
end Algorithms
end PovuLean
