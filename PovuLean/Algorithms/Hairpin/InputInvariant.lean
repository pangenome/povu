import PovuLean.Algorithms.Flubble.InputInvariant

/-!
Input-side contracts for the verified hairpin reference detector.

Hairpin detection in povu is performed during the same reverse DFS/bracket pass
that assigns flubble cycle-equivalence classes.  This module deliberately reuses
the flubble task's supported graph, GFA, and real-black tree-edge candidate
interfaces instead of redefining them here.

The current Lean scaffold does not yet construct the concrete C++ bracket list,
so the hairpin detector accepts a finite scan assignment over the supplied
reverse candidate stack.  The correctness module proves the executable scan
against that certified assignment.
-/

namespace PovuLean
namespace Algorithms
namespace Hairpin

open Core

/-- Hairpin boundaries use the same real black tree-edge candidate predicate as flubbles. -/
abbrev IsBoundaryCandidate {g : Graph} (frame : TraversalFrame g) (edge : Link) : Prop :=
  Flubble.IsBoundaryCandidate frame edge

/-- Canonical DFS stack inherited from the flubble intermediate interface. -/
abbrev candidateStack {g : Graph} (frame : TraversalFrame g) : List Link :=
  Flubble.candidateStack frame

/-- Supported core input for hairpin detection, reused from the flubble proof. -/
abbrev SupportedInput (g : Graph) (frame : TraversalFrame g) : Prop :=
  Flubble.SupportedInput g frame

/-- GFA-facing supported input, reused from the flubble proof. -/
abbrev SupportedGFAInput (doc : GFA.Document)
    (accepted : doc.Accepted)
    (frame : TraversalFrame doc.toGraph) : Prop :=
  Flubble.SupportedGFAInput doc accepted frame

/--
Hairpins are found while retreating through the DFS tree, so the reference scan
uses the reverse of the canonical candidate stack.
-/
def reverseCandidateStack {g : Graph} (frame : TraversalFrame g) : List Link :=
  (candidateStack frame).reverse

theorem mem_reverseCandidateStack {g : Graph} {frame : TraversalFrame g}
    {edge : Link} :
    edge ∈ reverseCandidateStack frame ↔ IsBoundaryCandidate frame edge := by
  simp [reverseCandidateStack, candidateStack, IsBoundaryCandidate,
    Flubble.mem_candidateStack]

theorem reverseCandidateStack_nodup {g : Graph} {frame : TraversalFrame g}
    (h : SupportedInput g frame) :
    (reverseCandidateStack frame).Nodup := by
  exact
    (List.Perm.nodup_iff
      (List.reverse_perm (candidateStack frame))).mpr h.candidateStack_nodup

/--
Certified per-edge events consumed by the Lean reference scan.

* `startsStem` corresponds to the bracket list becoming empty on a real
  non-dummy tree vertex.
* `extendsStem` corresponds to the simplifying backedge remaining on top while
  inside a hairpin.
* `terminatesStem` cuts off a candidate if a leaf/root boundary is reached
  before an extension is seen.
-/
structure HairpinScanAssignment (g : Graph) where
  startsStem : Link → Bool
  extendsStem : Link → Bool
  terminatesStem : Link → Bool

namespace HairpinScanAssignment

def StartsStem {g : Graph} (scan : HairpinScanAssignment g) (edge : Link) : Prop :=
  scan.startsStem edge = true

def ExtendsStem {g : Graph} (scan : HairpinScanAssignment g) (edge : Link) : Prop :=
  scan.extendsStem edge = true

def TerminatesStem {g : Graph} (scan : HairpinScanAssignment g) (edge : Link) : Prop :=
  scan.terminatesStem edge = true

instance decidableStartsStem {g : Graph}
    (scan : HairpinScanAssignment g) (edge : Link) :
    Decidable (scan.StartsStem edge) :=
  inferInstanceAs (Decidable (scan.startsStem edge = true))

instance decidableExtendsStem {g : Graph}
    (scan : HairpinScanAssignment g) (edge : Link) :
    Decidable (scan.ExtendsStem edge) :=
  inferInstanceAs (Decidable (scan.extendsStem edge = true))

instance decidableTerminatesStem {g : Graph}
    (scan : HairpinScanAssignment g) (edge : Link) :
    Decidable (scan.TerminatesStem edge) :=
  inferInstanceAs (Decidable (scan.terminatesStem edge = true))

end HairpinScanAssignment

end Hairpin
end Algorithms
end PovuLean
