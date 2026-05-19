import PovuLean.Core.Basic

/-!
Walks, paths, reachability, cycles, and traversal-frame hooks over the core
graph model.
-/

namespace PovuLean
namespace Core

/--
A supported walk from one oriented segment to another, represented by the exact
oriented links traversed.  The constructors encode link continuity.
-/
inductive IsWalk (g : Graph) :
    OrientedSegment → OrientedSegment → List Link → Prop where
  | nil {v : OrientedSegment}
      (hv : g.HasOriented v) :
      IsWalk g v v []
  | cons {a b c : OrientedSegment} {rest : List Link}
      (edge : Link)
      (hSource : edge.source = a)
      (hTarget : edge.target = b)
      (ha : g.HasOriented a)
      (he : g.HasLink edge)
      (tail : IsWalk g b c rest) :
      IsWalk g a c (edge :: rest)

namespace IsWalk

theorem append {g : Graph} {a b c : OrientedSegment}
    {left right : List Link}
    (first : IsWalk g a b left)
    (second : IsWalk g b c right) :
    IsWalk g a c (left ++ right) := by
  induction first generalizing c right with
  | nil hv =>
      simpa using second
  | cons edge hSource hTarget ha he tail ih =>
      simp
      exact IsWalk.cons edge hSource hTarget ha he (ih second)

theorem snoc {g : Graph} {a b c : OrientedSegment} {edges : List Link}
    (walk : IsWalk g a b edges)
    (edge : Link)
    (hSource : edge.source = b)
    (hTarget : edge.target = c)
    (hc : g.HasOriented c)
    (he : g.HasLink edge) :
    IsWalk g a c (edges ++ [edge]) := by
  induction walk generalizing c edge with
  | nil hv =>
      simpa using IsWalk.cons edge hSource hTarget hv he (IsWalk.nil hc)
  | cons head hHeadSource hHeadTarget ha hHead tail ih =>
      simp
      exact IsWalk.cons head hHeadSource hHeadTarget ha hHead
        (ih edge hSource hTarget hc he)

/-- Edge-list reversal with per-link reverse-complement traversal. -/
def reverseEdges : List Link → List Link
  | [] => []
  | edge :: rest => reverseEdges rest ++ [edge.reverse]

@[simp] theorem reverseEdges_nil : reverseEdges [] = [] := rfl

@[simp] theorem reverseEdges_cons (edge : Link) (rest : List Link) :
    reverseEdges (edge :: rest) = reverseEdges rest ++ [edge.reverse] := rfl

@[simp] theorem reverseEdges_length (edges : List Link) :
    (reverseEdges edges).length = edges.length := by
  induction edges with
  | nil => rfl
  | cons edge rest ih =>
      simp [reverseEdges, ih]

theorem reverse {g : Graph} {a b : OrientedSegment} {edges : List Link}
    (hSym : g.EdgeSymmetric)
    (walk : IsWalk g a b edges) :
    IsWalk g b.reverse a.reverse (reverseEdges edges) := by
  induction walk with
  | nil hv =>
      exact IsWalk.nil ((Graph.hasOriented_reverse).2 hv)
  | cons edge hSource hTarget ha he tail ih =>
      have reverseSource := by
        simpa [Link.reverse] using congrArg OrientedSegment.reverse hTarget
      have reverseTarget := by
        simpa [Link.reverse] using congrArg OrientedSegment.reverse hSource
      exact IsWalk.snoc ih edge.reverse reverseSource reverseTarget
        ((Graph.hasOriented_reverse).2 ha) (hSym _ he)

end IsWalk

/-- A concrete supported traversal. -/
structure Walk (g : Graph) where
  start : OrientedSegment
  finish : OrientedSegment
  edges : List Link
  valid : IsWalk g start finish edges

/--
Pangenome paths are supported walks.  `SimplePath` below adds the stronger
no-repeated-oriented-segment invariant when an algorithm needs it.
-/
abbrev Path (g : Graph) := Walk g

namespace Walk

/-- Oriented vertices visited by this walk, including the starting vertex. -/
def orientedVertices {g : Graph} (walk : Walk g) : List OrientedSegment :=
  walk.start :: walk.edges.map Link.target

@[simp] theorem orientedVertices_length {g : Graph} (walk : Walk g) :
    walk.orientedVertices.length = walk.edges.length + 1 := by
  simp [orientedVertices]

/-- Reverse-complement traversal of a walk in a bidirected graph. -/
def reverse {g : Graph} (walk : Walk g) (hSym : g.EdgeSymmetric) : Walk g :=
  { start := walk.finish.reverse
    finish := walk.start.reverse
    edges := IsWalk.reverseEdges walk.edges
    valid := IsWalk.reverse hSym walk.valid }

@[simp] theorem reverse_edges {g : Graph} (walk : Walk g)
    (hSym : g.EdgeSymmetric) :
    (walk.reverse hSym).edges = IsWalk.reverseEdges walk.edges := rfl

end Walk

/-- Metadata for a reference-path step before it is connected to graph support. -/
structure PathStep where
  path : PathId
  index : Nat
  segment : SegmentId
  orientation : Orientation
  deriving Repr, DecidableEq, BEq

/-- A named graph path, matching the shape later GFA path semantics will refine. -/
structure NamedPath (g : Graph) where
  id : PathId
  name : String
  sample : Option SampleId
  haplotype : Option HaplotypeId
  walk : Path g
  circular : Bool

/-- A simple path is a walk with no repeated oriented segment visits. -/
structure SimplePath (g : Graph) where
  walk : Walk g
  uniqueVertices : walk.orientedVertices.Nodup

/-- A nonempty closed walk. -/
structure Cycle (g : Graph) where
  root : OrientedSegment
  edges : List Link
  valid : IsWalk g root root edges
  nonempty : edges ≠ []

abbrev CarrierCycle (g : Graph) := Cycle g

def EdgeOnCycle (g : Graph) (edge : Link) : Prop :=
  ∃ cycle : CarrierCycle g, edge ∈ cycle.edges

def CycleEquivalentEdges (g : Graph) (left right : Link) : Prop :=
  ∀ cycle : CarrierCycle g, left ∈ cycle.edges ↔ right ∈ cycle.edges

namespace Graph

def Reachable (g : Graph) (source target : OrientedSegment) : Prop :=
  ∃ edges, IsWalk g source target edges

theorem reachable_refl {g : Graph} {v : OrientedSegment}
    (hv : g.HasOriented v) :
    g.Reachable v v :=
  ⟨[], IsWalk.nil hv⟩

theorem reachable_trans {g : Graph} {a b c : OrientedSegment}
    (hab : g.Reachable a b)
    (hbc : g.Reachable b c) :
    g.Reachable a c := by
  rcases hab with ⟨left, hleft⟩
  rcases hbc with ⟨right, hright⟩
  exact ⟨left ++ right, IsWalk.append hleft hright⟩

end Graph

/-- A semantic connected-component view over oriented vertices. -/
structure ComponentView (g : Graph) where
  id : ComponentId
  vertices : List OrientedSegment
  valid_vertices : ∀ v, v ∈ vertices → g.HasOriented v
  connected : ∀ a, a ∈ vertices → ∀ b, b ∈ vertices → g.Reachable a b

/-- A supported edge cut; separation predicates are supplied by algorithm modules. -/
structure EdgeCut (g : Graph) where
  edges : List Link
  supported : ∀ e, e ∈ edges → g.HasLink e

namespace EdgeCut

def AllBlack {g : Graph} (cut : EdgeCut g) : Prop :=
  ∀ e, e ∈ cut.edges → e.color = EdgeColor.black

def AllReal {g : Graph} (cut : EdgeCut g) : Prop :=
  ∀ e, e ∈ cut.edges → e.provenance = EdgeProvenance.real

end EdgeCut

/--
Minimal hook for algorithmic traversals such as spanning-tree or DFS
intermediate states.  It records only graph support; algorithm modules can add
ordering, dominance, and tree-specific invariants on top.
-/
structure TraversalFrame (g : Graph) where
  roots : List OrientedSegment
  discovered : List OrientedSegment
  treeLinks : List Link
  roots_valid : ∀ v, v ∈ roots → g.HasOriented v
  discovered_valid : ∀ v, v ∈ discovered → g.HasOriented v
  treeLinks_supported : ∀ e, e ∈ treeLinks → g.HasLink e

namespace TraversalFrame

def Spans {g : Graph} (frame : TraversalFrame g) : Prop :=
  ∀ v, g.HasOriented v → v ∈ frame.discovered

def TreeLinksAcyclic {g : Graph} (frame : TraversalFrame g) : Prop :=
  ∀ cycle : Cycle g, ∃ e, e ∈ cycle.edges ∧ e ∉ frame.treeLinks

end TraversalFrame

namespace TODO

/--
Obligation for later algorithm modules: prove that a supported traversal frame
is a true spanning forest for the graph region under consideration.
-/
def TraversalFrameSpanningForest {g : Graph} (frame : TraversalFrame g) : Prop :=
  frame.Spans ∧ frame.TreeLinksAcyclic

end TODO

end Core
end PovuLean
