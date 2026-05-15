/-!
Foundational graph objects for the povu Lean model.

The core graph is intentionally small: it records segments and oriented links,
then exposes predicates for the invariants that downstream GFA, algorithm, and
VCF modules need to preserve.
-/

namespace PovuLean
namespace Core

abbrev SegmentId := Nat
abbrev NodeId := SegmentId
abbrev EdgeId := Nat
abbrev PathId := Nat
abbrev SampleId := Nat
abbrev HaplotypeId := Nat
abbrev ComponentId := Nat

/-- A side of a segment in the biedged carrier view. -/
inductive Side where
  | left
  | right
  deriving Repr, DecidableEq, BEq

namespace Side

def flip : Side → Side
  | left => right
  | right => left

@[simp] theorem flip_flip (side : Side) : flip (flip side) = side := by
  cases side <;> rfl

end Side

/-- A vertex in the biedged carrier graph: one side of one segment. -/
structure SideVertex where
  segment : SegmentId
  side : Side
  deriving Repr, DecidableEq, BEq

/-- GFA-style segment orientation. -/
inductive Orientation where
  | forward
  | reverse
  deriving Repr, DecidableEq, BEq

namespace Orientation

/-- Reverse-complement orientation. -/
def flip : Orientation → Orientation
  | forward => reverse
  | reverse => forward

@[simp] theorem flip_forward : flip forward = reverse := rfl
@[simp] theorem flip_reverse : flip reverse = forward := rfl

@[simp] theorem flip_flip (o : Orientation) : flip (flip o) = o := by
  cases o <;> rfl

end Orientation

/-- A graph segment/node with the sequence string carried for later semantics. -/
structure Segment where
  id : SegmentId
  sequence : String
  deriving Repr, DecidableEq, BEq

/-- A segment as traversed in one orientation. -/
structure OrientedSegment where
  id : SegmentId
  orientation : Orientation
  deriving Repr, DecidableEq, BEq

namespace OrientedSegment

/-- The same segment traversed in the opposite orientation. -/
def reverse (v : OrientedSegment) : OrientedSegment :=
  { id := v.id, orientation := v.orientation.flip }

@[simp] theorem reverse_id (v : OrientedSegment) : v.reverse.id = v.id := rfl

@[simp] theorem reverse_orientation (v : OrientedSegment) :
    v.reverse.orientation = v.orientation.flip := rfl

@[simp] theorem reverse_reverse (v : OrientedSegment) : v.reverse.reverse = v := by
  cases v with
  | mk id orientation =>
      cases orientation <;> rfl

/-- Side where this oriented traversal enters the segment. -/
def startSide (v : OrientedSegment) : SideVertex :=
  match v.orientation with
  | Orientation.forward => { segment := v.id, side := Side.left }
  | Orientation.reverse => { segment := v.id, side := Side.right }

/-- Side where this oriented traversal leaves the segment. -/
def endSide (v : OrientedSegment) : SideVertex :=
  match v.orientation with
  | Orientation.forward => { segment := v.id, side := Side.right }
  | Orientation.reverse => { segment := v.id, side := Side.left }

@[simp] theorem reverse_startSide (v : OrientedSegment) :
    v.reverse.startSide = v.endSide := by
  cases v with
  | mk id orientation =>
      cases orientation <;> rfl

@[simp] theorem reverse_endSide (v : OrientedSegment) :
    v.reverse.endSide = v.startSide := by
  cases v with
  | mk id orientation =>
      cases orientation <;> rfl

end OrientedSegment

/-- Edge color in the biedged view: segment/sequence edges or GFA-link edges. -/
inductive EdgeColor where
  | black
  | grey
  deriving Repr, DecidableEq, BEq

/-- Whether an edge is part of the real input graph or an algorithmic artifact. -/
inductive EdgeProvenance where
  | real
  | dummy
  | virtual
  deriving Repr, DecidableEq, BEq

/--
An oriented edge/link.  A bidirected pangenome graph stores enough links that
`e.reverse` is also traversable whenever `e` is traversable.  `id` is stable for
the directed record; `reverseId` names the record expected for reverse
traversal.
-/
structure Link where
  id : EdgeId
  reverseId : EdgeId
  source : OrientedSegment
  target : OrientedSegment
  color : EdgeColor
  provenance : EdgeProvenance
  deriving Repr, DecidableEq, BEq

/-- Algorithm-facing synonym for a GFA-style oriented link. -/
abbrev Edge := Link

namespace Link

def between (id reverseId : EdgeId)
    (source target : OrientedSegment)
    (color : EdgeColor := EdgeColor.grey)
    (provenance : EdgeProvenance := EdgeProvenance.real) : Link :=
  { id, reverseId, source, target, color, provenance }

/-- Reverse-complement traversal of an oriented link. -/
def reverse (e : Link) : Link :=
  { id := e.reverseId
    reverseId := e.id
    source := e.target.reverse
    target := e.source.reverse
    color := e.color
    provenance := e.provenance }

@[simp] theorem reverse_id (e : Link) : e.reverse.id = e.reverseId := rfl
@[simp] theorem reverse_reverseId (e : Link) : e.reverse.reverseId = e.id := rfl

@[simp] theorem reverse_source (e : Link) : e.reverse.source = e.target.reverse := rfl
@[simp] theorem reverse_target (e : Link) : e.reverse.target = e.source.reverse := rfl
@[simp] theorem reverse_color (e : Link) : e.reverse.color = e.color := rfl
@[simp] theorem reverse_provenance (e : Link) : e.reverse.provenance = e.provenance := rfl

@[simp] theorem reverse_reverse (e : Link) : e.reverse.reverse = e := by
  cases e with
  | mk id reverseId source target color provenance =>
      simp [reverse]

end Link

/-- Directed graph data: segments plus oriented links between segment sides. -/
structure Graph where
  segments : List Segment
  links : List Link
  deriving Repr, DecidableEq, BEq

namespace Graph

def hasSegment (g : Graph) (id : SegmentId) : Bool :=
  g.segments.any (fun segment => segment.id == id)

def HasSegment (g : Graph) (id : SegmentId) : Prop :=
  g.hasSegment id = true

instance decidableHasSegment (g : Graph) (id : SegmentId) :
    Decidable (g.HasSegment id) :=
  inferInstanceAs (Decidable (g.hasSegment id = true))

def hasOriented (g : Graph) (v : OrientedSegment) : Bool :=
  g.hasSegment v.id

def HasOriented (g : Graph) (v : OrientedSegment) : Prop :=
  g.hasOriented v = true

instance decidableHasOriented (g : Graph) (v : OrientedSegment) :
    Decidable (g.HasOriented v) :=
  inferInstanceAs (Decidable (g.hasOriented v = true))

def hasLink (g : Graph) (e : Link) : Bool :=
  g.links.contains e

def HasLink (g : Graph) (e : Link) : Prop :=
  g.hasLink e = true

instance decidableHasLink (g : Graph) (e : Link) :
    Decidable (g.HasLink e) :=
  inferInstanceAs (Decidable (g.hasLink e = true))

def edgeCount (g : Graph) : Nat :=
  g.links.length

def blackEdgeCount (g : Graph) : Nat :=
  (g.links.filter (fun e => e.color == EdgeColor.black)).length

def realEdgeCount (g : Graph) : Nat :=
  (g.links.filter (fun e => e.provenance == EdgeProvenance.real)).length

/-- Segment identifiers are unique. -/
def UniqueSegments (g : Graph) : Prop :=
  (g.segments.map Segment.id).Nodup

instance decidableUniqueSegments (g : Graph) :
    Decidable g.UniqueSegments :=
  inferInstanceAs (Decidable ((g.segments.map Segment.id).Nodup))

/-- Directed edge records have stable, unique ids. -/
def UniqueEdgeIds (g : Graph) : Prop :=
  (g.links.map Link.id).Nodup

instance decidableUniqueEdgeIds (g : Graph) :
    Decidable g.UniqueEdgeIds :=
  inferInstanceAs (Decidable ((g.links.map Link.id).Nodup))

/-- Every listed link points at segments that exist in the graph. -/
def EndpointsValid (g : Graph) : Prop :=
  g.links.all (fun e => g.hasOriented e.source && g.hasOriented e.target) = true

instance decidableEndpointsValid (g : Graph) :
    Decidable g.EndpointsValid :=
  inferInstanceAs
    (Decidable
      (g.links.all (fun e => g.hasOriented e.source && g.hasOriented e.target) = true))

/-- Executable check that every listed link has its reverse-complement link. -/
def BidirectedChecked (g : Graph) : Prop :=
  g.links.all (fun e => g.hasLink e.reverse) = true

instance decidableBidirectedChecked (g : Graph) :
    Decidable g.BidirectedChecked :=
  inferInstanceAs
    (Decidable (g.links.all (fun e => g.hasLink e.reverse) = true))

/-- Logical form used by path reversal lemmas. -/
def EdgeSymmetric (g : Graph) : Prop :=
  ∀ e, g.HasLink e → g.HasLink e.reverse

/-- Directed graph well-formedness independent of bidirected closure. -/
def WellFormed (g : Graph) : Prop :=
  g.UniqueSegments ∧ g.UniqueEdgeIds ∧ g.EndpointsValid

instance decidableWellFormed (g : Graph) :
    Decidable g.WellFormed :=
  inferInstanceAs (Decidable (g.UniqueSegments ∧ g.UniqueEdgeIds ∧ g.EndpointsValid))

/-- Bidirected pangenome graph well-formedness. -/
def BidirectedWellFormed (g : Graph) : Prop :=
  g.WellFormed ∧ g.EdgeSymmetric

@[simp] theorem hasOriented_reverse {g : Graph} {v : OrientedSegment} :
    g.HasOriented v.reverse ↔ g.HasOriented v := by
  cases v with
  | mk id orientation =>
      cases orientation <;> rfl

end Graph

end Core
end PovuLean
