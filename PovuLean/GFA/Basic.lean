import PovuLean.Core.Walk

/-!
Semantic GFA boundary accepted by the Lean povu model.

This module intentionally does not parse bytes.  A trusted parser/normalizer is
expected to resolve textual GFA names to stable natural ids and then supply the
records below.  The verified part begins at those normalized records and checks
that they construct the core graph model with the invariants required by later
proofs.

Accepted subset and normalization choices:

* `H` lines are accepted as metadata and ignored by this semantic graph model.
* `S` lines become `SegmentRecord`s.  Segment names and normalized ids must be
  unique, names must be nonempty, and sequences must be concrete strings rather
  than the unknown-sequence marker `*`.
* `L` lines become `LinkRecord`s.  Only exact zero-overlap links (`0M`) are in
  scope.  Each semantic link is normalized into two real grey core links: the
  listed traversal and its reverse-complement traversal.
* `P` and `W` lines are normalized to `PathRecord`s.  A path is accepted only
  when its steps are nonempty, refer to declared segments, and a corresponding
  `Core.Walk` exists in the constructed graph.
* Optional tags are outside this model.  They may be retained by a parser, but
  they do not affect graph/path construction here.
* Containment (`C`), jump (`J`), and all other record kinds are rejected at this
  semantic boundary.  Low-level syntax errors are also outside this model and
  should be rejected before a `Document` is produced.
-/

namespace PovuLean
namespace GFA

open Core

abbrev SegmentName := String
abbrev PathName := String
abbrev Cigar := String

/-- GFA record classes relevant to the semantic boundary. -/
inductive LineKind where
  | header
  | segment
  | link
  | path
  | walk
  | containment
  | jump
  | other
  deriving Repr, DecidableEq, BEq

namespace LineKind

/-- Record kinds accepted before normalization. -/
def Supported : LineKind → Prop
  | header => True
  | segment => True
  | link => True
  | path => True
  | walk => True
  | containment => False
  | jump => False
  | other => False

instance decidableSupported (kind : LineKind) : Decidable kind.Supported := by
  cases kind
  · exact isTrue trivial
  · exact isTrue trivial
  · exact isTrue trivial
  · exact isTrue trivial
  · exact isTrue trivial
  · exact isFalse (by intro h; exact h)
  · exact isFalse (by intro h; exact h)
  · exact isFalse (by intro h; exact h)

end LineKind

/-- Local lawful-BEq instances needed to connect executable graph checks to
logical membership.  They are kept in the GFA layer instead of refactoring the
core graph files. -/
instance lawfulBEqOrientation : LawfulBEq Orientation where
  rfl := by
    intro orientation
    cases orientation <;> rfl
  eq_of_beq := by
    intro a b h
    cases a <;> cases b
    · rfl
    · cases h
    · cases h
    · rfl

instance lawfulBEqOrientedSegment : LawfulBEq OrientedSegment where
  rfl := by
    intro v
    cases v with
    | mk id orientation =>
        change (id == id && orientation == orientation) = true
        simp
  eq_of_beq := by
    intro a b h
    cases a with
    | mk aId aOrientation =>
    cases b with
    | mk bId bOrientation =>
        change (aId == bId && aOrientation == bOrientation) = true at h
        simp at h
        rcases h with ⟨hId, hOrientation⟩
        cases hId
        cases hOrientation
        rfl

instance lawfulBEqEdgeColor : LawfulBEq EdgeColor where
  rfl := by
    intro color
    cases color <;> rfl
  eq_of_beq := by
    intro a b h
    cases a <;> cases b
    · rfl
    · cases h
    · cases h
    · rfl

instance lawfulBEqEdgeProvenance : LawfulBEq EdgeProvenance where
  rfl := by
    intro provenance
    cases provenance <;> rfl
  eq_of_beq := by
    intro a b h
    cases a <;> cases b
    · rfl
    · cases h
    · cases h
    · cases h
    · rfl
    · cases h
    · cases h
    · cases h
    · rfl

instance lawfulBEqLink : LawfulBEq Link where
  rfl := by
    intro e
    cases e with
    | mk id reverseId source target color provenance =>
        change
          ( id == id
            && ( reverseId == reverseId
              && ( source == source
                && ( target == target
                  && (color == color && provenance == provenance))))) = true
        simp
  eq_of_beq := by
    intro a b h
    cases a with
    | mk aId aReverseId aSource aTarget aColor aProvenance =>
    cases b with
    | mk bId bReverseId bSource bTarget bColor bProvenance =>
        change
          ( aId == bId
            && ( aReverseId == bReverseId
              && ( aSource == bSource
                && ( aTarget == bTarget
                  && (aColor == bColor && aProvenance == bProvenance))))) = true at h
        simp at h
        rcases h with ⟨hId, hReverseId, hSource, hTarget, hColor, hProvenance⟩
        cases hId
        cases hReverseId
        cases hSource
        cases hTarget
        cases hColor
        cases hProvenance
        rfl

/-- Executable bidirected closure implies the logical symmetry predicate used
by core walk-reversal lemmas. -/
theorem edgeSymmetric_of_bidirectedChecked {g : Graph}
    (h : g.BidirectedChecked) :
    g.EdgeSymmetric := by
  intro e he
  have hmem : e ∈ g.links := by
    exact List.contains_iff_mem.mp he
  exact (List.all_eq_true.mp h) e hmem

/-- A normalized GFA segment. -/
structure SegmentRecord where
  id : SegmentId
  name : SegmentName
  sequence : String
  deriving Repr, DecidableEq, BEq

namespace SegmentRecord

def toCore (record : SegmentRecord) : Segment :=
  { id := record.id, sequence := record.sequence }

/-- Segment-level subset constraints independent of other records. -/
def Supported (record : SegmentRecord) : Prop :=
  record.name ≠ "" ∧ record.sequence ≠ "*"

instance decidableSupported (record : SegmentRecord) :
    Decidable record.Supported :=
  inferInstanceAs (Decidable (record.name ≠ "" ∧ record.sequence ≠ "*"))

end SegmentRecord

/-- A normalized GFA link with parser-assigned directed edge ids. -/
structure LinkRecord where
  id : EdgeId
  reverseId : EdgeId
  source : OrientedSegment
  target : OrientedSegment
  overlap : Cigar
  deriving Repr, DecidableEq, BEq

namespace LinkRecord

/-- Link-level subset constraints independent of the segment table. -/
def Supported (record : LinkRecord) : Prop :=
  record.overlap = "0M"

instance decidableSupported (record : LinkRecord) :
    Decidable record.Supported :=
  inferInstanceAs (Decidable (record.overlap = "0M"))

/-- The listed oriented traversal as a real grey core link. -/
def toCoreLink (record : LinkRecord) : Link :=
  Link.between record.id record.reverseId record.source record.target
    EdgeColor.grey EdgeProvenance.real

/-- The bidirected pair contributed by one normalized GFA `L` line. -/
def toCoreLinks (record : LinkRecord) : List Link :=
  [record.toCoreLink, record.toCoreLink.reverse]

@[simp] theorem toCoreLinks_length (record : LinkRecord) :
    record.toCoreLinks.length = 2 := rfl

@[simp] theorem reverse_mem_toCoreLinks (record : LinkRecord) :
    record.toCoreLink.reverse ∈ record.toCoreLinks := by
  simp [toCoreLinks]

end LinkRecord

/-- A normalized GFA `P`/`W` path. -/
structure PathRecord where
  id : PathId
  name : PathName
  sample : Option SampleId
  haplotype : Option HaplotypeId
  steps : List OrientedSegment
  circular : Bool
  deriving Repr, DecidableEq, BEq

namespace PathRecord

/-- Path-level subset constraints independent of graph support. -/
def WellNamed (record : PathRecord) : Prop :=
  record.name ≠ "" ∧ record.steps ≠ []

instance decidableWellNamed (record : PathRecord) :
    Decidable record.WellNamed :=
  inferInstanceAs (Decidable (record.name ≠ "" ∧ record.steps ≠ []))

end PathRecord

/-- Normalized semantic GFA document consumed by the Lean model. -/
structure Document where
  segments : List SegmentRecord
  links : List LinkRecord
  paths : List PathRecord
  deriving Repr, DecidableEq, BEq

namespace Document

def segmentIds (doc : Document) : List SegmentId :=
  doc.segments.map SegmentRecord.id

def segmentNames (doc : Document) : List SegmentName :=
  doc.segments.map SegmentRecord.name

def pathIds (doc : Document) : List PathId :=
  doc.paths.map PathRecord.id

def pathNames (doc : Document) : List PathName :=
  doc.paths.map PathRecord.name

def hasSegmentId (doc : Document) (id : SegmentId) : Bool :=
  doc.segmentIds.contains id

def HasSegmentId (doc : Document) (id : SegmentId) : Prop :=
  doc.hasSegmentId id = true

instance decidableHasSegmentId (doc : Document) (id : SegmentId) :
    Decidable (doc.HasSegmentId id) :=
  inferInstanceAs (Decidable (doc.hasSegmentId id = true))

def SegmentsSupported (doc : Document) : Prop :=
  doc.segments.all (fun segment => decide segment.Supported) = true

def LinksSupported (doc : Document) : Prop :=
  doc.links.all (fun link => decide link.Supported) = true

def LinkEndpointsDeclared (doc : Document) : Prop :=
  doc.links.all
    (fun link => doc.hasSegmentId link.source.id && doc.hasSegmentId link.target.id) = true

def PathsWellNamed (doc : Document) : Prop :=
  doc.paths.all (fun path => decide path.WellNamed) = true

def PathStepsDeclared (doc : Document) : Prop :=
  doc.paths.all
    (fun path => path.steps.all (fun step => doc.hasSegmentId step.id)) = true

def SegmentIdsUnique (doc : Document) : Prop :=
  doc.segmentIds.Nodup

def SegmentNamesUnique (doc : Document) : Prop :=
  doc.segmentNames.Nodup

def PathIdsUnique (doc : Document) : Prop :=
  doc.pathIds.Nodup

def PathNamesUnique (doc : Document) : Prop :=
  doc.pathNames.Nodup

instance decidableSegmentsSupported (doc : Document) :
    Decidable doc.SegmentsSupported :=
  inferInstanceAs
    (Decidable (doc.segments.all (fun segment => decide segment.Supported) = true))

instance decidableLinksSupported (doc : Document) :
    Decidable doc.LinksSupported :=
  inferInstanceAs
    (Decidable (doc.links.all (fun link => decide link.Supported) = true))

instance decidableLinkEndpointsDeclared (doc : Document) :
    Decidable doc.LinkEndpointsDeclared :=
  inferInstanceAs
    (Decidable
      (doc.links.all
        (fun link => doc.hasSegmentId link.source.id && doc.hasSegmentId link.target.id) = true))

instance decidablePathsWellNamed (doc : Document) :
    Decidable doc.PathsWellNamed :=
  inferInstanceAs
    (Decidable (doc.paths.all (fun path => decide path.WellNamed) = true))

instance decidablePathStepsDeclared (doc : Document) :
    Decidable doc.PathStepsDeclared :=
  inferInstanceAs
    (Decidable
      (doc.paths.all
        (fun path => path.steps.all (fun step => doc.hasSegmentId step.id)) = true))

instance decidableSegmentIdsUnique (doc : Document) :
    Decidable doc.SegmentIdsUnique :=
  inferInstanceAs (Decidable doc.segmentIds.Nodup)

instance decidableSegmentNamesUnique (doc : Document) :
    Decidable doc.SegmentNamesUnique :=
  inferInstanceAs (Decidable doc.segmentNames.Nodup)

instance decidablePathIdsUnique (doc : Document) :
    Decidable doc.PathIdsUnique :=
  inferInstanceAs (Decidable doc.pathIds.Nodup)

instance decidablePathNamesUnique (doc : Document) :
    Decidable doc.PathNamesUnique :=
  inferInstanceAs (Decidable doc.pathNames.Nodup)

/--
Semantic subset accepted before graph construction.  This captures rejected
semantic cases such as duplicate segment names, unknown sequences, nonzero link
overlap, dangling link endpoints, unnamed/empty paths, and path steps that
reference undeclared segments.
-/
def SubsetAccepted (doc : Document) : Prop :=
  doc.SegmentsSupported
    ∧ doc.SegmentIdsUnique
    ∧ doc.SegmentNamesUnique
    ∧ doc.LinksSupported
    ∧ doc.LinkEndpointsDeclared
    ∧ doc.PathsWellNamed
    ∧ doc.PathStepsDeclared
    ∧ doc.PathIdsUnique
    ∧ doc.PathNamesUnique

instance decidableSubsetAccepted (doc : Document) :
    Decidable doc.SubsetAccepted :=
  inferInstanceAs
    (Decidable
      ( doc.SegmentsSupported
        ∧ doc.SegmentIdsUnique
        ∧ doc.SegmentNamesUnique
        ∧ doc.LinksSupported
        ∧ doc.LinkEndpointsDeclared
        ∧ doc.PathsWellNamed
        ∧ doc.PathStepsDeclared
        ∧ doc.PathIdsUnique
        ∧ doc.PathNamesUnique))

def toGraph (doc : Document) : Graph :=
  { segments := doc.segments.map SegmentRecord.toCore
    links := doc.links.flatMap LinkRecord.toCoreLinks }

/--
Executable construction check for the graph portion of a normalized document.
Path support is proof-producing and appears in `Accepted` below.
-/
def CheckedCore (doc : Document) : Prop :=
  doc.SubsetAccepted ∧ doc.toGraph.WellFormed ∧ doc.toGraph.BidirectedChecked

instance decidableCheckedCore (doc : Document) :
    Decidable doc.CheckedCore :=
  inferInstanceAs
    (Decidable (doc.SubsetAccepted ∧ doc.toGraph.WellFormed ∧ doc.toGraph.BidirectedChecked))

/-- A normalized path is supported when it is exactly the vertex trace of a core
walk over the constructed graph. -/
def PathSupported (doc : Document) (path : PathRecord) : Prop :=
  ∃ walk : Walk doc.toGraph, walk.orientedVertices = path.steps

/-- Full accepted semantic document: executable graph checks plus supported
paths. -/
structure Accepted (doc : Document) : Prop where
  checkedCore : doc.CheckedCore
  pathsSupported : ∀ path, path ∈ doc.paths → doc.PathSupported path

/--
Main construction theorem: a semantic GFA document accepted by this boundary
constructs a bidirected core graph satisfying the graph invariants consumed by
walk/path and algorithm proofs.
-/
theorem checkedCore_toGraph_bidirectedWellFormed {doc : Document}
    (h : doc.CheckedCore) :
    doc.toGraph.BidirectedWellFormed :=
  ⟨h.2.1, edgeSymmetric_of_bidirectedChecked h.2.2⟩

theorem accepted_toGraph_bidirectedWellFormed {doc : Document}
    (h : doc.Accepted) :
    doc.toGraph.BidirectedWellFormed :=
  checkedCore_toGraph_bidirectedWellFormed h.checkedCore

/-- Accepted documents preserve path support as named core paths. -/
noncomputable def toNamedPath {doc : Document} (path : PathRecord)
    (support : doc.PathSupported path) :
    NamedPath doc.toGraph :=
  { id := path.id
    name := path.name
    sample := path.sample
    haplotype := path.haplotype
    walk := support.choose
    circular := path.circular }

@[simp] theorem toNamedPath_walk_vertices {doc : Document} (path : PathRecord)
    (support : doc.PathSupported path) :
    (doc.toNamedPath path support).walk.orientedVertices = path.steps :=
  support.choose_spec

end Document

end GFA
end PovuLean
