import PovuLean.GFA.Basic

/-!
Checked semantic-boundary examples for the GFA model.

The accepted examples cover ordinary forward links, reverse-oriented links, and
a supported path.  The rejected examples cover malformed or unsupported inputs
that should not cross from a byte parser into the verified graph construction
without being rejected first.
-/

namespace PovuLean
namespace GFA
namespace Examples

open Core

def segment (id : SegmentId) (name sequence : String) : SegmentRecord :=
  { id, name, sequence }

def fwd (id : SegmentId) : OrientedSegment :=
  { id, orientation := Orientation.forward }

def rev (id : SegmentId) : OrientedSegment :=
  { id, orientation := Orientation.reverse }

def link (id reverseId : EdgeId)
    (source target : OrientedSegment)
    (overlap : Cigar := "0M") : LinkRecord :=
  { id, reverseId, source, target, overlap }

example : LineKind.segment.Supported := by
  native_decide

example : LineKind.link.Supported := by
  native_decide

example : ¬ LineKind.containment.Supported := by
  native_decide

example : ¬ LineKind.jump.Supported := by
  native_decide

def tinySegment0 : SegmentRecord := segment 0 "s0" "A"
def tinySegment1 : SegmentRecord := segment 1 "s1" "C"
def tinyLink01 : LinkRecord := link 0 1 (fwd 0) (fwd 1)

def tinyPath : PathRecord :=
  { id := 0
    name := "ref"
    sample := none
    haplotype := none
    steps := [fwd 0, fwd 1]
    circular := false }

/-- Representative accepted S/L/P semantic records. -/
def tinyDoc : Document :=
  { segments := [tinySegment0, tinySegment1]
    links := [tinyLink01]
    paths := [tinyPath] }

example : tinyDoc.SubsetAccepted := by
  native_decide

example : tinyDoc.CheckedCore := by
  native_decide

example : tinyDoc.toGraph.WellFormed := by
  native_decide

example : tinyDoc.toGraph.BidirectedChecked := by
  native_decide

example : tinyDoc.toGraph.edgeCount = 2 := rfl

example : tinyDoc.toGraph.realEdgeCount = 2 := rfl

example : tinyDoc.toGraph.BidirectedWellFormed :=
  Document.checkedCore_toGraph_bidirectedWellFormed (by native_decide : tinyDoc.CheckedCore)

def tinyForwardLink : Link := tinyLink01.toCoreLink

def tinyForwardWalk : Walk tinyDoc.toGraph :=
  { start := fwd 0
    finish := fwd 1
    edges := [tinyForwardLink]
    valid :=
      IsWalk.cons tinyForwardLink rfl rfl
        (by native_decide)
        (by native_decide)
        (IsWalk.nil (by native_decide)) }

example : tinyForwardWalk.orientedVertices = [fwd 0, fwd 1] := rfl

def tinyPathSupport : tinyDoc.PathSupported tinyPath :=
  ⟨tinyForwardWalk, rfl⟩

def tinyAccepted : tinyDoc.Accepted :=
  { checkedCore := by native_decide
    pathsSupported := by
      intro path hPath
      simp [tinyDoc] at hPath
      cases hPath
      exact tinyPathSupport }

example : tinyDoc.toGraph.BidirectedWellFormed :=
  Document.accepted_toGraph_bidirectedWellFormed tinyAccepted

example :
    (tinyDoc.toNamedPath tinyPath tinyPathSupport).walk.orientedVertices =
      [fwd 0, fwd 1] := by
  simp [tinyPath]

#eval tinyDoc.toGraph.edgeCount
#eval tinyForwardWalk.orientedVertices

def reverseLink : LinkRecord := link 2 3 (fwd 0) (rev 1)

/-- Accepted reverse-oriented link, used by inversion-like traversals. -/
def reverseOrientedDoc : Document :=
  { segments := [tinySegment0, tinySegment1]
    links := [reverseLink]
    paths := [] }

example : reverseOrientedDoc.CheckedCore := by
  native_decide

example : reverseOrientedDoc.toGraph.BidirectedWellFormed :=
  Document.checkedCore_toGraph_bidirectedWellFormed
    (by native_decide : reverseOrientedDoc.CheckedCore)

def duplicateSegmentDoc : Document :=
  { segments := [tinySegment0, segment 0 "s0-duplicate" "G"]
    links := []
    paths := [] }

example : ¬ duplicateSegmentDoc.SubsetAccepted := by
  native_decide

def unknownSequenceDoc : Document :=
  { segments := [segment 0 "s0" "*"]
    links := []
    paths := [] }

example : ¬ unknownSequenceDoc.SubsetAccepted := by
  native_decide

def nonZeroOverlapDoc : Document :=
  { segments := [tinySegment0, tinySegment1]
    links := [link 0 1 (fwd 0) (fwd 1) "1M"]
    paths := [] }

example : ¬ nonZeroOverlapDoc.SubsetAccepted := by
  native_decide

def danglingLinkDoc : Document :=
  { segments := [tinySegment0]
    links := [tinyLink01]
    paths := [] }

example : ¬ danglingLinkDoc.SubsetAccepted := by
  native_decide

example : ¬ danglingLinkDoc.CheckedCore := by
  native_decide

def emptyPathDoc : Document :=
  { segments := [tinySegment0]
    links := []
    paths :=
      [ { id := 0
          name := "empty"
          sample := none
          haplotype := none
          steps := []
          circular := false } ] }

example : ¬ emptyPathDoc.SubsetAccepted := by
  native_decide

end Examples
end GFA
end PovuLean
