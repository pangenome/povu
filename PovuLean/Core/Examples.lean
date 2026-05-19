import PovuLean.Core.Walk

/-!
Small checked graph shapes for downstream modules to import while developing:
simple path, cycle, branch/rejoin bubble, and inversion-like traversal.
-/

namespace PovuLean
namespace Core
namespace Examples

def seg (id : SegmentId) (sequence : String) : Segment :=
  { id, sequence }

def fwd (id : SegmentId) : OrientedSegment :=
  { id, orientation := Orientation.forward }

def rev (id : SegmentId) : OrientedSegment :=
  { id, orientation := Orientation.reverse }

def edge (id reverseId : EdgeId) (source target : OrientedSegment) : Link :=
  Link.between id reverseId source target

def simple01 : Link := edge 0 2 (fwd 0) (fwd 1)
def simple12 : Link := edge 1 3 (fwd 1) (fwd 2)
def simple10rc : Link := edge 2 0 (rev 1) (rev 0)
def simple21rc : Link := edge 3 1 (rev 2) (rev 1)

def simplePathGraph : Graph :=
  { segments := [seg 0 "A", seg 1 "C", seg 2 "G"]
    links := [simple01, simple12, simple10rc, simple21rc] }

example : simplePathGraph.WellFormed := by
  native_decide

example : simplePathGraph.BidirectedChecked := by
  native_decide

def simplePathWalk : Walk simplePathGraph :=
  { start := fwd 0
    finish := fwd 2
    edges := [simple01, simple12]
    valid :=
      IsWalk.cons simple01 rfl rfl (by native_decide) (by native_decide)
        (IsWalk.cons simple12 rfl rfl (by native_decide) (by native_decide)
          (IsWalk.nil (by native_decide))) }

def simpleNamedPath : NamedPath simplePathGraph :=
  { id := 0
    name := "sample0"
    sample := none
    haplotype := none
    walk := simplePathWalk
    circular := false }

example : simplePathWalk.orientedVertices = [fwd 0, fwd 1, fwd 2] := rfl
example : simplePathWalk.edges.length = 2 := rfl

#eval simplePathWalk.edges.length

def cycle01 : Link := edge 0 3 (fwd 0) (fwd 1)
def cycle12 : Link := edge 1 4 (fwd 1) (fwd 2)
def cycle20 : Link := edge 2 5 (fwd 2) (fwd 0)
def cycle10rc : Link := edge 3 0 (rev 1) (rev 0)
def cycle21rc : Link := edge 4 1 (rev 2) (rev 1)
def cycle02rc : Link := edge 5 2 (rev 0) (rev 2)

def cycleGraph : Graph :=
  { segments := [seg 0 "A", seg 1 "C", seg 2 "G"]
    links := [cycle01, cycle12, cycle20, cycle10rc, cycle21rc, cycle02rc] }

example : cycleGraph.WellFormed := by
  native_decide

example : cycleGraph.BidirectedChecked := by
  native_decide

def threeCycle : Cycle cycleGraph :=
  { root := fwd 0
    edges := [cycle01, cycle12, cycle20]
    valid :=
      IsWalk.cons cycle01 rfl rfl (by native_decide) (by native_decide)
        (IsWalk.cons cycle12 rfl rfl (by native_decide) (by native_decide)
          (IsWalk.cons cycle20 rfl rfl (by native_decide) (by native_decide)
            (IsWalk.nil (by native_decide))))
    nonempty := by native_decide }

example : threeCycle.edges.length = 3 := rfl

#eval threeCycle.edges.length

def bubble01 : Link := edge 0 4 (fwd 0) (fwd 1)
def bubble13 : Link := edge 1 5 (fwd 1) (fwd 3)
def bubble02 : Link := edge 2 6 (fwd 0) (fwd 2)
def bubble23 : Link := edge 3 7 (fwd 2) (fwd 3)
def bubble10rc : Link := edge 4 0 (rev 1) (rev 0)
def bubble31rc : Link := edge 5 1 (rev 3) (rev 1)
def bubble20rc : Link := edge 6 2 (rev 2) (rev 0)
def bubble32rc : Link := edge 7 3 (rev 3) (rev 2)

def bubbleGraph : Graph :=
  { segments := [seg 0 "S", seg 1 "A", seg 2 "T", seg 3 "J"]
    links :=
      [ bubble01
      , bubble13
      , bubble02
      , bubble23
      , bubble10rc
      , bubble31rc
      , bubble20rc
      , bubble32rc
      ] }

example : bubbleGraph.WellFormed := by
  native_decide

example : bubbleGraph.BidirectedChecked := by
  native_decide

def bubbleLeftWalk : Walk bubbleGraph :=
  { start := fwd 0
    finish := fwd 3
    edges := [bubble01, bubble13]
    valid :=
      IsWalk.cons bubble01 rfl rfl (by native_decide) (by native_decide)
        (IsWalk.cons bubble13 rfl rfl (by native_decide) (by native_decide)
          (IsWalk.nil (by native_decide))) }

def bubbleRightWalk : Walk bubbleGraph :=
  { start := fwd 0
    finish := fwd 3
    edges := [bubble02, bubble23]
    valid :=
      IsWalk.cons bubble02 rfl rfl (by native_decide) (by native_decide)
        (IsWalk.cons bubble23 rfl rfl (by native_decide) (by native_decide)
          (IsWalk.nil (by native_decide))) }

example : bubbleLeftWalk.start = bubbleRightWalk.start := rfl
example : bubbleLeftWalk.finish = bubbleRightWalk.finish := rfl
example : bubbleLeftWalk.edges ≠ bubbleRightWalk.edges := by native_decide

#eval bubbleLeftWalk.edges.length
#eval bubbleRightWalk.edges.length

def inversion0r1 : Link := edge 0 2 (fwd 0) (rev 1)
def inversionr12 : Link := edge 1 3 (rev 1) (fwd 2)
def inversion1r0rc : Link := edge 2 0 (fwd 1) (rev 0)
def inversionr21rc : Link := edge 3 1 (rev 2) (fwd 1)

def inversionGraph : Graph :=
  { segments := [seg 0 "L", seg 1 "INV", seg 2 "R"]
    links := [inversion0r1, inversionr12, inversion1r0rc, inversionr21rc] }

example : inversionGraph.WellFormed := by
  native_decide

example : inversionGraph.BidirectedChecked := by
  native_decide

def inversionLikeWalk : Walk inversionGraph :=
  { start := fwd 0
    finish := fwd 2
    edges := [inversion0r1, inversionr12]
    valid :=
      IsWalk.cons inversion0r1 rfl rfl (by native_decide) (by native_decide)
        (IsWalk.cons inversionr12 rfl rfl (by native_decide) (by native_decide)
          (IsWalk.nil (by native_decide))) }

example : inversionLikeWalk.orientedVertices = [fwd 0, rev 1, fwd 2] := rfl

#eval inversionLikeWalk.orientedVertices

end Examples
end Core
end PovuLean
