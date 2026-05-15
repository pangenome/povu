import PovuLean.VCF.Correctness
import PovuLean.Algorithms.FlubbleTree.Examples
import PovuLean.Algorithms.Hairpin.Examples

/-!
Executable examples for the semantic VCF layer.

These examples cover SNP-like substitution, insertion, deletion, nested flubble
records, hairpin inversion (`SUBR`), no-variant output, and deterministic
ordering.  They exercise the Lean reference emitter over already verified
variant sources; end-to-end Rust/Lean fixture conformance is intentionally left
to the downstream harness task.
-/

namespace PovuLean
namespace VCF
namespace Examples

def sampleHG1Ref : GenotypeColumn :=
  { sample := "HG1"
    alleles := [GenotypeAllele.ref]
  }

def sampleHG2Alt : GenotypeColumn :=
  { sample := "HG2"
    alleles := [GenotypeAllele.alt 1]
  }

def sampleMissing : GenotypeColumn :=
  { sample := "HG3"
    alleles := [GenotypeAllele.missing]
  }

def simpleFlubbleNode : Algorithms.FlubbleTree.Node :=
  { boundary := Algorithms.Flubble.Examples.simpleBoundary
    parent? := none }

def nestedOuterNode : Algorithms.FlubbleTree.Node :=
  { boundary := Algorithms.Flubble.Examples.nestedOuterBoundary
    parent? := none }

def nestedInnerNode : Algorithms.FlubbleTree.Node :=
  { boundary := Algorithms.Flubble.Examples.nestedInnerBoundary
    parent? := some Algorithms.Flubble.Examples.nestedOuterBoundary }

def snpAlt : AlternateAllele :=
  { construction := AlleleConstruction.substitution "C" "G"
    traversal := ">1"
    count := 1 }

def snpCall : VariantCall :=
  { source := VariantSource.flubble simpleFlubbleNode
    chrom := "chr1"
    contigOrder := 0
    pos := 150
    id := ">0>3"
    ref := "C"
    refTraversal := ">0"
    alternates := [snpAlt]
    variantType := VariantType.sub
    tangled := false
    enclosingSite? := some "flubble:0-2"
    level? := some 0
    referenceAlleleCount := 1
    genotypes := [sampleHG1Ref, sampleHG2Alt, sampleMissing] }

example : snpCall.WellFormed := by
  native_decide

example : (recordOfCall snpCall).WellFormed := by
  exact recordOfCall_wellFormed (by native_decide)

example : (recordOfCall snpCall).alternateAlleles = ["G"] := rfl

example : (recordOfCall snpCall).info.ac = [1] := rfl

example : (recordOfCall snpCall).info.an = 2 := rfl

def insertionAlt : AlternateAllele :=
  { construction := AlleleConstruction.insertion "A" "CCC"
    traversal := ">10>11"
    count := 2 }

def insertionCall : VariantCall :=
  { source := VariantSource.flubble simpleFlubbleNode
    chrom := "chr1"
    contigOrder := 0
    pos := 100
    id := ">10>11"
    ref := "A"
    refTraversal := ">10"
    alternates := [insertionAlt]
    variantType := VariantType.ins
    tangled := false
    enclosingSite? := some "flubble:ins"
    level? := some 0
    referenceAlleleCount := 1
    genotypes := [sampleHG1Ref, sampleHG2Alt] }

example : insertionCall.alternateAlleles = ["ACCC"] := rfl

example : insertionCall.WellFormed := by
  native_decide

def deletionAlt : AlternateAllele :=
  { construction := AlleleConstruction.deletion "A" "CG"
    traversal := ">20"
    count := 1 }

def deletionCall : VariantCall :=
  { source := VariantSource.flubble simpleFlubbleNode
    chrom := "chr1"
    contigOrder := 0
    pos := 200
    id := ">20>23"
    ref := "ACG"
    refTraversal := ">20>21>22"
    alternates := [deletionAlt]
    variantType := VariantType.del
    tangled := false
    enclosingSite? := some "flubble:del"
    level? := some 0
    referenceAlleleCount := 1
    genotypes := [sampleHG1Ref, sampleHG2Alt] }

example : deletionCall.alternateAlleles = ["A"] := rfl

example : deletionCall.WellFormed := by
  native_decide

def nestedOuterAlt : AlternateAllele :=
  { construction := AlleleConstruction.substitution "CGT" "TAA"
    traversal := ">30>31>32"
    count := 1 }

def nestedInnerAlt : AlternateAllele :=
  { construction := AlleleConstruction.insertion "T" "G"
    traversal := ">31>32"
    count := 1 }

def nestedOuterCall : VariantCall :=
  { source := VariantSource.flubble nestedOuterNode
    chrom := "chr2"
    contigOrder := 1
    pos := 200
    id := ">30>35"
    ref := "CGT"
    refTraversal := ">30>31>32"
    alternates := [nestedOuterAlt]
    variantType := VariantType.sub
    tangled := true
    enclosingSite? := some "outer"
    level? := some 0
    referenceAlleleCount := 1
    genotypes := [sampleHG1Ref, sampleHG2Alt] }

def nestedInnerCall : VariantCall :=
  { source := VariantSource.flubble nestedInnerNode
    chrom := "chr2"
    contigOrder := 1
    pos := 210
    id := ">31>33"
    ref := "T"
    refTraversal := ">31"
    alternates := [nestedInnerAlt]
    variantType := VariantType.ins
    tangled := false
    enclosingSite? := some "inner"
    level? := some 1
    referenceAlleleCount := 1
    genotypes := [sampleHG1Ref, sampleHG2Alt] }

def nestedCalls : List VariantCall :=
  [nestedOuterCall, nestedInnerCall]

example : nestedOuterCall.WellFormed := by
  native_decide

example : nestedInnerCall.WellFormed := by
  native_decide

example : OrderedBy VariantCall.orderKey nestedCalls := by
  native_decide

example :
    emitRecords nestedCalls =
      [recordOfCall nestedOuterCall, recordOfCall nestedInnerCall] := rfl

example : OrderedBy Record.orderKey (emitRecords nestedCalls) := by
  exact emitRecords_ordered (by native_decide)

def hairpinAlt : AlternateAllele :=
  { construction := AlleleConstruction.reverseSubstitution "ATGC" "GCAT"
    traversal := "<5<4<3<2"
    count := 1 }

def hairpinCall : VariantCall :=
  { source := VariantSource.hairpin Algorithms.Hairpin.Examples.positiveBoundary
    chrom := "chr3"
    contigOrder := 2
    pos := 300
    id := ">0>1"
    ref := "ATGC"
    refTraversal := ">2>3>4>5"
    alternates := [hairpinAlt]
    variantType := VariantType.subr
    tangled := false
    enclosingSite? := none
    level? := none
    referenceAlleleCount := 1
    genotypes := [sampleHG1Ref, sampleHG2Alt] }

example : hairpinCall.WellFormed := by
  native_decide

example : (recordOfCall hairpinCall).info.variantType = VariantType.subr := rfl

example : (recordOfCall hairpinCall).info.enclosingSite? = none := rfl

def noVariantCalls : List VariantCall := []

example : emitRecords noVariantCalls = [] := rfl

def deterministicCalls : List VariantCall :=
  [insertionCall, snpCall, deletionCall, nestedOuterCall, hairpinCall]

example : OrderedBy VariantCall.orderKey deterministicCalls := by
  native_decide

example : OrderedBy Record.orderKey (emitRecords deterministicCalls) := by
  exact emitRecords_ordered (by native_decide)

end Examples
end VCF
end PovuLean
