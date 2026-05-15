import PovuLean.Pipeline
import PovuLean.VCF.Emit
import PovuLean.Algorithms.FlubbleTree.Examples
import PovuLean.Algorithms.Hairpin.Examples

/-!
Lean-side semantic VCF fixture emitter for the Rust conformance harness.

The executable keeps fixture expectations inside the Lean semantic boundary:
each row is built as a `VCF.VariantCall` and serialized only after
`VCF.emitRecords` maps it to semantic `VCF.Record`s.  The Rust harness compares
normalized semantics rather than trusting byte-for-byte VCF header formatting.
-/

namespace Lean4Conformance

open PovuLean

def joinSep (sep : String) : List String → String
  | [] => ""
  | first :: rest => rest.foldl (fun acc value => acc ++ sep ++ value) first

def joinComma (values : List String) : String :=
  joinSep "," values

def joinTab (values : List String) : String :=
  joinSep "\t" values

def natList (values : List Nat) : String :=
  joinComma (values.map toString)

def formatAf (freq : VCF.AlleleFrequency) : String :=
  if freq.totalAlleles == 0 then
    "NaN"
  else if freq.alternateCount == 1 && freq.totalAlleles == 2 then
    "0.5"
  else
    toString freq.alternateCount ++ "/" ++ toString freq.totalAlleles

def formatBoolInfo (value : Bool) : String :=
  if value then "T" else "F"

def formatGenotypeAllele : VCF.GenotypeAllele → String
  | .missing => "."
  | .ref => "0"
  | .alt index => toString index

def formatGenotypeColumn (column : VCF.GenotypeColumn) : String :=
  joinSep "|" (column.alleles.map formatGenotypeAllele)

def infoFields (info : VCF.Info) : List String :=
  let base :=
    [ "AC=" ++ natList info.ac
    , "AF=" ++ joinComma (info.af.map formatAf)
    , "AN=" ++ toString info.an
    , "NS=" ++ toString info.ns
    , "AT=" ++ joinComma info.atField
    , "VARTYPE=" ++ VCF.VariantType.vcfLabel info.variantType
    , "TANGLED=" ++ formatBoolInfo info.tangled
    ]
  let withEs :=
    match info.enclosingSite? with
    | none => base
    | some site => base ++ ["ES=" ++ site]
  match info.level? with
  | none => withEs
  | some level => withEs ++ ["LV=" ++ toString level]

def recordLine (record : VCF.Record) : String :=
  joinTab
    ([ record.chrom
     , toString record.pos
     , record.id
     , record.ref
     , joinComma (VCF.Record.alternateAlleles record)
     , record.qual
     , record.filter
     , joinSep ";" (infoFields record.info)
     , record.format
     ] ++ record.genotypes.map formatGenotypeColumn)

def vcfText (samples : List String) (records : List VCF.Record) : String :=
  let header :=
    joinTab
      (["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        ++ samples)
  joinSep "\n" (header :: records.map recordLine) ++ "\n"

def refSample (sample : String) : VCF.GenotypeColumn :=
  { sample := sample
    alleles := [VCF.GenotypeAllele.ref] }

def altSample (sample : String) : VCF.GenotypeColumn :=
  { sample := sample
    alleles := [VCF.GenotypeAllele.alt 1] }

def missingSample (sample : String) : VCF.GenotypeColumn :=
  { sample := sample
    alleles := [VCF.GenotypeAllele.missing] }

def sampleHG1Ref : VCF.GenotypeColumn :=
  refSample "HG1"

def sampleHG2Alt : VCF.GenotypeColumn :=
  altSample "HG2"

def sampleHG3Missing : VCF.GenotypeColumn :=
  missingSample "HG3"

def sampleRefRef : VCF.GenotypeColumn :=
  refSample "ref"

def sampleAltAlt : VCF.GenotypeColumn :=
  altSample "alt"

def minimalSubstitutionNode : Algorithms.FlubbleTree.Node :=
  { boundary := Algorithms.Flubble.Examples.simpleBoundary
    parent? := none }

def nestedInnerNode : Algorithms.FlubbleTree.Node :=
  { boundary := Algorithms.Flubble.Examples.nestedInnerBoundary
    parent? := some Algorithms.Flubble.Examples.nestedOuterBoundary }

def minimalSubstitutionAlt : VCF.AlternateAllele :=
  { construction := VCF.AlleleConstruction.substitution "C" "G"
    traversal := ">2"
    count := 1 }

/--
Semantic expectation for `fixtures/minimal_substitution.gfa`.

The byte-level GFA describes two haploid paths, `HG1#1#chr1` and
`HG2#1#chr1`, that diverge from segment `0` to either `1` or `2` and rejoin at
segment `3`.  Current povu emits the `HG1` path as reference and `HG2` as the
alternate substitution.
-/
def minimalSubstitutionCall : VCF.VariantCall :=
  { source := VCF.VariantSource.flubble minimalSubstitutionNode
    chrom := "HG1#1#chr1"
    contigOrder := 0
    pos := 2
    id := ">0>3"
    ref := "C"
    refTraversal := ">1"
    alternates := [minimalSubstitutionAlt]
    variantType := VCF.VariantType.sub
    tangled := false
    enclosingSite? := some ">0>3"
    level? := some 0
    referenceAlleleCount := 1
    genotypes := [sampleHG1Ref, sampleHG2Alt] }

example : minimalSubstitutionCall.WellFormed := by
  native_decide

/--
Semantic expectation for `fixtures/insertion_flubble.gfa`.

The alternate path traverses segment `2` between the reference anchor segment
`0` and rejoin segment `1`, so current povu emits an anchored `INS`.
-/
def insertionFlubbleAlt : VCF.AlternateAllele :=
  { construction := VCF.AlleleConstruction.insertion "A" "G"
    traversal := ">0>2"
    count := 1 }

def insertionFlubbleCall : VCF.VariantCall :=
  { source := VCF.VariantSource.flubble minimalSubstitutionNode
    chrom := "HG1#1#chr1"
    contigOrder := 0
    pos := 1
    id := ">0>1"
    ref := "A"
    refTraversal := ">0"
    alternates := [insertionFlubbleAlt]
    variantType := VCF.VariantType.ins
    tangled := false
    enclosingSite? := some ">0>1"
    level? := some 0
    referenceAlleleCount := 1
    genotypes := [sampleHG1Ref, sampleHG2Alt] }

example : insertionFlubbleCall.WellFormed := by
  native_decide

/--
Semantic expectation for `fixtures/deletion_flubble.gfa`.

The reference path contains segment `2` after anchor segment `0`; the alternate
path skips it, so current povu emits an anchored `DEL`.
-/
def deletionFlubbleAlt : VCF.AlternateAllele :=
  { construction := VCF.AlleleConstruction.deletion "A" "G"
    traversal := ">0"
    count := 1 }

def deletionFlubbleCall : VCF.VariantCall :=
  { source := VCF.VariantSource.flubble minimalSubstitutionNode
    chrom := "HG1#1#chr1"
    contigOrder := 0
    pos := 1
    id := ">0>1"
    ref := "AG"
    refTraversal := ">0>2"
    alternates := [deletionFlubbleAlt]
    variantType := VCF.VariantType.del
    tangled := false
    enclosingSite? := some ">0>1"
    level? := some 0
    referenceAlleleCount := 1
    genotypes := [sampleHG1Ref, sampleHG2Alt] }

example : deletionFlubbleCall.WellFormed := by
  native_decide

/--
Semantic expectation for `fixtures/nested_deletion.gfa`.

The graph has an outer branch and an inner branch on the reference side.  The
current end-to-end implementation emits the inner deletion at `LV=1`; the outer
sample that takes the sibling branch is represented as missing for this record.
-/
def nestedDeletionAlt : VCF.AlternateAllele :=
  { construction := VCF.AlleleConstruction.deletion "C" "T"
    traversal := ">1"
    count := 1 }

def nestedDeletionCall : VCF.VariantCall :=
  { source := VCF.VariantSource.flubble nestedInnerNode
    chrom := "HG1#1#chr1"
    contigOrder := 0
    pos := 2
    id := ">1>4"
    ref := "CT"
    refTraversal := ">1>3"
    alternates := [nestedDeletionAlt]
    variantType := VCF.VariantType.del
    tangled := false
    enclosingSite? := some ">1>4"
    level? := some 1
    referenceAlleleCount := 1
    genotypes := [sampleHG1Ref, sampleHG2Alt, sampleHG3Missing] }

example : nestedDeletionCall.WellFormed := by
  native_decide

/--
Semantic expectation for `fixtures/hairpin_inversion_subr.gfa`.

The alternate path traverses the same five labelled segments in reverse
orientation through the `1` to `5` region, matching povu's `SUBR` convention.
-/
def hairpinInversionAlt : VCF.AlternateAllele :=
  { construction := VCF.AlleleConstruction.reverseSubstitution "ACGTA" "TACGT"
    traversal := "<5<4<3<2<1"
    count := 1 }

def hairpinInversionCall : VCF.VariantCall :=
  { source := VCF.VariantSource.hairpin Algorithms.Hairpin.Examples.positiveBoundary
    chrom := "ref"
    contigOrder := 0
    pos := 2
    id := ">1>5"
    ref := "ACGTA"
    refTraversal := ">1>2>3>4>5"
    alternates := [hairpinInversionAlt]
    variantType := VCF.VariantType.subr
    tangled := false
    enclosingSite? := none
    level? := none
    referenceAlleleCount := 1
    genotypes := [sampleRefRef, sampleAltAlt] }

example : hairpinInversionCall.WellFormed := by
  native_decide

/--
Semantic expectation for `fixtures/two_ordered_substitutions.gfa`.

The graph contains two independent top-level flubbles on one contig.  The call
list is intentionally ordered by position so the Rust harness can require raw
record-order agreement for this fixture.
-/
def secondOrderedSubstitutionAlt : VCF.AlternateAllele :=
  { construction := VCF.AlleleConstruction.substitution "A" "C"
    traversal := ">5"
    count := 1 }

def secondOrderedSubstitutionCall : VCF.VariantCall :=
  { source := VCF.VariantSource.flubble minimalSubstitutionNode
    chrom := "HG1#1#chr1"
    contigOrder := 0
    pos := 4
    id := ">3>6"
    ref := "A"
    refTraversal := ">4"
    alternates := [secondOrderedSubstitutionAlt]
    variantType := VCF.VariantType.sub
    tangled := false
    enclosingSite? := some ">3>6"
    level? := some 0
    referenceAlleleCount := 1
    genotypes := [sampleHG1Ref, sampleHG2Alt] }

example : secondOrderedSubstitutionCall.WellFormed := by
  native_decide

def twoOrderedSubstitutionCalls : List VCF.VariantCall :=
  [minimalSubstitutionCall, secondOrderedSubstitutionCall]

def fixtureOutput? : String → Option String
  | "minimal-substitution" =>
      some (vcfText ["HG1", "HG2"] (VCF.emitRecords [minimalSubstitutionCall]))
  | "insertion-flubble" =>
      some (vcfText ["HG1", "HG2"] (VCF.emitRecords [insertionFlubbleCall]))
  | "deletion-flubble" =>
      some (vcfText ["HG1", "HG2"] (VCF.emitRecords [deletionFlubbleCall]))
  | "nested-deletion" =>
      some (vcfText ["HG1", "HG2", "HG3"] (VCF.emitRecords [nestedDeletionCall]))
  | "hairpin-inversion-subr" =>
      some (vcfText ["ref", "alt"] (VCF.emitRecords [hairpinInversionCall]))
  | "linear-no-variant" =>
      some (vcfText ["HG1", "HG2"] (VCF.emitRecords []))
  | "two-ordered-substitutions" =>
      some (vcfText ["HG1", "HG2"] (VCF.emitRecords twoOrderedSubstitutionCalls))
  | _ => none

def runMain (args : List String) : IO UInt32 := do
  match args with
  | [fixtureId] =>
      match fixtureOutput? fixtureId with
      | some output =>
          IO.print output
          return 0
      | none =>
          IO.eprintln s!"unknown Lean conformance fixture: {fixtureId}"
          return 2
  | _ =>
      IO.eprintln "usage: lake env lean --run tests/lean4_conformance/lean_reference.lean <fixture-id>"
      return 2

end Lean4Conformance

def main (args : List String) : IO UInt32 :=
  Lean4Conformance.runMain args
