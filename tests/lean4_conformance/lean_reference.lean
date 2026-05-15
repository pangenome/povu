import PovuLean.Pipeline
import PovuLean.VCF.Emit
import PovuLean.Algorithms.FlubbleTree.Examples

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

def sampleHG1Ref : VCF.GenotypeColumn :=
  { sample := "HG1"
    alleles := [VCF.GenotypeAllele.ref] }

def sampleHG2Alt : VCF.GenotypeColumn :=
  { sample := "HG2"
    alleles := [VCF.GenotypeAllele.alt 1] }

def minimalSubstitutionNode : Algorithms.FlubbleTree.Node :=
  { boundary := Algorithms.Flubble.Examples.simpleBoundary
    parent? := none }

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

def fixtureOutput? : String → Option String
  | "minimal-substitution" =>
      some (vcfText ["HG1", "HG2"] (VCF.emitRecords [minimalSubstitutionCall]))
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
