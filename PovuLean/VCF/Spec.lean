import PovuLean.Algorithms.FlubbleTree.Correctness
import PovuLean.Algorithms.Hairpin.Correctness

/-!
Semantic VCF subset used by the Lean povu reference pipeline.

This module models records after graph variants have already been verified and
after a caller has supplied semantic allele spellings.  It intentionally does
not prove byte-level VCF serialization, decimal formatting for `AF`, PanSN name
parsing, or Rust/Lean conformance for extracting DNA strings from graph slices.
Those are downstream conformance obligations.  The trusted boundary here is a
`VariantCall`: it carries a verified flubble-tree or hairpin source together
with constructed REF/ALT alleles, traversal strings, allele counts, and genotype
columns.

Supported VCF subset:

* VCF 4.2-style records with `CHROM`, one-based positive `POS`, nonempty `ID`,
  nonempty `REF`, and at least one nonempty `ALT`.
* Fixed `QUAL = 60`, `FILTER = PASS`, and `FORMAT = GT`.
* Variant types `DEL`, `INS`, `SUB`, and `SUBR`.  Flubble-tree sources may emit
  `DEL`, `INS`, or `SUB`; hairpin sources emit `SUBR`.
* INFO fields `AC`, semantic `AF` ratios, `AN`, `NS`, `AT`, `VARTYPE`,
  `TANGLED`, and optional `ES`/`LV`.  Non-`SUBR` records carry `ES` and `LV`;
  `SUBR` hairpin inversion records omit them.
* Genotypes are phased allele calls where `0` is REF, alternate calls are
  one-based, and missing phases are represented semantically by `missing`.

Allele construction is explicit.  Insertions and deletions retain an anchor so
neither serialized allele is empty; substitutions and reverse substitutions use
the supplied reference and alternate sequence directly.
-/

namespace PovuLean
namespace VCF

open Core

abbrev ChromName := String
abbrev RecordId := String
abbrev TraversalString := String
abbrev SampleName := String

/-- Variant classes emitted by the supported povu VCF subset. -/
inductive VariantType where
  | del
  | ins
  | sub
  | subr
  deriving Repr, DecidableEq, BEq

namespace VariantType

def vcfLabel : VariantType → String
  | del => "DEL"
  | ins => "INS"
  | sub => "SUB"
  | subr => "SUBR"

def IsFlubbleType : VariantType → Prop
  | del => True
  | ins => True
  | sub => True
  | subr => False

instance decidableIsFlubbleType (ty : VariantType) :
    Decidable ty.IsFlubbleType := by
  cases ty
  · exact isTrue trivial
  · exact isTrue trivial
  · exact isTrue trivial
  · exact isFalse (by intro h; exact h)

end VariantType

/-- Construction certificate for one alternate allele. -/
inductive AlleleConstruction where
  | deletion (anchor deleted : String)
  | insertion (anchor inserted : String)
  | substitution (reference alternate : String)
  | reverseSubstitution (reference alternate : String)
  deriving Repr, DecidableEq, BEq

namespace AlleleConstruction

def ref : AlleleConstruction → String
  | deletion anchor deleted => anchor ++ deleted
  | insertion anchor _inserted => anchor
  | substitution reference _alternate => reference
  | reverseSubstitution reference _alternate => reference

def alt : AlleleConstruction → String
  | deletion anchor _deleted => anchor
  | insertion anchor inserted => anchor ++ inserted
  | substitution _reference alternate => alternate
  | reverseSubstitution _reference alternate => alternate

def variantType : AlleleConstruction → VariantType
  | deletion _ _ => VariantType.del
  | insertion _ _ => VariantType.ins
  | substitution _ _ => VariantType.sub
  | reverseSubstitution _ _ => VariantType.subr

/--
The construction is usable in VCF: both emitted alleles are nonempty, and
anchored events record nonempty anchor and changed sequence components.
-/
def WellFormed (construction : AlleleConstruction) : Prop :=
  construction.ref ≠ ""
    ∧ construction.alt ≠ ""
    ∧
      match construction with
      | deletion anchor deleted => anchor ≠ "" ∧ deleted ≠ ""
      | insertion anchor inserted => anchor ≠ "" ∧ inserted ≠ ""
      | substitution reference alternate =>
          reference ≠ "" ∧ alternate ≠ ""
      | reverseSubstitution reference alternate =>
          reference ≠ "" ∧ alternate ≠ ""

instance decidableWellFormed (construction : AlleleConstruction) :
    Decidable construction.WellFormed := by
  cases construction <;> unfold WellFormed ref alt <;> infer_instance

end AlleleConstruction

/-- One alternate allele group and its supporting graph traversal. -/
structure AlternateAllele where
  construction : AlleleConstruction
  traversal : TraversalString
  count : Nat
  deriving Repr, DecidableEq, BEq

namespace AlternateAllele

def ref (alt : AlternateAllele) : String :=
  alt.construction.ref

def allele (alt : AlternateAllele) : String :=
  alt.construction.alt

def variantType (alt : AlternateAllele) : VariantType :=
  alt.construction.variantType

def WellFormed (alt : AlternateAllele) : Prop :=
  alt.construction.WellFormed ∧ alt.traversal ≠ ""

instance decidableWellFormed (alt : AlternateAllele) :
    Decidable alt.WellFormed :=
  inferInstanceAs
    (Decidable (alt.construction.WellFormed ∧ alt.traversal ≠ ""))

end AlternateAllele

/-- Source variant carried through to VCF semantics. -/
inductive VariantSource where
  | flubble (node : Algorithms.FlubbleTree.Node)
  | hairpin (boundary : Algorithms.Hairpin.Boundary)
  deriving Repr, DecidableEq

namespace VariantSource

def primaryKey : VariantSource → Nat
  | flubble node => node.boundary.openEdge.id
  | hairpin boundary => boundary.outerEdge.id

def secondaryKey : VariantSource → Nat
  | flubble node => node.boundary.closeEdge.id
  | hairpin boundary => boundary.innerEdge.id

def SupportsVariantType (source : VariantSource) (ty : VariantType) : Prop :=
  match source with
  | flubble _ => ty.IsFlubbleType
  | hairpin _ => ty = VariantType.subr

instance decidableSupportsVariantType (source : VariantSource)
    (ty : VariantType) :
    Decidable (source.SupportsVariantType ty) := by
  cases source <;> cases ty
  · exact isTrue trivial
  · exact isTrue trivial
  · exact isTrue trivial
  · exact isFalse (by intro h; exact h)
  · exact isFalse (by intro h; cases h)
  · exact isFalse (by intro h; cases h)
  · exact isFalse (by intro h; cases h)
  · exact isTrue rfl

end VariantSource

/-- One phased genotype allele in a sample column. -/
inductive GenotypeAllele where
  | missing
  | ref
  | alt (index : Nat)
  deriving Repr, DecidableEq, BEq

namespace GenotypeAllele

def isCalled : GenotypeAllele → Bool
  | missing => false
  | ref => true
  | alt _ => true

def Valid (altCount : Nat) : GenotypeAllele → Prop
  | missing => True
  | ref => True
  | alt index => 1 ≤ index ∧ index ≤ altCount

instance decidableValid (altCount : Nat) (allele : GenotypeAllele) :
    Decidable (allele.Valid altCount) := by
  cases allele
  · exact isTrue trivial
  · exact isTrue trivial
  · unfold Valid
    infer_instance

end GenotypeAllele

/-- One VCF sample column with phased allele entries. -/
structure GenotypeColumn where
  sample : SampleName
  alleles : List GenotypeAllele
  deriving Repr, DecidableEq, BEq

namespace GenotypeColumn

def hasData (column : GenotypeColumn) : Bool :=
  column.alleles.any GenotypeAllele.isCalled

def WellFormed (altCount : Nat) (column : GenotypeColumn) : Prop :=
  column.sample ≠ ""
    ∧ column.alleles ≠ []
    ∧ ∀ allele, allele ∈ column.alleles → allele.Valid altCount

instance decidableWellFormed (altCount : Nat) (column : GenotypeColumn) :
    Decidable (column.WellFormed altCount) :=
  inferInstanceAs
    (Decidable
      ( column.sample ≠ ""
        ∧ column.alleles ≠ []
        ∧ ∀ allele, allele ∈ column.alleles → allele.Valid altCount))

end GenotypeColumn

def genotypeColumnsWellFormed
    (altCount : Nat) (columns : List GenotypeColumn) : Prop :=
  ∀ column, column ∈ columns → column.WellFormed altCount

instance decidableGenotypeColumnsWellFormed
    (altCount : Nat) (columns : List GenotypeColumn) :
    Decidable (genotypeColumnsWellFormed altCount columns) :=
  inferInstanceAs
    (Decidable
      (∀ column, column ∈ columns → column.WellFormed altCount))

def nsOfGenotypes (columns : List GenotypeColumn) : Nat :=
  (columns.filter GenotypeColumn.hasData).length

/-- Semantic AF ratio before the intentionally trusted decimal formatting step. -/
structure AlleleFrequency where
  alternateCount : Nat
  totalAlleles : Nat
  deriving Repr, DecidableEq, BEq

def sumNat : List Nat → Nat
  | [] => 0
  | value :: rest => value + sumNat rest

/-- INFO fields in the supported semantic subset. -/
structure Info where
  ac : List Nat
  af : List AlleleFrequency
  an : Nat
  ns : Nat
  atField : List TraversalString
  variantType : VariantType
  tangled : Bool
  enclosingSite? : Option String
  level? : Option Nat
  deriving Repr, DecidableEq, BEq

/-- Deterministic semantic record ordering key. -/
structure OrderKey where
  contigOrder : Nat
  pos : Nat
  sourcePrimary : Nat
  sourceSecondary : Nat
  deriving Repr, DecidableEq, BEq

namespace OrderKey

def le (left right : OrderKey) : Prop :=
  left.contigOrder < right.contigOrder
    ∨ (left.contigOrder = right.contigOrder
      ∧ (left.pos < right.pos
        ∨ (left.pos = right.pos
          ∧ (left.sourcePrimary < right.sourcePrimary
            ∨ (left.sourcePrimary = right.sourcePrimary
              ∧ left.sourceSecondary ≤ right.sourceSecondary)))))

instance decidableLe (left right : OrderKey) : Decidable (le left right) := by
  unfold le
  infer_instance

end OrderKey

def OrderedBy (key : α → OrderKey) : List α → Prop
  | [] => True
  | [_] => True
  | first :: second :: rest =>
      OrderKey.le (key first) (key second) ∧ OrderedBy key (second :: rest)

def orderedByDecidable (key : α → OrderKey) :
    (values : List α) → Decidable (OrderedBy key values)
  | [] => isTrue trivial
  | [_] => isTrue trivial
  | first :: second :: rest =>
      match inferInstanceAs (Decidable (OrderKey.le (key first) (key second))),
          orderedByDecidable key (second :: rest) with
      | isTrue hLe, isTrue hRest => isTrue ⟨hLe, hRest⟩
      | isFalse hLe, _ => isFalse (fun h => hLe h.1)
      | _, isFalse hRest => isFalse (fun h => hRest h.2)

instance decidableOrderedBy (key : α → OrderKey) (values : List α) :
    Decidable (OrderedBy key values) :=
  orderedByDecidable key values

/-- INFO `ES`/`LV` presence policy for the supported subset. -/
def InfoFieldsForType
    (ty : VariantType) (enclosingSite? : Option String) (level? : Option Nat) :
    Prop :=
  (ty = VariantType.subr ∧ enclosingSite? = none ∧ level? = none)
    ∨ (ty.IsFlubbleType
      ∧ enclosingSite?.isSome = true
      ∧ level?.isSome = true)

instance decidableInfoFieldsForType
    (ty : VariantType) (enclosingSite? : Option String) (level? : Option Nat) :
    Decidable (InfoFieldsForType ty enclosingSite? level?) := by
  unfold InfoFieldsForType
  infer_instance

/--
Semantic call input consumed by the Lean VCF emitter.  A call is a verified
variant source plus all REF/ALT, traversal, count, and genotype facts needed to
produce the supported semantic VCF record.
-/
structure VariantCall where
  source : VariantSource
  chrom : ChromName
  contigOrder : Nat
  pos : Nat
  id : RecordId
  ref : String
  refTraversal : TraversalString
  alternates : List AlternateAllele
  variantType : VariantType
  tangled : Bool
  enclosingSite? : Option String
  level? : Option Nat
  referenceAlleleCount : Nat
  genotypes : List GenotypeColumn
  deriving Repr, DecidableEq

namespace VariantCall

def alternateAlleles (call : VariantCall) : List String :=
  call.alternates.map AlternateAllele.allele

def alternateCounts (call : VariantCall) : List Nat :=
  call.alternates.map AlternateAllele.count

def alternateTraversals (call : VariantCall) : List TraversalString :=
  call.alternates.map AlternateAllele.traversal

def alleleNumberCount (call : VariantCall) : Nat :=
  call.alternates.length

def totalAlleles (call : VariantCall) : Nat :=
  call.referenceAlleleCount + sumNat call.alternateCounts

def alleleFrequencies (call : VariantCall) : List AlleleFrequency :=
  call.alternateCounts.map
    (fun count =>
      { alternateCount := count
        totalAlleles := call.totalAlleles })

def orderKey (call : VariantCall) : OrderKey :=
  { contigOrder := call.contigOrder
    pos := call.pos
    sourcePrimary := call.source.primaryKey
    sourceSecondary := call.source.secondaryKey }

def SourceInfoFields (call : VariantCall) : Prop :=
  InfoFieldsForType call.variantType call.enclosingSite? call.level?

def AlternatesWellFormed (call : VariantCall) : Prop :=
  ∀ alt, alt ∈ call.alternates →
    alt.WellFormed
      ∧ alt.ref = call.ref
      ∧ alt.variantType = call.variantType

def WellFormed (call : VariantCall) : Prop :=
  call.chrom ≠ ""
    ∧ call.pos > 0
    ∧ call.id ≠ ""
    ∧ call.ref ≠ ""
    ∧ call.refTraversal ≠ ""
    ∧ call.alternates ≠ []
    ∧ call.AlternatesWellFormed
    ∧ call.totalAlleles > 0
    ∧ call.source.SupportsVariantType call.variantType
    ∧ call.SourceInfoFields
    ∧ genotypeColumnsWellFormed call.alleleNumberCount call.genotypes

instance decidableSourceInfoFields (call : VariantCall) :
    Decidable call.SourceInfoFields := by
  unfold SourceInfoFields
  infer_instance

instance decidableAlternatesWellFormed (call : VariantCall) :
    Decidable call.AlternatesWellFormed :=
  inferInstanceAs
    (Decidable
      (∀ alt, alt ∈ call.alternates →
        alt.WellFormed
          ∧ alt.ref = call.ref
          ∧ alt.variantType = call.variantType))

instance decidableWellFormed (call : VariantCall) :
    Decidable call.WellFormed :=
  inferInstanceAs
    (Decidable
      ( call.chrom ≠ ""
        ∧ call.pos > 0
        ∧ call.id ≠ ""
        ∧ call.ref ≠ ""
        ∧ call.refTraversal ≠ ""
        ∧ call.alternates ≠ []
        ∧ call.AlternatesWellFormed
        ∧ call.totalAlleles > 0
        ∧ call.source.SupportsVariantType call.variantType
        ∧ call.SourceInfoFields
        ∧ genotypeColumnsWellFormed call.alleleNumberCount call.genotypes))

end VariantCall

/-- A semantic VCF record with hidden proof-oriented fields retained. -/
structure Record where
  source : VariantSource
  chrom : ChromName
  contigOrder : Nat
  pos : Nat
  id : RecordId
  ref : String
  alternates : List AlternateAllele
  qual : String
  filter : String
  info : Info
  format : String
  referenceAlleleCount : Nat
  refTraversal : TraversalString
  genotypes : List GenotypeColumn
  deriving Repr, DecidableEq

namespace Record

def alternateAlleles (record : Record) : List String :=
  record.alternates.map AlternateAllele.allele

def alternateCounts (record : Record) : List Nat :=
  record.alternates.map AlternateAllele.count

def alternateTraversals (record : Record) : List TraversalString :=
  record.alternates.map AlternateAllele.traversal

def alleleNumberCount (record : Record) : Nat :=
  record.alternates.length

def totalAlleles (record : Record) : Nat :=
  record.referenceAlleleCount + sumNat record.alternateCounts

def alleleFrequencies (record : Record) : List AlleleFrequency :=
  record.alternateCounts.map
    (fun count =>
      { alternateCount := count
        totalAlleles := record.totalAlleles })

def orderKey (record : Record) : OrderKey :=
  { contigOrder := record.contigOrder
    pos := record.pos
    sourcePrimary := record.source.primaryKey
    sourceSecondary := record.source.secondaryKey }

def InfoConsistent (record : Record) : Prop :=
  record.info.ac = record.alternateCounts
    ∧ record.info.af = record.alleleFrequencies
    ∧ record.info.an = record.totalAlleles
    ∧ record.info.ns = nsOfGenotypes record.genotypes
    ∧ record.info.atField = record.refTraversal :: record.alternateTraversals
    ∧ record.source.SupportsVariantType record.info.variantType
    ∧ InfoFieldsForType
        record.info.variantType
        record.info.enclosingSite?
        record.info.level?

def AlternatesWellFormed (record : Record) : Prop :=
  ∀ alt, alt ∈ record.alternates →
    alt.WellFormed
      ∧ alt.ref = record.ref
      ∧ alt.variantType = record.info.variantType

def WellFormed (record : Record) : Prop :=
  record.chrom ≠ ""
    ∧ record.pos > 0
    ∧ record.id ≠ ""
    ∧ record.ref ≠ ""
    ∧ record.refTraversal ≠ ""
    ∧ record.alternates ≠ []
    ∧ record.AlternatesWellFormed
    ∧ record.totalAlleles > 0
    ∧ record.qual = "60"
    ∧ record.filter = "PASS"
    ∧ record.format = "GT"
    ∧ record.InfoConsistent
    ∧ genotypeColumnsWellFormed record.alleleNumberCount record.genotypes

instance decidableInfoConsistent (record : Record) :
    Decidable record.InfoConsistent := by
  unfold InfoConsistent
  infer_instance

instance decidableAlternatesWellFormed (record : Record) :
    Decidable record.AlternatesWellFormed :=
  inferInstanceAs
    (Decidable
      (∀ alt, alt ∈ record.alternates →
        alt.WellFormed
          ∧ alt.ref = record.ref
          ∧ alt.variantType = record.info.variantType))

instance decidableWellFormed (record : Record) :
    Decidable record.WellFormed :=
  inferInstanceAs
    (Decidable
      ( record.chrom ≠ ""
        ∧ record.pos > 0
        ∧ record.id ≠ ""
        ∧ record.ref ≠ ""
        ∧ record.refTraversal ≠ ""
        ∧ record.alternates ≠ []
        ∧ record.AlternatesWellFormed
        ∧ record.totalAlleles > 0
        ∧ record.qual = "60"
        ∧ record.filter = "PASS"
        ∧ record.format = "GT"
        ∧ record.InfoConsistent
        ∧ genotypeColumnsWellFormed record.alleleNumberCount record.genotypes))

end Record

end VCF
end PovuLean
