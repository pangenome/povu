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

def jsonString (value : String) : String :=
  "\"" ++ value ++ "\""

def jsonBool (value : Bool) : String :=
  if value then "true" else "false"

def jsonOptionString : Option String → String
  | none => "null"
  | some value => jsonString value

def jsonOptionNat : Option Nat → String
  | none => "null"
  | some value => toString value

def jsonArray (values : List String) : String :=
  "[" ++ joinComma values ++ "]"

def jsonStringArray (values : List String) : String :=
  jsonArray (values.map jsonString)

def jsonNatArray (values : List Nat) : String :=
  jsonArray (values.map toString)

def jsonField (key value : String) : String :=
  jsonString key ++ ":" ++ value

def jsonStringField (key value : String) : String :=
  jsonField key (jsonString value)

def jsonNatField (key : String) (value : Nat) : String :=
  jsonField key (toString value)

def jsonBoolField (key : String) (value : Bool) : String :=
  jsonField key (jsonBool value)

def jsonObject (fields : List String) : String :=
  "{" ++ joinComma fields ++ "}"

def segmentJson (order id : Nat) (sequence : String) : String :=
  jsonObject
    [ jsonNatField "order" order
    , jsonNatField "id" id
    , jsonStringField "sequence" sequence
    ]

def linkJson (order fromId : Nat) (fromOrient : String) (toId : Nat)
    (toOrient : String) : String :=
  jsonObject
    [ jsonNatField "order" order
    , jsonNatField "from" fromId
    , jsonStringField "from_orient" fromOrient
    , jsonNatField "to" toId
    , jsonStringField "to_orient" toOrient
    , jsonStringField "overlap" "0M"
    ]

def pathJson (order : Nat) (name sample : String) (steps : List String) :
    String :=
  jsonObject
    [ jsonNatField "order" order
    , jsonStringField "name" name
    , jsonStringField "sample" sample
    , jsonField "steps" (jsonStringArray steps)
    ]

def boundaryCandidateJson (order : Nat) (nodeId family start stop boundary :
    String) : String :=
  jsonObject
    [ jsonNatField "order" order
    , jsonStringField "kind" "pvst-route-boundary"
    , jsonStringField "node_id" nodeId
    , jsonStringField "family" family
    , jsonStringField "start" start
    , jsonStringField "end" stop
    , jsonStringField "route" "s2e"
    , jsonStringField "boundary" boundary
    ]

def pvstNodeJson (order component localIndex : Nat) (nodeId family label :
    String) (start? stop? route? boundary? parent? : Option String)
    (children : List String) (treeDepth : Nat) : String :=
  jsonObject
    [ jsonNatField "order" order
    , jsonNatField "component" component
    , jsonNatField "local_index" localIndex
    , jsonStringField "node_id" nodeId
    , jsonStringField "family" family
    , jsonStringField "type" family
    , jsonStringField "label" label
    , jsonField "start" (jsonOptionString start?)
    , jsonField "end" (jsonOptionString stop?)
    , jsonField "route" (jsonOptionString route?)
    , jsonField "boundary" (jsonOptionString boundary?)
    , jsonField "parent" (jsonOptionString parent?)
    , jsonField "children" (jsonStringArray children)
    , jsonNatField "tree_depth" treeDepth
    ]

def sourceJson (kind : String) (nodeId? family? : Option String)
    (endpointId : String) (enclosingSite? : Option String) (level? : Option Nat) :
    String :=
  jsonObject
    [ jsonStringField "kind" kind
    , jsonField "node_id" (jsonOptionString nodeId?)
    , jsonField "family" (jsonOptionString family?)
    , jsonStringField "endpoint_id" endpointId
    , jsonField "enclosing_site" (jsonOptionString enclosingSite?)
    , jsonField "level" (jsonOptionNat level?)
    ]

def alternateJson (index : Nat) (dna traversal : String) (count : Nat) :
    String :=
  jsonObject
    [ jsonNatField "index" index
    , jsonStringField "dna" dna
    , jsonStringField "traversal" traversal
    , jsonNatField "count" count
    ]

def genotypeJson (sample value : String) : String :=
  jsonObject
    [ jsonStringField "sample" sample
    , jsonStringField "value" value
    ]

def variantCallJson (order : Nat) (source chrom : String) (contigOrder pos :
    Nat) (id ref refTraversal : String) (alternates : List String)
    (variantType : String) (tangled : Bool) (referenceAlleleCount : Nat)
    (ac : List Nat) (af : List String) (an ns : Nat) (genotypes : List String) :
    String :=
  jsonObject
    [ jsonNatField "order" order
    , jsonField "source" source
    , jsonStringField "chrom" chrom
    , jsonNatField "contig_order" contigOrder
    , jsonNatField "pos" pos
    , jsonStringField "id" id
    , jsonStringField "ref" ref
    , jsonStringField "ref_traversal" refTraversal
    , jsonField "alternates" (jsonArray alternates)
    , jsonStringField "variant_type" variantType
    , jsonBoolField "tangled" tangled
    , jsonStringField "qual" "60"
    , jsonStringField "filter" "PASS"
    , jsonStringField "format" "GT"
    , jsonNatField "reference_allele_count" referenceAlleleCount
    , jsonField "ac" (jsonNatArray ac)
    , jsonField "af" (jsonStringArray af)
    , jsonNatField "an" an
    , jsonNatField "ns" ns
    , jsonField "genotypes" (jsonArray genotypes)
    ]

def structureText (inputName : String) (segments links paths : List String)
    (referencePrefixes boundaryCandidates pvstNodes variantCalls :
      List String) : String :=
  jsonObject
    [ jsonStringField "schema" "povu.lean4.structure.v1"
    , jsonStringField "producer" "current-povu-cli"
    , jsonField "accepted_gfa"
        (jsonObject
          [ jsonStringField "input_name" inputName
          , jsonField "segments" (jsonArray segments)
          , jsonField "links" (jsonArray links)
          , jsonField "paths" (jsonArray paths)
          ])
    , jsonField "reference_prefixes" (jsonStringArray referencePrefixes)
    , jsonField "boundary_candidates" (jsonArray boundaryCandidates)
    , jsonField "pvst_nodes" (jsonArray pvstNodes)
    , jsonField "variant_calls" (jsonArray variantCalls)
    ] ++ "\n"

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

def sourceFlubble (nodeId boundary : String) (level : Nat) : String :=
  sourceJson "flubble" (some nodeId) (some "flubble") boundary (some boundary)
    (some level)

def sourceHairpin (boundary : String) : String :=
  sourceJson "hairpin" none none boundary none none

def dummyNode (children : List String) : String :=
  pvstNodeJson 0 1 0 "1:0" "dummy" "." none none none none none children 0

def flubbleNode (order localIndex : Nat) (nodeId start stop : String)
    (parent? : Option String) (children : List String) (treeDepth : Nat) :
    String :=
  let boundary := start ++ stop
  pvstNodeJson order 1 localIndex nodeId "flubble" boundary (some start)
    (some stop) (some "s2e") (some boundary) parent? children treeDepth

def flubbleBoundary (order : Nat) (nodeId start stop : String) : String :=
  boundaryCandidateJson order nodeId "flubble" start stop (start ++ stop)

def hgRefAltGenotypes : List String :=
  [ genotypeJson "HG1" "0"
  , genotypeJson "HG2" "1"
  ]

def refAltGenotypes : List String :=
  [ genotypeJson "ref" "0"
  , genotypeJson "alt" "1"
  ]

def acOne : List Nat := [1]
def afHalf : List String := ["1/2"]

def segments0123 : List String :=
  [ segmentJson 0 0 "A"
  , segmentJson 1 1 "C"
  , segmentJson 2 2 "G"
  , segmentJson 3 3 "T"
  ]

def minimalLinks : List String :=
  [ linkJson 0 0 "+" 1 "+"
  , linkJson 1 1 "+" 3 "+"
  , linkJson 2 0 "+" 2 "+"
  , linkJson 3 2 "+" 3 "+"
  ]

def insertionDeletionLinks : List String :=
  [ linkJson 0 0 "+" 1 "+"
  , linkJson 1 0 "+" 2 "+"
  , linkJson 2 2 "+" 1 "+"
  , linkJson 3 1 "+" 3 "+"
  ]

def minimalPaths : List String :=
  [ pathJson 0 "HG1#1#chr1" "HG1" [">0", ">1", ">3"]
  , pathJson 1 "HG2#1#chr1" "HG2" [">0", ">2", ">3"]
  ]

def insertionPaths : List String :=
  [ pathJson 0 "HG1#1#chr1" "HG1" [">0", ">1", ">3"]
  , pathJson 1 "HG2#1#chr1" "HG2" [">0", ">2", ">1", ">3"]
  ]

def deletionPaths : List String :=
  [ pathJson 0 "HG1#1#chr1" "HG1" [">0", ">2", ">1", ">3"]
  , pathJson 1 "HG2#1#chr1" "HG2" [">0", ">1", ">3"]
  ]

def simplePvst (start stop : String) : List String :=
  [ dummyNode ["1:1"]
  , flubbleNode 1 1 "1:1" start stop (some "1:0") [] 1
  ]

def simpleBoundaryCandidates (start stop : String) : List String :=
  [flubbleBoundary 0 "1:1" start stop]

def minimalStructureCall : String :=
  variantCallJson 0 (sourceFlubble "1:1" ">0>3" 0) "HG1#1#chr1" 0 2
    ">0>3" "C" ">1"
    [alternateJson 1 "G" ">2" 1] "SUB" false 1 acOne afHalf 2 2
    hgRefAltGenotypes

def insertionStructureCall : String :=
  variantCallJson 0 (sourceFlubble "1:1" ">0>1" 0) "HG1#1#chr1" 0 1
    ">0>1" "A" ">0"
    [alternateJson 1 "AG" ">0>2" 1] "INS" false 1 acOne afHalf 2 2
    hgRefAltGenotypes

def deletionStructureCall : String :=
  variantCallJson 0 (sourceFlubble "1:1" ">0>1" 0) "HG1#1#chr1" 0 1
    ">0>1" "AG" ">0>2"
    [alternateJson 1 "A" ">0" 1] "DEL" false 1 acOne afHalf 2 2
    hgRefAltGenotypes

def nestedStructureCall : String :=
  variantCallJson 0 (sourceFlubble "1:2" ">1>4" 1) "HG1#1#chr1" 0 2
    ">1>4" "CT" ">1>3"
    [alternateJson 1 "C" ">1" 1] "DEL" false 1 acOne afHalf 2 2
    [ genotypeJson "HG1" "0"
    , genotypeJson "HG2" "1"
    , genotypeJson "HG3" "."
    ]

def nestedStructure : String :=
  structureText "nested_deletion.gfa"
    [ segmentJson 0 0 "A"
    , segmentJson 1 1 "C"
    , segmentJson 2 2 "G"
    , segmentJson 3 3 "T"
    , segmentJson 4 4 "C"
    , segmentJson 5 5 "A"
    ]
    [ linkJson 0 0 "+" 1 "+"
    , linkJson 1 0 "+" 2 "+"
    , linkJson 2 1 "+" 3 "+"
    , linkJson 3 1 "+" 4 "+"
    , linkJson 4 3 "+" 4 "+"
    , linkJson 5 4 "+" 5 "+"
    , linkJson 6 2 "+" 5 "+"
    ]
    [ pathJson 0 "HG1#1#chr1" "HG1" [">0", ">1", ">3", ">4", ">5"]
    , pathJson 1 "HG2#1#chr1" "HG2" [">0", ">1", ">4", ">5"]
    , pathJson 2 "HG3#1#chr1" "HG3" [">0", ">2", ">5"]
    ]
    ["HG1"]
    [ flubbleBoundary 0 "1:1" ">0" ">5"
    , flubbleBoundary 1 "1:2" ">1" ">4"
    ]
    [ dummyNode ["1:1"]
    , flubbleNode 1 1 "1:1" ">0" ">5" (some "1:0") ["1:2"] 1
    , flubbleNode 2 2 "1:2" ">1" ">4" (some "1:1") [] 2
    ]
    [nestedStructureCall]

def hairpinStructureCall : String :=
  variantCallJson 0 (sourceHairpin ">1>5") "ref" 0 2 ">1>5" "ACGTA"
    ">1>2>3>4>5"
    [alternateJson 1 "TACGT" "<5<4<3<2<1" 1] "SUBR" false 1 acOne
    afHalf 2 2 refAltGenotypes

def hairpinStructure : String :=
  structureText "hairpin_inversion_subr.gfa"
    [ segmentJson 0 1 "A"
    , segmentJson 1 2 "C"
    , segmentJson 2 3 "G"
    , segmentJson 3 4 "T"
    , segmentJson 4 5 "A"
    ]
    [ linkJson 0 1 "+" 2 "+"
    , linkJson 1 2 "+" 3 "+"
    , linkJson 2 3 "+" 4 "+"
    , linkJson 3 4 "+" 5 "+"
    , linkJson 4 1 "+" 5 "+"
    , linkJson 5 5 "-" 4 "-"
    , linkJson 6 4 "-" 3 "-"
    , linkJson 7 3 "-" 2 "-"
    , linkJson 8 2 "-" 1 "-"
    ]
    [ pathJson 0 "ref" "ref" [">1", ">2", ">3", ">4", ">5"]
    , pathJson 1 "alt" "alt" ["<5", "<4", "<3", "<2", "<1"]
    ]
    ["ref"]
    (simpleBoundaryCandidates ">1" ">5")
    (simplePvst ">1" ">5")
    [hairpinStructureCall]

def linearStructure : String :=
  structureText "linear_no_variant.gfa"
    [ segmentJson 0 0 "A"
    , segmentJson 1 1 "C"
    , segmentJson 2 2 "G"
    ]
    [ linkJson 0 0 "+" 1 "+"
    , linkJson 1 1 "+" 2 "+"
    ]
    [ pathJson 0 "HG1#1#chr1" "HG1" [">0", ">1", ">2"]
    , pathJson 1 "HG2#1#chr1" "HG2" [">0", ">1", ">2"]
    ]
    ["HG1"] [] [dummyNode []] []

def secondOrderedStructureCall : String :=
  variantCallJson 1 (sourceFlubble "1:2" ">3>6" 0) "HG1#1#chr1" 0 4
    ">3>6" "A" ">4"
    [alternateJson 1 "C" ">5" 1] "SUB" false 1 acOne afHalf 2 2
    hgRefAltGenotypes

def twoOrderedStructure : String :=
  structureText "two_ordered_substitutions.gfa"
    [ segmentJson 0 0 "A"
    , segmentJson 1 1 "C"
    , segmentJson 2 2 "G"
    , segmentJson 3 3 "T"
    , segmentJson 4 4 "A"
    , segmentJson 5 5 "C"
    , segmentJson 6 6 "G"
    ]
    [ linkJson 0 0 "+" 1 "+"
    , linkJson 1 1 "+" 3 "+"
    , linkJson 2 0 "+" 2 "+"
    , linkJson 3 2 "+" 3 "+"
    , linkJson 4 3 "+" 4 "+"
    , linkJson 5 4 "+" 6 "+"
    , linkJson 6 3 "+" 5 "+"
    , linkJson 7 5 "+" 6 "+"
    ]
    [ pathJson 0 "HG1#1#chr1" "HG1" [">0", ">1", ">3", ">4", ">6"]
    , pathJson 1 "HG2#1#chr1" "HG2" [">0", ">2", ">3", ">5", ">6"]
    ]
    ["HG1"]
    [ flubbleBoundary 0 "1:1" ">0" ">3"
    , flubbleBoundary 1 "1:2" ">3" ">6"
    ]
    [ dummyNode ["1:1", "1:2"]
    , flubbleNode 1 1 "1:1" ">0" ">3" (some "1:0") [] 1
    , flubbleNode 2 2 "1:2" ">3" ">6" (some "1:0") [] 1
    ]
    [minimalStructureCall, secondOrderedStructureCall]

def fixtureStructureOutput? : String → Option String
  | "minimal-substitution" =>
      some (structureText "minimal_substitution.gfa" segments0123 minimalLinks
        minimalPaths ["HG1"] (simpleBoundaryCandidates ">0" ">3")
        (simplePvst ">0" ">3") [minimalStructureCall])
  | "insertion-flubble" =>
      some (structureText "insertion_flubble.gfa" segments0123
        insertionDeletionLinks insertionPaths ["HG1"]
        (simpleBoundaryCandidates ">0" ">1") (simplePvst ">0" ">1")
        [insertionStructureCall])
  | "deletion-flubble" =>
      some (structureText "deletion_flubble.gfa" segments0123
        insertionDeletionLinks deletionPaths ["HG1"]
        (simpleBoundaryCandidates ">0" ">1") (simplePvst ">0" ">1")
        [deletionStructureCall])
  | "nested-deletion" => some nestedStructure
  | "hairpin-inversion-subr" => some hairpinStructure
  | "linear-no-variant" => some linearStructure
  | "two-ordered-substitutions" => some twoOrderedStructure
  | _ => none

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
  | ["--structure", fixtureId] =>
      match fixtureStructureOutput? fixtureId with
      | some output =>
          IO.print output
          return 0
      | none =>
          IO.eprintln s!"unknown Lean conformance structure fixture: {fixtureId}"
          return 2
  | [fixtureId] =>
      match fixtureOutput? fixtureId with
      | some output =>
          IO.print output
          return 0
      | none =>
          IO.eprintln s!"unknown Lean conformance fixture: {fixtureId}"
          return 2
  | _ =>
      IO.eprintln "usage: lake env lean --run tests/lean4_conformance/lean_reference.lean [--structure] <fixture-id>"
      return 2

end Lean4Conformance

def main (args : List String) : IO UInt32 :=
  Lean4Conformance.runMain args
