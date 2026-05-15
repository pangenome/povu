import PovuLean.VCF.Spec

/-!
Pure semantic VCF emitter.

The emitter maps certified `VariantCall`s to semantic VCF records.  It fixes
the povu-supported row constants (`QUAL`, `FILTER`, `FORMAT`) and computes the
INFO payload from the semantic call fields.  Converting these records to tabular
text, choosing the concrete `fileDate`, and decimal rendering of AF ratios are
outside the verified layer documented in `Spec.lean`.
-/

namespace PovuLean
namespace VCF

namespace VariantCall

def info (call : VariantCall) : Info :=
  { ac := call.alternateCounts
    af := call.alleleFrequencies
    an := call.totalAlleles
    ns := nsOfGenotypes call.genotypes
    atField := call.refTraversal :: call.alternateTraversals
    variantType := call.variantType
    tangled := call.tangled
    enclosingSite? := call.enclosingSite?
    level? := call.level? }

end VariantCall

/-- Emit one supported semantic VCF record from one verified variant call. -/
def recordOfCall (call : VariantCall) : Record :=
  { source := call.source
    chrom := call.chrom
    contigOrder := call.contigOrder
    pos := call.pos
    id := call.id
    ref := call.ref
    alternates := call.alternates
    qual := "60"
    filter := "PASS"
    info := call.info
    format := "GT"
    referenceAlleleCount := call.referenceAlleleCount
    refTraversal := call.refTraversal
    genotypes := call.genotypes }

def emitRecords (calls : List VariantCall) : List Record :=
  calls.map recordOfCall

@[simp] theorem recordOfCall_orderKey (call : VariantCall) :
    (recordOfCall call).orderKey = call.orderKey := rfl

theorem recordOfCall_infoConsistent {call : VariantCall}
    (hSource : call.source.SupportsVariantType call.variantType)
    (hSourceInfo : call.SourceInfoFields) :
    (recordOfCall call).InfoConsistent := by
  unfold Record.InfoConsistent recordOfCall VariantCall.info
  exact ⟨rfl, rfl, rfl, rfl, rfl, hSource, hSourceInfo⟩

theorem recordOfCall_wellFormed {call : VariantCall}
    (hCall : call.WellFormed) :
    (recordOfCall call).WellFormed := by
  rcases hCall with
    ⟨hChrom, hPos, hId, hRef, hRefTraversal, hAlternates,
      hAltWf, hAn, hSource, hSourceInfo, hGenotypes⟩
  unfold Record.WellFormed
  refine
    ⟨ hChrom
    , hPos
    , hId
    , hRef
    , hRefTraversal
    , hAlternates
    , ?_
    , hAn
    , rfl
    , rfl
    , rfl
    , ?_
    , hGenotypes ⟩
  · intro alt hMem
    exact hAltWf alt hMem
  · exact recordOfCall_infoConsistent hSource hSourceInfo

theorem orderedBy_map_of_key_eq {α β : Type}
    {keyA : α → OrderKey} {keyB : β → OrderKey} {f : α → β}
    (hKey : ∀ value, keyB (f value) = keyA value) :
    ∀ {values : List α},
      OrderedBy keyA values → OrderedBy keyB (values.map f)
  | [], _ => trivial
  | [_], _ => trivial
  | first :: second :: rest, hOrdered => by
      dsimp [OrderedBy] at hOrdered ⊢
      exact
        ⟨ by simpa [hKey first, hKey second] using hOrdered.1
        , orderedBy_map_of_key_eq (α := α) (β := β)
            (keyA := keyA) (keyB := keyB) (f := f) hKey hOrdered.2 ⟩

theorem emitRecords_ordered {calls : List VariantCall}
    (hOrdered : OrderedBy VariantCall.orderKey calls) :
    OrderedBy Record.orderKey (emitRecords calls) := by
  unfold emitRecords
  exact
    orderedBy_map_of_key_eq
      (keyA := VariantCall.orderKey)
      (keyB := Record.orderKey)
      (f := recordOfCall)
      (fun call => recordOfCall_orderKey call)
      hOrdered

end VCF
end PovuLean
