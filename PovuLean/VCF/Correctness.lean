import PovuLean.VCF.Emit

/-!
Correctness of semantic VCF emission from verified graph variants.

The main theorem, `emitVcfRecords_correct_for_gfa`, composes the accepted-GFA
input contract, verified flubble-tree construction, verified hairpin detection,
and semantic VCF emitter.  It says that every emitted record is well formed in
the supported VCF subset, is ordered deterministically, and is faithful to a
verified source reported by the Lean reference pipeline.
-/

namespace PovuLean
namespace VCF

open Core

namespace VariantSource

/-- The source was reported by the Lean reference pipeline before VCF emission. -/
def ReportedByPipeline {g : Graph}
    (frame : TraversalFrame g)
    (scan : Algorithms.Hairpin.HairpinScanAssignment g)
    (hierarchy : Algorithms.FlubbleTree.Hierarchy) :
    VariantSource → Prop
  | flubble node => node ∈ hierarchy.nodes
  | hairpin boundary => boundary ∈ Algorithms.Hairpin.detectHairpins frame scan

/-- The reported source satisfies the verified graph-variant specification. -/
def Derived {g : Graph}
    (frame : TraversalFrame g)
    (classes : Algorithms.Flubble.CycleClassAssignment g)
    (scan : Algorithms.Hairpin.HairpinScanAssignment g)
    (hierarchy : Algorithms.FlubbleTree.Hierarchy) :
    VariantSource → Prop
  | flubble node =>
      node ∈ hierarchy.nodes
        ∧ Algorithms.Flubble.IsFlubbleBoundary frame classes node.boundary
  | hairpin boundary =>
      boundary ∈ Algorithms.Hairpin.detectHairpins frame scan
        ∧ Algorithms.Hairpin.IsPaperHairpinBoundary frame boundary

end VariantSource

/--
Calls supplied to the emitter are semantically valid, pipeline-reported, and in
deterministic order.
-/
structure ReferenceCallSet {g : Graph}
    (frame : TraversalFrame g)
    (scan : Algorithms.Hairpin.HairpinScanAssignment g)
    (hierarchy : Algorithms.FlubbleTree.Hierarchy)
    (calls : List VariantCall) : Prop where
  wellFormed : ∀ call, call ∈ calls → call.WellFormed
  reported :
    ∀ call, call ∈ calls →
      call.source.ReportedByPipeline frame scan hierarchy
  ordered : OrderedBy VariantCall.orderKey calls

/-- Semantic correctness contract for emitted records. -/
structure EmissionCorrect {g : Graph}
    (frame : TraversalFrame g)
    (classes : Algorithms.Flubble.CycleClassAssignment g)
    (scan : Algorithms.Hairpin.HairpinScanAssignment g)
    (hierarchy : Algorithms.FlubbleTree.Hierarchy)
    (calls : List VariantCall)
    (records : List Record) : Prop where
  emitted : records = emitRecords calls
  wellFormed : ∀ record, record ∈ records → record.WellFormed
  derived :
    ∀ record, record ∈ records →
      record.source.Derived frame classes scan hierarchy
  ordered : OrderedBy Record.orderKey records

theorem flubbleSource_derived_of_correct {g : Graph}
    {frame : TraversalFrame g}
    {classes : Algorithms.Flubble.CycleClassAssignment g}
    {scan : Algorithms.Hairpin.HairpinScanAssignment g}
    {hierarchy : Algorithms.FlubbleTree.Hierarchy}
    {node : Algorithms.FlubbleTree.Node}
    (hHierarchy : Algorithms.FlubbleTree.IsCorrectHierarchy frame classes hierarchy)
    (hNode : node ∈ hierarchy.nodes) :
    (VariantSource.flubble node).Derived frame classes scan hierarchy := by
  unfold VariantSource.Derived
  constructor
  · exact hNode
  · have hBoundaryMem : node.boundary ∈ hierarchy.boundaries := by
      unfold Algorithms.FlubbleTree.Hierarchy.boundaries
      exact List.mem_map.mpr ⟨node, hNode, rfl⟩
    exact (hHierarchy.2.1 node.boundary).mp hBoundaryMem

theorem hairpinSource_derived_of_sound {g : Graph}
    {frame : TraversalFrame g}
    {classes : Algorithms.Flubble.CycleClassAssignment g}
    {scan : Algorithms.Hairpin.HairpinScanAssignment g}
    {hierarchy : Algorithms.FlubbleTree.Hierarchy}
    {boundary : Algorithms.Hairpin.Boundary}
    (hSound :
      ∀ boundary,
        boundary ∈ Algorithms.Hairpin.detectHairpins frame scan →
          Algorithms.Hairpin.IsPaperHairpinBoundary frame boundary)
    (hBoundary : boundary ∈ Algorithms.Hairpin.detectHairpins frame scan) :
    (VariantSource.hairpin boundary).Derived frame classes scan hierarchy := by
  unfold VariantSource.Derived
  exact ⟨hBoundary, hSound boundary hBoundary⟩

theorem reportedSource_derived_of_correct {g : Graph}
    {frame : TraversalFrame g}
    {classes : Algorithms.Flubble.CycleClassAssignment g}
    {scan : Algorithms.Hairpin.HairpinScanAssignment g}
    {hierarchy : Algorithms.FlubbleTree.Hierarchy}
    {source : VariantSource}
    (hHierarchy : Algorithms.FlubbleTree.IsCorrectHierarchy frame classes hierarchy)
    (hHairpinSound :
      ∀ boundary,
        boundary ∈ Algorithms.Hairpin.detectHairpins frame scan →
          Algorithms.Hairpin.IsPaperHairpinBoundary frame boundary)
    (hReported : source.ReportedByPipeline frame scan hierarchy) :
    source.Derived frame classes scan hierarchy := by
  cases source with
  | flubble node =>
      exact flubbleSource_derived_of_correct hHierarchy hReported
  | hairpin boundary =>
      exact hairpinSource_derived_of_sound hHairpinSound hReported

theorem emitRecords_correct {g : Graph}
    {frame : TraversalFrame g}
    {classes : Algorithms.Flubble.CycleClassAssignment g}
    {scan : Algorithms.Hairpin.HairpinScanAssignment g}
    {hierarchy : Algorithms.FlubbleTree.Hierarchy}
    {calls : List VariantCall}
    (hCalls : ReferenceCallSet frame scan hierarchy calls)
    (hHierarchy : Algorithms.FlubbleTree.IsCorrectHierarchy frame classes hierarchy)
    (hHairpinSound :
      ∀ boundary,
        boundary ∈ Algorithms.Hairpin.detectHairpins frame scan →
          Algorithms.Hairpin.IsPaperHairpinBoundary frame boundary) :
    EmissionCorrect frame classes scan hierarchy calls (emitRecords calls) := by
  refine
    { emitted := rfl
      wellFormed := ?_
      derived := ?_
      ordered := emitRecords_ordered hCalls.ordered }
  · intro record hRecord
    rcases List.mem_map.mp hRecord with ⟨call, hCallMem, hRecordEq⟩
    rw [← hRecordEq]
    exact recordOfCall_wellFormed (hCalls.wellFormed call hCallMem)
  · intro record hRecord
    rcases List.mem_map.mp hRecord with ⟨call, hCallMem, hRecordEq⟩
    rw [← hRecordEq]
    exact
      reportedSource_derived_of_correct
        hHierarchy
        hHairpinSound
        (hCalls.reported call hCallMem)

/--
Main GFA-facing theorem: accepted semantic GFA input plus verified upstream
flubble-tree and hairpin interfaces yield semantically correct VCF records.
-/
theorem emitVcfRecords_correct_for_gfa {doc : GFA.Document}
    {accepted : doc.Accepted}
    {frame : TraversalFrame doc.toGraph}
    {classes : Algorithms.Flubble.CycleClassAssignment doc.toGraph}
    {scan : Algorithms.Hairpin.HairpinScanAssignment doc.toGraph}
    {calls : List VariantCall}
    (hInput : Algorithms.Flubble.SupportedGFAInput doc accepted frame)
    (hClasses : Algorithms.Flubble.CycleClassAssignment.Correct frame classes)
    (hScan : Algorithms.Hairpin.HairpinScanAssignment.Correct frame scan)
    (hHierarchySupported :
      Algorithms.FlubbleTree.SupportedHierarchyInput
        (Algorithms.Flubble.candidateStack frame)
        (Algorithms.Flubble.detectFlubbles frame classes))
    (hCalls :
      ReferenceCallSet
        frame
        scan
        (Algorithms.FlubbleTree.buildHierarchy frame classes)
        calls) :
    EmissionCorrect
      frame
      classes
      scan
      (Algorithms.FlubbleTree.buildHierarchy frame classes)
      calls
      (emitRecords calls) := by
  have hHierarchy :
      Algorithms.FlubbleTree.IsCorrectHierarchy frame classes
        (Algorithms.FlubbleTree.buildHierarchy frame classes) :=
    Algorithms.FlubbleTree.buildHierarchy_correct
      (Algorithms.Flubble.SupportedGFAInput.toSupportedInput hInput)
      hClasses
      hHierarchySupported
  have hHairpin :=
    Algorithms.Hairpin.detectHairpins_correct_for_gfa
      (doc := doc)
      (accepted := accepted)
      (frame := frame)
      (scan := scan)
      hInput
      hScan
  exact emitRecords_correct hCalls hHierarchy hHairpin.1

end VCF
end PovuLean
