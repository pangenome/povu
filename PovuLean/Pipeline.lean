import PovuLean.Core
import PovuLean.GFA
import PovuLean.Algorithms
import PovuLean.VCF.Correctness

/-!
Integrated trusted proof path for the semantic povu reference pipeline.

This module is the synthesis-owned join point for the completed Lean proof
families.  It does not introduce a Rust/Lean conformance harness or byte-level
parsers.  Instead, it exposes the final semantic theorem boundary: accepted GFA
records, certified traversal/class/scan witnesses, supported flubble hierarchy
inputs, and well-formed pipeline-reported variant calls emit semantically
correct VCF records.
-/

namespace PovuLean
namespace Pipeline

open Core

/--
Integrated semantic GFA-to-VCF correctness theorem.

The remaining external/conformance work is precisely in supplying the
certificates consumed here from the current povu implementation and comparing
the implementation's serialized VCF output against these semantic records.
-/
theorem semanticGfaToVcf_correct {doc : GFA.Document}
    {accepted : doc.Accepted}
    {frame : TraversalFrame doc.toGraph}
    {classes : Algorithms.Flubble.CycleClassAssignment doc.toGraph}
    {scan : Algorithms.Hairpin.HairpinScanAssignment doc.toGraph}
    {calls : List VCF.VariantCall}
    (hInput : Algorithms.Flubble.SupportedGFAInput doc accepted frame)
    (hClasses : Algorithms.Flubble.CycleClassAssignment.Correct frame classes)
    (hScan : Algorithms.Hairpin.HairpinScanAssignment.Correct frame scan)
    (hHierarchySupported :
      Algorithms.FlubbleTree.SupportedHierarchyInput
        (Algorithms.Flubble.candidateStack frame)
        (Algorithms.Flubble.detectFlubbles frame classes))
    (hCalls :
      VCF.ReferenceCallSet
        frame
        scan
        (Algorithms.FlubbleTree.buildHierarchy frame classes)
        calls) :
    VCF.EmissionCorrect
      frame
      classes
      scan
      (Algorithms.FlubbleTree.buildHierarchy frame classes)
      calls
      (VCF.emitRecords calls) :=
  VCF.emitVcfRecords_correct_for_gfa
    hInput
    hClasses
    hScan
    hHierarchySupported
    hCalls

end Pipeline
end PovuLean
