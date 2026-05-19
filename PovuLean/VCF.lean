import PovuLean.VCF.Spec
import PovuLean.VCF.Emit
import PovuLean.VCF.Correctness
import PovuLean.VCF.Examples

/-!
Import aggregator for semantic VCF modules.

The VCF family defines the supported semantic record subset, the pure reference
emitter from verified variant calls, and correctness theorems connecting emitted
records to verified flubble-tree and hairpin sources.  Byte-level formatting and
Rust/Lean fixture conformance are deliberately outside this trusted import path.
-/

namespace PovuLean
namespace VCF

end VCF
end PovuLean
