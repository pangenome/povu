/-!
Small cost and asymptotic vocabulary for proof-side runtime statements.

This module deliberately contains only arithmetic predicates and lemmas.  It
does not connect any C++, Rust, or executable Lean implementation to a runtime
claim; later bridge modules must supply concrete cost theorems as explicit
parameters before using this vocabulary for an implementation statement.
-/

namespace PovuLean
namespace Complexity

/-- A concrete linear bound with named constants. -/
def LinearBoundWith (cost n c k : Nat) : Prop :=
  cost ≤ c * n + k

/--
Existential linearity over a chosen size parameter.  Use `LinearBoundWith` when
the constants themselves should remain visible in a theorem statement.
-/
def LinearBound (cost n : Nat) : Prop :=
  ∃ c k, LinearBoundWith cost n c k

namespace LinearBoundWith

theorem weaken {cost n c k c' k' : Nat}
    (h : LinearBoundWith cost n c k)
    (hc : c ≤ c') (hk : k ≤ k') :
    LinearBoundWith cost n c' k' := by
  unfold LinearBoundWith at h ⊢
  exact Nat.le_trans h (Nat.add_le_add (Nat.mul_le_mul_right n hc) hk)

theorem trans_cost {cost cost' n c k : Nat}
    (hCost : cost ≤ cost')
    (hBound : LinearBoundWith cost' n c k) :
    LinearBoundWith cost n c k :=
  Nat.le_trans hCost hBound

/--
Sequential composition of two stages over the same size parameter.

If stage `A` costs at most `ca * n + ka` and stage `B` costs at most
`cb * n + kb`, the combined stage is bounded by
`(ca + cb) * n + (ka + kb)`.
-/
theorem seq {costA costB n ca ka cb kb : Nat}
    (hA : LinearBoundWith costA n ca ka)
    (hB : LinearBoundWith costB n cb kb) :
    LinearBoundWith (costA + costB) n (ca + cb) (ka + kb) := by
  unfold LinearBoundWith at hA hB ⊢
  calc
    costA + costB ≤ (ca * n + ka) + (cb * n + kb) :=
      Nat.add_le_add hA hB
    _ = (ca * n + cb * n) + (ka + kb) := by
      omega
    _ = (ca + cb) * n + (ka + kb) := by
      rw [Nat.add_mul]

/--
Add a stage whose cost is bounded by a constant independent of `n`.
-/
theorem seq_const {cost stage n c k stageK : Nat}
    (hCost : LinearBoundWith cost n c k)
    (hStage : stage ≤ stageK) :
    LinearBoundWith (cost + stage) n c (k + stageK) := by
  unfold LinearBoundWith at hCost ⊢
  calc
    cost + stage ≤ (c * n + k) + stageK :=
      Nat.add_le_add hCost hStage
    _ = c * n + (k + stageK) := by
      omega

/--
Reuse an output-size bound for a later proof obligation whose measured cost is
no larger than that output size.
-/
theorem reuse_output_size {cost outputSize n c k : Nat}
    (hCost : cost ≤ outputSize)
    (hOutput : LinearBoundWith outputSize n c k) :
    LinearBoundWith cost n c k :=
  trans_cost hCost hOutput

end LinearBoundWith

namespace LinearBound

theorem with_constants {cost n c k : Nat}
    (h : LinearBoundWith cost n c k) :
    LinearBound cost n :=
  ⟨c, k, h⟩

theorem trans_cost {cost cost' n : Nat}
    (hCost : cost ≤ cost')
    (hBound : LinearBound cost' n) :
    LinearBound cost n := by
  rcases hBound with ⟨c, k, h⟩
  exact ⟨c, k, LinearBoundWith.trans_cost hCost h⟩

theorem seq {costA costB n : Nat}
    (hA : LinearBound costA n)
    (hB : LinearBound costB n) :
    LinearBound (costA + costB) n := by
  rcases hA with ⟨ca, ka, hA⟩
  rcases hB with ⟨cb, kb, hB⟩
  exact ⟨ca + cb, ka + kb, LinearBoundWith.seq hA hB⟩

theorem seq_const {cost stage n : Nat}
    (hCost : LinearBound cost n)
    (hStage : ∃ stageK, stage ≤ stageK) :
    LinearBound (cost + stage) n := by
  rcases hCost with ⟨c, k, hCost⟩
  rcases hStage with ⟨stageK, hStage⟩
  exact ⟨c, k + stageK, LinearBoundWith.seq_const hCost hStage⟩

theorem reuse_output_size {cost outputSize n : Nat}
    (hCost : cost ≤ outputSize)
    (hOutput : LinearBound outputSize n) :
    LinearBound cost n :=
  trans_cost hCost hOutput

end LinearBound

/-- Named costs for two sequential proof stages. -/
structure SequentialCosts where
  first : Nat
  second : Nat

namespace SequentialCosts

def total (costs : SequentialCosts) : Nat :=
  costs.first + costs.second

theorem linear_total_with {costs : SequentialCosts}
    {n cFirst kFirst cSecond kSecond : Nat}
    (hFirst : LinearBoundWith costs.first n cFirst kFirst)
    (hSecond : LinearBoundWith costs.second n cSecond kSecond) :
    LinearBoundWith costs.total n (cFirst + cSecond) (kFirst + kSecond) :=
  LinearBoundWith.seq hFirst hSecond

end SequentialCosts

end Complexity
end PovuLean
