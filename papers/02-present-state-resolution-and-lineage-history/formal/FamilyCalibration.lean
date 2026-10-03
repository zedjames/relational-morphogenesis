import Mathlib.Data.Real.Basic
import Mathlib.Data.Finset.Card

/-! Logical controls for the finite family calibration. These establish algebra
and multiplicity relationships, not biological exchangeability or effect size. -/
namespace Tier30.Publication.RMMOPaper2.Controls

def threeFoldMinimum (r₁ r₂ r₃ : ℝ) : ℝ := min r₁ (min r₂ r₃)

theorem threeFoldMinimum_positive_iff (r₁ r₂ r₃ : ℝ) :
    0 < threeFoldMinimum r₁ r₂ r₃ ↔ 0 < r₁ ∧ 0 < r₂ ∧ 0 < r₃ := by
  simp only [threeFoldMinimum, lt_min_iff]

theorem material_threeFoldMinimum_iff (τ r₁ r₂ r₃ : ℝ) :
    τ < threeFoldMinimum r₁ r₂ r₃ ↔ τ < r₁ ∧ τ < r₂ ∧ τ < r₃ := by
  simp only [threeFoldMinimum, lt_min_iff]

theorem family_exceedance_card_dominates_member
    {ι : Type*} [DecidableEq ι] (bank : Finset ι)
    (member family : ι → ℝ) (threshold : ℝ)
    (dominates : ∀ b ∈ bank, member b ≤ family b) :
    (bank.filter (fun b => threshold ≤ member b)).card ≤
      (bank.filter (fun b => threshold ≤ family b)).card := by
  classical
  apply Finset.card_le_card
  intro b hb
  obtain ⟨hbank, hmember⟩ := Finset.mem_filter.mp hb
  exact Finset.mem_filter.mpr ⟨hbank, le_trans hmember (dominates b hbank)⟩

end Tier30.Publication.RMMOPaper2.Controls
