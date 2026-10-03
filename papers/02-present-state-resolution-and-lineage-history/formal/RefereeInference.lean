import Mathlib.Data.Real.Basic
import Mathlib.Tactic.NormNum

/-! Referee inference controls are countermodels and input-ownership results,
not axioms about MELA, cell identity, history irrelevance, or development. -/
namespace Tier30.Publication.RMMOPaper2.Referee

/-- A large survival collapse can be exactly the control's ordinary granularity. -/
theorem collapse_not_preferential :
    (3 / 100 : ℝ) < 1 / 10 ∧ ¬ ((3 / 100 : ℝ) < 3 / 100) := by
  constructor <;> norm_num

noncomputable def sharedNormalizedPredictor (stateCount targetCount : ℝ) : ℝ :=
  stateCount / (stateCount + targetCount)

/-- Disjoint feature identifiers do not remove target dependence from a shared
normalizer: changing only the target count changes the alleged state feature. -/
theorem disjoint_not_target_independent :
    (0 : Nat) ≠ 1 ∧ sharedNormalizedPredictor 1 1 ≠
      sharedNormalizedPredictor 1 3 := by
  constructor
  · decide
  · norm_num [sharedNormalizedPredictor]

def restrictedPrediction (_history : Bool) : Bool := false
def jointPrediction (history : Bool) : Bool := history

/-- No change for an unexposed restricted predictor is compatible with perfect
history prediction by another estimator on the same two-point carrier. -/
theorem restricted_null_not_general_null :
    restrictedPrediction false = restrictedPrediction true ∧
    jointPrediction false = false ∧ jointPrediction true = true ∧
    jointPrediction false ≠ jointPrediction true := by decide

structure InductiveFold (Cell : Type*) where
  training : Cell → Prop
  test : Cell → Prop
  disjoint : ∀ c, training c → ¬ test c
  discoveryInputs : Cell → Prop
  fitInputs : Cell → Prop
  tuningInputs : Cell → Prop
  discoveryOwned : ∀ c, discoveryInputs c → training c
  fitOwned : ∀ c, fitInputs c → training c
  tuningOwned : ∀ c, tuningInputs c → training c

/-- Ownership at all three stages excludes test observations at each stage.
The Python validator, not this abstract theorem, checks the recorded fit sets. -/
theorem inductive_excludes_test {Cell : Type*} (f : InductiveFold Cell) :
    ∀ c, (f.discoveryInputs c ∨ f.fitInputs c ∨ f.tuningInputs c) → ¬ f.test c := by
  intro c hc
  rcases hc with h | h | h
  · exact f.disjoint c (f.discoveryOwned c h)
  · exact f.disjoint c (f.fitOwned c h)
  · exact f.disjoint c (f.tuningOwned c h)

#print axioms collapse_not_preferential
#print axioms disjoint_not_target_independent
#print axioms restricted_null_not_general_null
#print axioms inductive_excludes_test
end Tier30.Publication.RMMOPaper2.Referee
