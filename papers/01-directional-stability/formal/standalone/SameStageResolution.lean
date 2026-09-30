import Mathlib

namespace Tier30.Publication.RMMOPaper1

def sameStageMagnitudeReference : ℚ := 6731361414157754 / 100000000000000000

def chronologicalAdjacentScores : List ℚ :=
  [-58526 / 1000000, 39793 / 1000000, -2965 / 1000000,
   -32655 / 1000000, -22993 / 1000000, -47473 / 1000000,
   -14959 / 1000000]

def resolvedAboveSameStageReference (score : ℚ) : Bool :=
  decide (|score| > sameStageMagnitudeReference)

def chronologyAlignedSign (score : ℚ) : Bool :=
  decide (score > 0)

theorem same_stage_resolution_rule_applies_uniformly :
    chronologicalAdjacentScores.map resolvedAboveSameStageReference =
      chronologicalAdjacentScores.map fun score =>
        decide (|score| > sameStageMagnitudeReference) := by
  rfl

theorem all_adjacent_below_same_stage_reference :
    ∀ score ∈ chronologicalAdjacentScores,
      |score| ≤ sameStageMagnitudeReference := by
  norm_num [chronologicalAdjacentScores, sameStageMagnitudeReference]

theorem chronology_aligned_sign_count_eq_one :
    (chronologicalAdjacentScores.filter fun score => chronologyAlignedSign score).length = 1 := by
  native_decide

theorem zero_resolved_adjacent_intervals :
    (chronologicalAdjacentScores.filter fun score =>
      resolvedAboveSameStageReference score).length = 0 := by
  native_decide

end Tier30.Publication.RMMOPaper1
