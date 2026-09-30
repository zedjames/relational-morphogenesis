import Mathlib

namespace Tier30.Publication.RMMOPaper1

structure ChainNullDecision where
  bootstrapSignStable : Bool
  centeredBootstrapNullSeparated : Bool
  sameStageAccumulationNullSeparated : Bool
  exactSignNullSeparated : Bool
  deriving DecidableEq, Repr

def referee02ChainNullDecision : ChainNullDecision :=
  ⟨true, true, false, false⟩

theorem bootstrap_sign_stability_does_not_imply_chain_null_separation :
    referee02ChainNullDecision.bootstrapSignStable = true ∧
    referee02ChainNullDecision.sameStageAccumulationNullSeparated = false := by
  decide

inductive CalibrationStage where
  | early
  | late
  deriving DecidableEq, Repr

def stageSpecificReference : CalibrationStage → ℚ
  | .early => 1 / 50
  | .late => 3 / 20

theorem same_stage_reference_may_depend_on_stage :
    stageSpecificReference .early ≠ stageSpecificReference .late := by
  norm_num [stageSpecificReference]

def pooledReference : ℚ := 1 / 15

theorem pooled_reference_and_stage_specific_reference_are_distinct :
    pooledReference ≠ stageSpecificReference .early := by
  norm_num [pooledReference, stageSpecificReference]

def conservativeEndpointReference (left right : ℚ) : ℚ := max left right

theorem edge_specific_conservative_reference_is_endpoint_owned (left right : ℚ) :
    conservativeEndpointReference left right = left ∨
      conservativeEndpointReference left right = right := by
  exact max_choice left right

structure ConditionedScoreFamily where
  mappedSupport : ℚ
  sinkOnly : ℚ
  full : ℚ

def nonadditiveConditioningWitness : ConditionedScoreFamily := ⟨2, 1, 2⟩

theorem mapped_conditioning_is_not_additive_decomposition :
    nonadditiveConditioningWitness.full ≠
      nonadditiveConditioningWitness.mappedSupport +
        nonadditiveConditioningWitness.sinkOnly := by
  norm_num [nonadditiveConditioningWitness]

theorem conditioning_preserving_sign_does_not_establish_fractional_ownership :
    nonadditiveConditioningWitness.mappedSupport > 0 ∧
    nonadditiveConditioningWitness.full > 0 ∧
    nonadditiveConditioningWitness.mappedSupport /
      nonadditiveConditioningWitness.full = 1 ∧
    nonadditiveConditioningWitness.sinkOnly ≠ 0 := by
  norm_num [nonadditiveConditioningWitness]

inductive GraphNullMechanism where
  | endpointLabelPermutation
  | degreePreservingCurveball
  deriving DecidableEq, Repr

theorem endpoint_label_permutation_is_not_edge_switch_randomization :
    GraphNullMechanism.endpointLabelPermutation ≠
      GraphNullMechanism.degreePreservingCurveball := by
  decide

theorem common_weight_scaling_with_threshold_preserves_relation
    {distance epsilon scale : ℝ} (positive : 0 < scale) :
    distance ≤ epsilon ↔ scale * distance ≤ scale * epsilon := by
  constructor <;> intro h <;> nlinarith

theorem composition_effect_ratio_near_one_does_not_imply_large_excess_information :
    let observed : ℚ := 10001 / 10000
    let nullMean : ℚ := 1
    observed / nullMean > 1 ∧ observed - nullMean = 1 / 10000 := by
  norm_num

structure ToleranceAuthority where
  primary : ℚ
  sensitivity : List ℚ

def referee02ToleranceAuthority : ToleranceAuthority :=
  ⟨1 / 4, [3 / 20, 7 / 20]⟩

theorem tolerance_sensitivity_does_not_redefine_primary_authority :
    referee02ToleranceAuthority.primary = 1 / 4 ∧
      referee02ToleranceAuthority.sensitivity = [3 / 20, 7 / 20] := by
  norm_num [referee02ToleranceAuthority]

end Tier30.Publication.RMMOPaper1
