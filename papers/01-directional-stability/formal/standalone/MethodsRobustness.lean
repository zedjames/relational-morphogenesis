import Mathlib

namespace Tier30.Publication.RMMOPaper1

/-- The registered estimator and a neighboring-model control have distinct authority. -/
inductive EstimatorAuthority where
  | registeredPrimary
  | neighboringModelSensitivity
  deriving DecidableEq, Repr

theorem neighboring_model_control_does_not_redefine_registered_primary :
    EstimatorAuthority.neighboringModelSensitivity ≠ EstimatorAuthority.registeredPrimary := by
  decide

/-- A finite methods audit separates the aggregate, local, and component conclusions. -/
structure MethodsRobustnessOutcome where
  negativeAggregateStable : Bool
  zeroLocalResolutionStable : Bool
  everyComponentSignStable : Bool
  deriving DecidableEq, Repr

def completedMethodsRobustness : MethodsRobustnessOutcome := ⟨true, true, false⟩

theorem aggregate_and_local_stability_do_not_force_component_sign_stability :
    completedMethodsRobustness.negativeAggregateStable = true ∧
      completedMethodsRobustness.zeroLocalResolutionStable = true ∧
      completedMethodsRobustness.everyComponentSignStable = false := by
  decide

/-- Section recurrence and embryo replication are different support authorities. -/
inductive SupportAuthority where
  | sectionRecurrent
  | embryoReplicated
  deriving DecidableEq, Repr

theorem section_recurrent_support_is_not_embryo_replicated_support :
    SupportAuthority.sectionRecurrent ≠ SupportAuthority.embryoReplicated := by
  decide

/-- The historical 32-field bank and the reported 1,024-field bank remain distinct. -/
inductive PotentialNullAuthority where
  | historicalThirtyTwo
  | independentOneThousandTwentyFour
  deriving DecidableEq, Repr

theorem reported_potential_null_is_not_the_historical_thirty_two_rewire_bank :
    PotentialNullAuthority.independentOneThousandTwentyFour ≠
      PotentialNullAuthority.historicalThirtyTwo := by
  decide

/-- Stable effect magnitude need not make a thresholded multiplicity count stable. -/
structure CompositionSensitivity where
  residualMagnitudeSmall : Bool
  holmCountInvariant : Bool
  deriving DecidableEq, Repr

def effortTwentyComposition : CompositionSensitivity := ⟨true, false⟩

theorem small_composition_residual_does_not_imply_holm_count_invariance :
    effortTwentyComposition.residualMagnitudeSmall = true ∧
      effortTwentyComposition.holmCountInvariant = false := by
  decide

end Tier30.Publication.RMMOPaper1
