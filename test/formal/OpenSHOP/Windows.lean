import Mathlib

namespace OpenSHOP

/-- A single level band describes exactly the pairwise range restrictions.
The band endpoints need not be forced to equal the attained extrema. -/
theorem level_window_band {ι : Type*} (s : Finset ι) (hne : s.Nonempty)
    (z : ι → ℝ) (limit : ℝ) :
    (∃ floor : ℝ, ∀ i ∈ s, floor ≤ z i ∧ z i ≤ floor+limit) ↔
      ∀ i ∈ s, ∀ j ∈ s, z i-z j ≤ limit := by
  constructor
  · rintro ⟨floor, h⟩ i hi j hj
    obtain ⟨_, hu⟩ := h i hi
    obtain ⟨hl, _⟩ := h j hj
    linarith
  · intro h
    obtain ⟨j, hj, hmin⟩ := s.exists_min_image z hne
    refine ⟨z j, ?_⟩
    intro i hi
    exact ⟨hmin i hi, by linarith [h i hi j hj]⟩

/-- An affine deviation on a knot interval is bounded by its endpoint bounds.
For a constant flow and an affine rolling integral this justifies endpoint
tests. Knots must include jumps of the flow and shifted jumps of its history. -/
theorem rolling_affine_interval (a b weight up down : ℝ)
    (hw : 0 ≤ weight) (hw1 : weight ≤ 1)
    (ha : -down ≤ a ∧ a ≤ up) (hb : -down ≤ b ∧ b ≤ up) :
    -down ≤ (1-weight)*a+weight*b ∧ (1-weight)*a+weight*b ≤ up := by
  constructor
  · have h₁ := mul_nonneg (sub_nonneg.mpr hw1) (by linarith [ha.1] : 0 ≤ a+down)
    have h₂ := mul_nonneg hw (by linarith [hb.1] : 0 ≤ b+down)
    nlinarith
  · have h₁ := mul_nonneg (sub_nonneg.mpr hw1) (by linarith [ha.2] : 0 ≤ up-a)
    have h₂ := mul_nonneg hw (by linarith [hb.2] : 0 ≤ up-b)
    nlinarith

/-- Splitting a measured water allowance at a restart counts the prefix once. -/
theorem window_restart_additivity (before atRestart after : ℝ) :
    (atRestart-before)+(after-atRestart)=after-before := by ring

/-- Removing a common cumulative-volume datum leaves rolling means unchanged. -/
theorem rolling_datum_invariance (before after datum duration scale : ℝ) :
    ((after-datum)-(before-datum))/(scale*duration)=
      (after-before)/(scale*duration) := by ring

end OpenSHOP
