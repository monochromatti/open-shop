import Mathlib

namespace OpenSHOP

/-- Each release retains its own mixture. Convex interpolation of reference
transfer fractions remains a fraction. -/
theorem transfer_fraction_bounds (a b w : ℝ)
    (ha : 0 ≤ a ∧ a ≤ 1) (hb : 0 ≤ b ∧ b ≤ 1) (hw : 0 ≤ w ∧ w ≤ 1) :
    0 ≤ (1-w)*a+w*b ∧ (1-w)*a+w*b ≤ 1 := by
  constructor
  · exact add_nonneg (mul_nonneg (sub_nonneg.mpr hw.2) ha.1) (mul_nonneg hw.1 hb.1)
  · nlinarith [mul_nonneg (sub_nonneg.mpr hw.2) (sub_nonneg.mpr ha.2),
      mul_nonneg hw.1 (sub_nonneg.mpr hb.2)]

/-- Flow times an affine reference blend is quadratic, including its exact
first and second derivatives' coefficients. -/
theorem released_transfer_quadratic (q reference a slope : ℝ) :
    q*(a+slope*(q-reference)) = slope*q^2+(a-slope*reference)*q := by ring

/-- Arrived and remaining volume partition an already released cohort. -/
theorem cohort_inventory (volume fraction : ℝ)
    (hv : 0 ≤ volume) (hf : 0 ≤ fraction ∧ fraction ≤ 1) :
    0 ≤ volume*fraction ∧ 0 ≤ volume*(1-fraction) ∧
      volume*fraction+volume*(1-fraction)=volume := by
  exact ⟨mul_nonneg hv hf.1, mul_nonneg hv (sub_nonneg.mpr hf.2), by ring⟩

/-- Adding adjacent observation windows preserves their integrated volume. -/
theorem adjacent_window_additivity (primitive : ℝ → ℝ) (a b c : ℝ) :
    (primitive b-primitive a)+(primitive c-primitive b)=primitive c-primitive a := by ring

/-- Splitting a constant release does not change its volume. This does not
assert invariance of a confluence's interval-averaged timing. -/
theorem constant_cohort_refinement (conversion rate a b c : ℝ) :
    conversion*(b-a)*rate+conversion*(c-b)*rate=conversion*(c-a)*rate := by ring

/-- Mixture weights sum exactly to one in the mathematical model. -/
theorem blend_unit_mass (w : ℝ) : (1-w)+w=1 := by ring


/-- A cumulative response gives nonnegative instantaneous contributions from
rectangular release cohorts. -/
theorem cohort_point_nonnegative (cdf : ℝ → ℝ) (monotone : Monotone cdf)
    (a b time q : ℝ) (ordered : a ≤ b) (hq : 0 ≤ q) :
    0 ≤ q*(cdf (time-a)-cdf (time-b)) := by
  exact mul_nonneg hq (sub_nonneg.mpr (monotone (by linarith)))

/-- No water arrives from a future cohort if the delay distribution is causal. -/
theorem future_cohort_zero (cdf : ℝ → ℝ) (causal : ∀ x ≤ 0, cdf x=0)
    (a b time q : ℝ) (ordered : a ≤ b) (future : time ≤ a) :
    q*(cdf (time-a)-cdf (time-b))=0 := by
  rw [causal (time-a) (by linarith), causal (time-b) (by linarith)]
  ring

/-- The cumulative primitive telescopes across a finite window partition.
Endpoint differences give integrated arrival volumes in real arithmetic. -/
theorem window_partition (primitive : ℕ → ℝ) (n : ℕ) :
    ∑ k ∈ Finset.range n, (primitive (k+1)-primitive k) = primitive n-primitive 0 := by
  induction n with
  | zero => simp
  | succ n ih => rw [Finset.sum_range_succ,ih]; ring

/-- Conservation on transport cells recovers the original control-window
balance. No new control decision is introduced by subdividing transport. -/
theorem transport_cell_aggregation {M : Type*} [Fintype M]
    (conversion duration localFlow share meanFlow : ℝ)
    (cellDuration cellFlow cellArrival : M → ℝ)
    (partition : duration = ∑ m, cellDuration m)
    (mean : duration*meanFlow = ∑ m, cellDuration m*cellFlow m)
    (balance : ∀ m, conversion*cellDuration m*cellFlow m =
      conversion*cellDuration m*localFlow+share*cellArrival m) :
    conversion*duration*meanFlow =
      conversion*duration*localFlow+share*∑ m, cellArrival m := by
  have total : conversion*(∑ m, cellDuration m*cellFlow m) =
      conversion*(∑ m, cellDuration m)*localFlow+share*∑ m, cellArrival m := by
    calc
      conversion*(∑ m, cellDuration m*cellFlow m) =
          ∑ m, conversion*cellDuration m*cellFlow m := by
        rw [Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro m _
        ring
      _ = ∑ m, (conversion*cellDuration m*localFlow+share*cellArrival m) := by
        apply Finset.sum_congr rfl
        intro m _
        exact balance m
      _ = conversion*(∑ m, cellDuration m)*localFlow+share*∑ m, cellArrival m := by
        simp only [Finset.sum_add_distrib, ← Finset.sum_mul, ← Finset.mul_sum]
  rw [← mean, ← partition] at total
  simpa only [mul_assoc] using total

end OpenSHOP
