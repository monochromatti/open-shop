import Mathlib

namespace OpenSHOP

/-- The difference from an affine slope on one release-reference interval is
a quadratic. Its extremum on the interval occurs at an endpoint or at the
stationary point, clamped to the interval. -/
theorem quadratic_flux_enclosure (a b lo hi q stationary : ℝ)
    (hlo : lo ≤ q) (hhi : q ≤ hi) (vertex : a = 0 ∨ 2*a*stationary + b = 0) :
    let f := fun z : ℝ => a*z^2 + b*z
    let z := max lo (min hi stationary)
    min (min (f lo) (f hi)) (f z) ≤ f q ∧
      f q ≤ max (max (f lo) (f hi)) (f z) := by
  intro f z
  rcases vertex with zero | critical
  · have hl : min (f lo) (f hi) ≤ f q := by
      by_cases hb : 0 ≤ b
      · exact (min_le_left _ _).trans (by simpa [f, zero] using mul_le_mul_of_nonneg_left hlo hb)
      · apply (min_le_right _ _).trans
        simpa [f, zero] using mul_le_mul_of_nonpos_left hhi (le_of_lt (lt_of_not_ge hb))
    have hu : f q ≤ max (f lo) (f hi) := by
      by_cases hb : 0 ≤ b
      · exact le_trans (by simpa [f, zero] using mul_le_mul_of_nonneg_left hhi hb) (le_max_right _ _)
      · apply le_trans _ (le_max_left _ _)
        simpa [f, zero] using mul_le_mul_of_nonpos_left hlo (le_of_lt (lt_of_not_ge hb))
    exact ⟨(min_le_left _ _).trans hl, hu.trans (le_max_left _ _)⟩
  have width : lo ≤ hi := hlo.trans hhi
  have zl : lo ≤ z := le_max_left _ _
  have zh : z ≤ hi := max_le width (min_le_left _ _)
  have closest : (z-stationary)^2 ≤ (q-stationary)^2 := by
    by_cases left : stationary ≤ lo
    · have eq : z = lo := max_eq_left ((min_le_right _ _).trans left)
      rw [eq]
      nlinarith [mul_nonneg (sub_nonneg.mpr hlo) (show 0 ≤ q+lo-2*stationary by linarith)]
    · by_cases right : hi ≤ stationary
      · have eq : z = hi := by simp [z, min_eq_left right, max_eq_right width]
        rw [eq]
        nlinarith [mul_nonneg (sub_nonneg.mpr hhi) (show 0 ≤ 2*stationary-q-hi by linarith)]
      · have eq : z = stationary := by
          simp [z, min_eq_right (le_of_lt (lt_of_not_ge right)),
            max_eq_right (le_of_lt (lt_of_not_ge left))]
        rw [eq]
        simpa using sq_nonneg (q-stationary)
  have farthest :
      (q-stationary)^2 ≤ max ((lo-stationary)^2) ((hi-stationary)^2) := by
    by_cases side : q ≤ stationary
    · apply le_trans _ (le_max_left _ _)
      nlinarith [mul_nonneg (sub_nonneg.mpr hlo) (show 0 ≤ 2*stationary-q-lo by linarith)]
    · apply le_trans _ (le_max_right _ _)
      nlinarith [mul_nonneg (sub_nonneg.mpr hhi)
        (show 0 ≤ hi+q-2*stationary by linarith [le_of_lt (lt_of_not_ge side)])]
  have diff (x y : ℝ) : f x - f y =
      a * ((x-stationary)^2 - (y-stationary)^2) := by
    have hb : b = -(2*a*stationary) := by linarith
    dsimp [f]
    rw [hb]
    ring
  by_cases convex : 0 ≤ a
  · constructor
    · apply (min_le_right _ _).trans
      have sign := mul_nonneg convex (sub_nonneg.mpr closest)
      linarith [diff q z]
    · apply le_trans _ (le_max_left _ _)
      rcases max_choice ((lo-stationary)^2) ((hi-stationary)^2) with h | h
      · rw [h] at farthest
        apply le_trans _ (le_max_left _ _)
        linarith [diff lo q, mul_nonneg convex (sub_nonneg.mpr farthest)]
      · rw [h] at farthest
        apply le_trans _ (le_max_right _ _)
        linarith [diff hi q, mul_nonneg convex (sub_nonneg.mpr farthest)]
  · have concave : a ≤ 0 := le_of_lt (lt_of_not_ge convex)
    constructor
    · apply (min_le_left _ _).trans
      rcases max_choice ((lo-stationary)^2) ((hi-stationary)^2) with h | h
      · rw [h] at farthest
        apply (min_le_left _ _).trans
        linarith [diff lo q, mul_nonpos_of_nonpos_of_nonneg concave (sub_nonneg.mpr farthest)]
      · rw [h] at farthest
        apply (min_le_right _ _).trans
        linarith [diff hi q, mul_nonpos_of_nonpos_of_nonneg concave (sub_nonneg.mpr farthest)]
    · apply le_trans _ (le_max_right _ _)
      linarith [diff q z, mul_nonpos_of_nonpos_of_nonneg concave (sub_nonneg.mpr closest)]

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

/-- Arrived and remaining volume partition an already released release. -/
theorem release_inventory (volume fraction : ℝ)
    (hv : 0 ≤ volume) (hf : 0 ≤ fraction ∧ fraction ≤ 1) :
    0 ≤ volume*fraction ∧ 0 ≤ volume*(1-fraction) ∧
      volume*fraction+volume*(1-fraction)=volume := by
  exact ⟨mul_nonneg hv hf.1, mul_nonneg hv (sub_nonneg.mpr hf.2), by ring⟩

/-- Adding adjacent observation windows preserves their integrated volume. -/
theorem adjacent_window_additivity (primitive : ℝ → ℝ) (a b c : ℝ) :
    (primitive b-primitive a)+(primitive c-primitive b)=primitive c-primitive a := by ring

/-- Splitting a constant release does not change its volume. This does not
assert invariance of a river junction's interval-averaged timing. -/
theorem constant_release_refinement (conversion rate a b c : ℝ) :
    conversion*(b-a)*rate+conversion*(c-b)*rate=conversion*(c-a)*rate := by ring

/-- Mixture weights sum exactly to one in the mathematical model. -/
theorem blend_unit_mass (w : ℝ) : (1-w)+w=1 := by ring


/-- A cumulative response gives nonnegative instantaneous contributions from
constant-rate releases. -/
theorem release_arrival_nonnegative (cdf : ℝ → ℝ) (monotone : Monotone cdf)
    (a b time q : ℝ) (ordered : a ≤ b) (hq : 0 ≤ q) :
    0 ≤ q*(cdf (time-a)-cdf (time-b)) := by
  exact mul_nonneg hq (sub_nonneg.mpr (monotone (by linarith)))

/-- No water arrives from a future release if the delay distribution is causal. -/
theorem future_release_zero (cdf : ℝ → ℝ) (causal : ∀ x ≤ 0, cdf x=0)
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

/-- Anchored integral recurrences uniquely recover the explicit finite sum.
The increments may be nonlinear reference flows evaluated at any fixed point. -/
theorem anchored_reference_integral (primitive increment : ℕ → ℝ)
    (anchor : primitive 0 = 0)
    (step : ∀ k, primitive (k+1)-primitive k = increment k) (n : ℕ) :
    primitive n = ∑ k ∈ Finset.range n, increment k := by
  induction n with
  | zero => simpa using anchor
  | succ n ih =>
      rw [Finset.sum_range_succ, ← ih]
      linarith [step n]

/-- An affine observation weight on a constant release is represented exactly
by its volume and normalized first time moment. No quadrature approximation
is introduced by this algebraic change of variables. -/
theorem affine_reference_moments (conversion rate origin duration a b wa wb : ℝ)
    (hab : b-a ≠ 0) (hd : duration ≠ 0) :
    conversion*rate*(b-a)*(wa+wb)/2 =
      (wa+(origin-a)*(wb-wa)/(b-a))*(conversion*rate*(b-a)) +
      duration*(wb-wa)/(b-a)*
        (conversion*rate*((b-origin)^2-(a-origin)^2)/(2*duration)) := by
  field_simp
  ring

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

/-- Fixed external inputs determine at most one state of an acyclic river
network. Each transfer map may be nonlinear, but depends only on earlier
reaches in topological order. This does not assert existence or uniqueness of
the coupled hydraulic problem. -/
theorem prescribed_routing_unique {n : ℕ} {A : Type*}
    (F : (Fin n → A) → Fin n → A)
    (dependency : ∀ i x y, (∀ j, j.val < i.val → x j = y j) → F x i = F y i)
    (x y : Fin n → A) (hx : ∀ i, x i = F x i) (hy : ∀ i, y i = F y i) : x = y := by
  have same : ∀ k : ℕ, ∀ hk : k < n, x ⟨k,hk⟩ = y ⟨k,hk⟩ := by
    intro k
    induction k using Nat.strong_induction_on with
    | h k ih =>
      intro hk
      rw [hx, hy]
      apply dependency
      intro j hj
      exact ih j.val hj j.isLt
  funext i
  exact same i.val i.isLt

/-- Substituting a prescribed routing solution preserves every remaining
physical restriction. The hypothesis states that the inputs are fixed. -/
theorem prescribed_routing_substitution {X R : Type*} (inputs : X → R)
    (route : R → R) (P : X → R → Prop) (fixed : R) (x : X)
    (hfixed : inputs x = fixed) :
    (∃ r, r = route (inputs x) ∧ P x r) ↔ P x (route fixed) := by
  rw [hfixed]
  simp

/-- Nonnegative released water and transfer weights already imply a zero
arrival minimum, for either volume or instantaneous-rate observations.
The history term is water already in transit. -/
theorem arrival_zero_bound_redundant {I : Type*} [Fintype I]
    (arrival history : ℝ) (weight release : I → ℝ)
    (nonnegativeHistory : 0 ≤ history) (nonnegativeWeight : ∀ i, 0 ≤ weight i)
    (nonnegativeRelease : ∀ i, 0 ≤ release i) :
    (arrival = history + ∑ i, weight i * release i ∧ 0 ≤ arrival) ↔
      arrival = history + ∑ i, weight i * release i := by
  constructor
  · exact And.left
  · intro law
    refine ⟨law, ?_⟩
    rw [law]
    exact add_nonneg nonnegativeHistory
      (Finset.sum_nonneg fun i _ => mul_nonneg (nonnegativeWeight i) (nonnegativeRelease i))

end OpenSHOP
