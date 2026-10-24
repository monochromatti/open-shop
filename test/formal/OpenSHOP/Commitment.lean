import Mathlib

namespace OpenSHOP

/-- A zero-or-minimum-flow disjunction has an unambiguous operating state. -/
theorem river_flow_disjunction (q minimum capacity state : ℝ)
    (positive : 0 < minimum) (binary : state=0 ∨ state=1)
    (lower : minimum*state ≤ q) (upper : q ≤ capacity*state) :
    (state=0 ∧ q=0) ∨ (state=1 ∧ minimum ≤ q ∧ 0 < q) := by
  rcases binary with h | h
  · left; constructor
    · exact h
    · rw [h] at lower upper; nlinarith
  · right; refine ⟨h, ?_, ?_⟩
    · simpa [h] using lower
    · rw [h] at lower; linarith

/-- The gross-head lower row imposes the stated limit precisely in the on
state; the off state retains the pre-existing physical lower bound. -/
theorem conditional_head_lower (head limit lower state : ℝ)
    (domain : lower ≤ head) (binary : state=0 ∨ state=1) :
    (limit-max 0 (limit-lower)*(1-state) ≤ head) ↔ (state=1 → limit ≤ head) := by
  rcases binary with h | h <;> subst state <;> simp only [sub_zero, sub_self, mul_one, mul_zero]
  · constructor
    · intro _ impossible; norm_num at impossible
    · intro _; have hm := le_max_right (0 : ℝ) (limit-lower); linarith
  · simp

/-- Restarted volume allowances count each executed interval exactly once. -/
theorem observed_volume_restart (initial consumed remaining : ℝ) :
    initial+(consumed+remaining)=(initial+consumed)+remaining := by ring

/-- The linear epigraph used for a rate-change charge encloses the exact
absolute change. Positive costs make its smallest feasible value exact. -/
theorem flow_change_epigraph (change variation : ℝ) :
    (change ≤ variation ∧ -change ≤ variation) ↔ |change| ≤ variation := by
  rw [abs_le]
  constructor <;> rintro ⟨h₁,h₂⟩ <;> constructor <;> linarith

/-- The four product rows are exact at either binary state, for signed x. -/
theorem binary_product_exact (u x z lo hi : ℝ)
    (binary : u = 0 ∨ u = 1) (bounds : lo ≤ x ∧ x ≤ hi) :
    (lo*u ≤ z ∧ z ≤ hi*u ∧ x-hi*(1-u) ≤ z ∧ z ≤ x-lo*(1-u)) ↔ z=u*x := by
  rcases binary with h | h <;> subst u <;> simp only [mul_zero, mul_one, zero_mul,
    one_mul, sub_zero, sub_self]
  · constructor
    · intro ⟨a,b,_,_⟩; exact le_antisymm b a
    · intro hz; subst z; constructor
      · exact le_refl 0
      constructor
      · exact le_refl 0
      constructor <;> linarith [bounds.1,bounds.2]
  · constructor
    · intro ⟨_,_,a,b⟩; exact le_antisymm b a
    · intro hz; subst z; exact ⟨bounds.1,bounds.2,le_refl _,le_refl _⟩

/-- Positive power implies positive head only with a nonnegative conversion
factor. Commitment alone supplies neither strict inequality. -/
theorem positive_power_requires_positive_head (p factor h : ℝ)
    (hp : 0 < p) (hf : 0 ≤ factor) (conversion : p = factor*h) : 0 < h := by
  by_contra hn
  have : factor*h ≤ 0 := mul_nonpos_of_nonneg_of_nonpos hf (le_of_not_gt hn)
  linarith

/-- A zero-power committed idle state need not have a nonnegative head. -/
theorem negative_head_idle_counterexample :
    ∃ u q p h eta : ℝ, u=1 ∧ q=0 ∧ p=0 ∧ h<0 ∧ eta=9/10 ∧
      p=(981/100000)*q*h*eta := by
  refine ⟨1,0,0,-5,9/10,?_,?_,?_,?_,?_,?_⟩ <;> norm_num

/-- The transition rows uniquely identify startup and shutdown with binary
adjacent states, even though transition variables need only be continuous. -/
theorem binary_transitions (previous current startup shutdown : ℝ)
    (hp : previous=0 ∨ previous=1) (hc : current=0 ∨ current=1)
    (hs : 0 ≤ startup) (hd : 0 ≤ shutdown)
    (balance : current-previous=startup-shutdown)
    (suCurrent : startup ≤ current) (suPrevious : startup ≤ 1-previous)
    (sdPrevious : shutdown ≤ previous) (sdCurrent : shutdown ≤ 1-current) :
    startup=max 0 (current-previous) ∧ shutdown=max 0 (previous-current) := by
  rcases hp with hp | hp <;> rcases hc with hc | hc <;>
    subst previous <;> subst current <;> norm_num at * <;> constructor <;> linarith

/-- A certificate requires a real enclosing upper bound on the full feasible
set. Neither a feasible schedule nor a local stationary point supplies it. -/
theorem absolute_gap_bound (incumbent upper candidate epsilon : ℝ)
    (enclosing : candidate ≤ upper) (gap : upper-incumbent ≤ epsilon) :
    candidate-incumbent ≤ epsilon := by linarith

/-- An equality is contained in every nonnegative band about its right side.
The band is used only for bounding; it does not redefine physical feasibility. -/
theorem equality_in_band (value rhs radius : ℝ) (nonnegative : 0 ≤ radius)
    (equation : value=rhs) : rhs-radius ≤ value ∧ value ≤ rhs+radius := by
  rw [equation]
  constructor <;> linarith

/-- An upper bound over banded equalities also bounds their exact solutions.
This proves the relaxation implication, not a floating-point solver's output. -/
theorem equality_band_upper_bound {X I : Type*}
    (objective : X → ℝ) (equation : I → X → ℝ) (rhs radius : I → ℝ) (upper : ℝ)
    (nonnegative : ∀ i, 0 ≤ radius i)
    (enclosing : ∀ x, (∀ i, rhs i-radius i ≤ equation i x ∧
      equation i x ≤ rhs i+radius i) → objective x ≤ upper) :
    ∀ x, (∀ i, equation i x=rhs i) → objective x ≤ upper := by
  intro x exactEquations
  apply enclosing x
  intro i
  exact equality_in_band _ _ _ (nonnegative i) (exactEquations i)


/-- Interval propagation of storage bounds follows directly from mass balance.
Backward propagation is the same equation solved for the preceding state. -/
theorem storage_interval_enclosure (a b delta alo ahi dlo dhi : ℝ)
    (state : alo ≤ a ∧ a ≤ ahi) (flow : dlo ≤ delta ∧ delta ≤ dhi)
    (balance : b=a+delta) :
    alo+dlo ≤ b ∧ b ≤ ahi+dhi ∧ b-dhi ≤ a ∧ a ≤ b-dlo := by
  rcases state with ⟨_,_⟩
  rcases flow with ⟨_,_⟩
  exact ⟨by linarith,by linarith,by linarith,by linarith⟩

/-- Signed on-head moments require bounds including both zero and the entire
running head interval, including its negative part. -/
theorem signed_onhead_bounds (u head lo hi : ℝ)
    (binary : u=0 ∨ u=1) (bounds : lo ≤ head ∧ head ≤ hi) :
    min 0 lo ≤ u*head ∧ u*head ≤ max 0 hi := by
  rcases binary with h | h <;> subst u
  · simp only [zero_mul]; exact ⟨min_le_left _ _,le_max_left _ _⟩
  · simp only [one_mul]
    exact ⟨le_trans (min_le_right _ _) bounds.1,le_trans bounds.2 (le_max_right _ _)⟩

end OpenSHOP
