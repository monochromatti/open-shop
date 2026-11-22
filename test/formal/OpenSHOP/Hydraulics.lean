import Mathlib

/-! Real-arithmetic identities for the declared discrete hydraulic model.
These theorems do not verify program implementations or floating-point evaluation. -/
namespace OpenSHOP

/-- Every signed internal edge credits exactly the water it withdraws. -/
theorem incidence_cancels {N E : Type*} [Fintype N] [Fintype E] [DecidableEq N]
    (source target : E → N) (q : E → ℝ) :
    ∑ n, ∑ e, ((if target e = n then q e else 0) -
      (if source e = n then q e else 0)) = 0 := by
  classical
  rw [Finset.sum_comm]
  apply Finset.sum_eq_zero
  intro e _
  rw [Finset.sum_sub_distrib]
  simp

/-- Reservoir and in-transit inventory together change only by external water.
Arrivals and releases use volumes, including the common 0.0036 Δt conversion. -/
theorem storage_and_transit_conservation {R : Type*} [Fintype R]
    (storageChange external : ℝ) (release arrival transitChange : R → ℝ)
    (storage : storageChange = external + ∑ r, arrival r - ∑ r, release r)
    (transit : ∀ r, transitChange r = release r - arrival r) :
    storageChange + ∑ r, transitChange r = external := by
  simp_rw [transit]
  rw [Finset.sum_sub_distrib, storage]
  ring

/-- Prescribed branch shares conserve upstream water; local injections remain
on their own branches and are not multiplied by the shares. -/
theorem junction_allocation {B : Type*} [Fintype B]
    (share localWater : B → ℝ) (upstream : ℝ) (shares : ∑ b, share b = 1) :
    ∑ b, (localWater b + share b * upstream) = (∑ b, localWater b) + upstream := by
  rw [Finset.sum_add_distrib, ← Finset.sum_mul, shares, one_mul]

/-- A positive opening multiplies head, equivalently dividing the resistance.
The closed-gate case is deliberately excluded. -/
theorem tunnel_opening_head_form (opening resistance dh q : ℝ)
    (ho : 0 < opening) :
    opening * dh - resistance * q * |q| = 0 ↔
      dh = (resistance / opening) * q * |q| := by
  have hn : opening ≠ 0 := ne_of_gt ho
  constructor
  · intro h
    calc
      dh = (resistance * q * |q|) / opening :=
        (eq_div_iff hn).mpr (by nlinarith [h])
      _ = (resistance / opening) * q * |q| := by ring
  · intro h
    have hdiv : dh = (resistance * q * |q|) / opening := by
      rw [h]
      ring
    have hmul := (eq_div_iff hn).mp hdiv
    nlinarith

/-- Dividing the scaled pressure residual by a positive head multiplier reports
head error in metres, without hiding errors at small openings. -/
theorem tunnel_pressure_error (opening resistance dh q : ℝ)
    (ho : 0 < opening) :
    abs (opening * dh - resistance * q * |q|) / opening =
      abs (dh - (resistance / opening) * q * |q|) := by
  calc
    abs (opening * dh - resistance * q * |q|) / opening =
        abs ((opening * dh - resistance * q * |q|) / opening) := by
      rw [abs_div, abs_of_pos ho]
    _ = abs (dh - (resistance / opening) * q * |q|) := by
      congr 1
      field_simp [ne_of_gt ho]

/-- A passive open tunnel dissipates hydraulic work for either flow direction.
A closed gate instead imposes q=0 and need not impose equality of heads. -/
theorem tunnel_dissipation (resistance q dh : ℝ) (hr : 0 ≤ resistance)
    (law : dh = resistance * q * |q|) :
    dh * q = resistance * |q| ^ 3 ∧ 0 ≤ dh * q := by
  have hsq : q ^ 2 = |q| ^ 2 := (sq_abs q).symm
  have heq : dh * q = resistance * |q| ^ 3 := by
    rw [law]
    calc
      resistance * q * |q| * q = resistance * |q| * q ^ 2 := by ring
      _ = resistance * |q| ^ 3 := by rw [hsq]; ring
  exact ⟨heq, heq ▸ mul_nonneg hr (pow_nonneg (abs_nonneg q) 3)⟩

/-- The work balance for a pump is the reciprocal-efficiency formula on its
physical domain. The balance itself needs no division at numerical trial points. -/
theorem pumping_work_balance (work power efficiency motor : ℝ)
    (he : 0 < efficiency) (hm : 0 < motor) :
    work = power * efficiency * motor ↔
      power = work / (efficiency * motor) := by
  have hn : efficiency * motor ≠ 0 := ne_of_gt (mul_pos he hm)
  constructor
  · intro h
    apply (eq_div_iff hn).mpr
    nlinarith [h]
  · intro h
    have hmul := (eq_div_iff hn).mp h
    nlinarith

/-- On a linear electrical-efficiency segment, the work equation has this
unique root when the denominator is positive. -/
theorem electrical_conversion_root (shaft intercept slope power : ℝ)
    (hd : 0 < 1 - shaft * slope) :
    power = shaft * (intercept + slope * power) ↔
      power = shaft * intercept / (1 - shaft * slope) := by
  have hn : 1 - shaft * slope ≠ 0 := ne_of_gt hd
  constructor
  · intro h
    apply (eq_div_iff hn).mpr
    nlinarith [h]
  · intro h
    have hmul := (eq_div_iff hn).mp h
    nlinarith

/-- A monotone delta-head curve through zero cannot send water uphill. -/
theorem channel_dissipation (flow : ℝ → ℝ) (monotone : Monotone flow)
    (atZero : flow 0 = 0) (dh : ℝ) : 0 ≤ dh * flow dh := by
  rcases le_total 0 dh with h | h
  · exact mul_nonneg h (by simpa [atZero] using monotone h)
  · exact mul_nonneg_of_nonpos_of_nonpos h (by simpa [atZero] using monotone h)

/-- Convex interpolation of two passive reference laws remains passive. -/
theorem reference_blend_passive (dh q₀ q₁ w : ℝ)
    (hw₀ : 0 ≤ w) (hw₁ : w ≤ 1)
    (h₀ : 0 ≤ dh * q₀) (h₁ : 0 ≤ dh * q₁) :
    0 ≤ dh * ((1-w)*q₀+w*q₁) := by
  nlinarith [mul_nonneg (sub_nonneg.mpr hw₁) h₀, mul_nonneg hw₀ h₁]

/-- At the same head and volume, turbine generation cannot exceed ideal work. -/
theorem generation_energy_bound (work eta electrical : ℝ)
    (hw : 0 ≤ work) (he₀ : 0 ≤ eta) (he₁ : eta ≤ 1)
    (hm₀ : 0 ≤ electrical) (hm₁ : electrical ≤ 1) :
    0 ≤ work * eta * electrical ∧ work * eta * electrical ≤ work := by
  constructor
  · positivity
  · calc
      work * eta * electrical ≤ work * eta := mul_le_of_le_one_right (mul_nonneg hw he₀) hm₁
      _ ≤ work := mul_le_of_le_one_right hw he₁

/-- The reciprocal pump law consumes at least ideal work, for positive
hydraulic and motor efficiency at most one. -/
theorem pumping_energy_bound (work eta motor : ℝ)
    (hw : 0 ≤ work) (he₀ : 0 < eta) (he₁ : eta ≤ 1)
    (hm₀ : 0 < motor) (hm₁ : motor ≤ 1) :
    work ≤ work / (eta * motor) := by
  have hp : 0 < eta * motor := mul_pos he₀ hm₀
  apply (le_div_iff₀ hp).2
  calc
    work * (eta * motor) ≤ work * 1 := mul_le_mul_of_nonneg_left
      (by nlinarith [mul_nonneg (le_of_lt he₀) (sub_nonneg.mpr hm₁)]) hw
    _ = work := mul_one work

/-- For a quadratic storage-level curve, midpoint work has a cubic remainder.
This quantifies a time-discretization error, rather than asserting exact energy
conservation for nonlinear storage-level curves. -/
theorem midpoint_potential_remainder (a b z slope curvature : ℝ) :
    (z*b+slope*b^2/2+curvature*b^3/3) -
      (z*a+slope*a^2/2+curvature*a^3/3) =
    (b-a)*(z+slope*((a+b)/2)+curvature*((a+b)/2)^2) +
      curvature*(b-a)^3/12 := by ring


/-- Unequal pumping/generating steps can create artificial cycle work under
midpoint quadrature for a curved storage-level relation, even with ideal units.
The path is 0 → 2 → 1 → 0; positive work means pumping into storage. -/
theorem midpoint_cycle_work (z slope curvature : ℝ) :
    2*(z+slope+curvature) -
      (z+slope*(3/2)+curvature*(3/2)^2) -
      (z+slope*(1/2)+curvature*(1/2)^2) = -curvature/2 := by ring

/-- Shifting all heads by one constant preserves passive tunnel equations.
This symmetry does not preserve absolute operating bounds or efficiency data. -/
theorem pressure_shift_invariance (source target shift : ℝ) :
    (source+shift)-(target+shift)=source-target := by ring

/-- A reference shift respects an absolute head bound exactly when it belongs
to that bound's translated interval. All such intervals must be intersected. -/
theorem pressure_shift_bounds (head shift lo hi : ℝ) :
    (lo ≤ head+shift ∧ head+shift ≤ hi) ↔ (lo-head ≤ shift ∧ shift ≤ hi-head) := by
  constructor <;> intro h <;> constructor <;> linarith [h.1,h.2]


/-- Positive conductance, linearly interpolated between three reference levels.
Flow here is measured as volume per interval, absorbing 0.0036 Δt. -/
noncomputable def referenceConductance (reference : ℝ) : ℝ :=
  if reference ≤ 9/4 then 3 else
  if reference ≤ 5/2 then 3-8*(reference-9/4) else
  if reference ≤ 11/4 then 1-(8/3)*(reference-5/2) else 1/3

theorem referenceConductance_positive (reference : ℝ) :
    0 < referenceConductance reference := by
  unfold referenceConductance
  split_ifs <;> linarith

theorem referenceConductance_passive (reference delta : ℝ) :
    0 ≤ delta*(referenceConductance reference*delta) := by
  have := referenceConductance_positive reference
  nlinarith [mul_nonneg (le_of_lt this) (sq_nonneg delta)]

/-- With initial volumes (3,1) and level Z(V)=V, the upper ending volume x
satisfies x=3-K((3+x)/2)*(x-1). The same passive table admits three distinct
ending volumes. This disproves uniqueness under passivity alone. -/
theorem passive_midpoint_multiple_roots :
    (3/2 : ℝ)=3-referenceConductance ((3+3/2)/2)*(3/2-1) ∧
    (2 : ℝ)=3-referenceConductance ((3+2)/2)*(2-1) ∧
    (5/2 : ℝ)=3-referenceConductance ((3+5/2)/2)*(5/2-1) ∧
    (3/2 : ℝ) ≠ 2 ∧ (2 : ℝ) ≠ 5/2 ∧ (3/2 : ℝ) ≠ 5/2 := by
  norm_num [referenceConductance]


/-- For two equal-slope reservoirs joined by constant conductance, the midpoint
head difference obeys (1+r)d_new=(1-r)d_old. A large positive step can reverse
the levels even though continuous passive equalization cannot cross equilibrium. -/
theorem midpoint_equalization_overshoot (r oldDifference newDifference : ℝ)
    (largeStep : 1 < r) (initialOrder : 0 < oldDifference)
    (step : (1+r)*newDifference=(1-r)*oldDifference) : newDifference < 0 := by
  have negative : (1-r)*oldDifference < 0 := mul_neg_of_neg_of_pos (by linarith) initialOrder
  nlinarith


/-- Summing the actual local balance equations, including zero-storage nodes,
and each reach's inventory equation leaves only external water. Rates have
already been converted to interval volumes. Boundary exchange is external. -/
theorem network_water_conservation {N E R : Type*}
    [Fintype N] [Fintype E] [Fintype R] [DecidableEq N]
    (edgeSource edgeTarget : E → N) (riverSource riverTarget : R → N)
    (nodeChange external : N → ℝ) (edgeVolume : E → ℝ)
    (release arrival transitChange : R → ℝ)
    (nodeBalance : ∀ n, nodeChange n = external n +
      (∑ e, ((if edgeTarget e=n then edgeVolume e else 0)-
        (if edgeSource e=n then edgeVolume e else 0))) +
      (∑ r, if riverTarget r=n then arrival r else 0) -
      (∑ r, if riverSource r=n then release r else 0))
    (transitBalance : ∀ r, transitChange r=release r-arrival r) :
    (∑ n, nodeChange n)+(∑ r, transitChange r)=∑ n, external n := by
  classical
  have arrivals : (∑ n, ∑ r, if riverTarget r=n then arrival r else 0)=∑ r, arrival r := by
    rw [Finset.sum_comm]
    apply Finset.sum_congr rfl
    intro r _
    simp
  have releases : (∑ n, ∑ r, if riverSource r=n then release r else 0)=∑ r, release r := by
    rw [Finset.sum_comm]
    apply Finset.sum_congr rfl
    intro r _
    simp
  have total : (∑ n, nodeChange n)=(∑ n, external n)+(∑ r, arrival r)-(∑ r, release r) := by
    simp_rw [nodeBalance]
    rw [Finset.sum_sub_distrib,Finset.sum_add_distrib,Finset.sum_add_distrib,
      incidence_cancels,arrivals,releases]
    ring
  exact storage_and_transit_conservation _ _ release arrival transitChange total transitBalance


/-- The positive/negative split in the global tunnel graph gives q|q| exactly
at binary direction states. It does not assert exactness of its relaxation. -/
theorem signed_square_split (direction positive negative capacity : ℝ)
    (binary : direction=0 ∨ direction=1)
    (hp : 0 ≤ positive) (hn : 0 ≤ negative)
    (posBound : positive ≤ capacity*direction)
    (negBound : negative ≤ capacity*(1-direction)) :
    (positive-negative)*|positive-negative|=positive^2-negative^2 := by
  rcases binary with h | h
  · subst direction
    have hz : positive=0 := by nlinarith
    rw [hz,abs_of_nonpos (by linarith)]
    ring
  · subst direction
    have hz : negative=0 := by nlinarith
    rw [hz,abs_of_nonneg (by linarith)]
    ring

/-- The binary positive-part rows used by dry/wet outlets force exactly the
nonnegative part of head, including the crest point. -/
theorem wet_branch_exact (direction head wet lo hi : ℝ)
    (binary : direction=0 ∨ direction=1) (nonnegative : 0 ≤ wet)
    (aboveHead : head ≤ wet) (upper : wet ≤ hi*direction)
    (branch : wet ≤ head-lo*(1-direction)) : wet=max head 0 := by
  rcases binary with h | h
  · subst direction
    have hz : wet=0 := by nlinarith
    have hh : head ≤ 0 := by linarith
    simp [hz,max_eq_right hh]
  · subst direction
    have he : wet=head := by nlinarith
    have hh : 0 ≤ head := by linarith
    simp [he,max_eq_left hh]

/-- Endpoint and mean-head power checks can all pass while interior power
exceeds capacity. Here eta(h)=(4-h)/3 and the held hydraulic coefficient is
109/1600 MW/m. The inequality is strict, so it is not a rounding example. -/
theorem endpoint_and_mean_power_miss_peak :
    let P := fun h : ℝ => (109/4800 : ℝ)*h*(4-h)
    P 1 < 9/100 ∧ P (7/2) < 9/100 ∧ P (9/4) < 9/100 ∧ 9/100 < P 2 := by
  norm_num

/-- A continuous head path reaches every head between its endpoints. Thus a
power bound along the path necessarily holds on that whole head span.
No converse is claimed: the path may turn beyond its endpoint heads. -/
theorem continuous_head_span_limit {h : ℝ → ℝ} {a b bound : ℝ}
    (continuous : ContinuousOn h (Set.uIcc a b)) (power : ℝ → ℝ)
    (limits : ∀ t ∈ Set.uIcc a b, power (h t) ≤ bound) :
    ∀ z ∈ Set.uIcc (h a) (h b), power z ≤ bound := by
  intro z hz
  obtain ⟨t, ht, rfl⟩ := intermediate_value_uIcc continuous hz
  exact limits t ht

/-- Multiplying a pump power bound by positive conversion efficiency is exact.
The coefficient includes both hydraulic and motor efficiency. -/
theorem pump_power_upper_iff (work efficiency capacity : ℝ)
    (positive : 0 < efficiency) :
    work / efficiency ≤ capacity ↔ work - capacity * efficiency ≤ 0 := by
  rw [div_le_iff₀ positive]
  constructor <;> intro h <;> linarith

/-- A strictly increasing electrical conversion map preserves power bounds.
For generation the map is P/e(P), and conversion equates it to shaft power. -/
theorem generator_power_upper_iff {conversion : ℝ → ℝ}
    (increasing : StrictMono conversion) (power shaft capacity : ℝ)
    (equation : conversion power = shaft) :
    power ≤ capacity ↔ shaft ≤ conversion capacity := by
  rw [← equation]
  exact increasing.le_iff_le.symm

/-- For increasing conversion, endpoint inequalities enforce the whole head
interval. The numerical monotonicity check is separate from this implication. -/
theorem monotone_head_limits (power : ℝ → ℝ) (a b lower upper : ℝ)
    (ordered : a ≤ b) (increasing : MonotoneOn power (Set.Icc a b))
    (at_left : lower ≤ power a) (at_right : power b ≤ upper) :
    ∀ h ∈ Set.Icc a b, lower ≤ power h ∧ power h ≤ upper := by
  intro h hh
  exact ⟨at_left.trans (increasing ⟨le_refl a, ordered⟩ hh hh.1),
    (increasing hh ⟨ordered, le_refl b⟩ hh.2).trans at_right⟩

/-- A reservoir outlet's release law already supplies its inflow lower bound.
The rating is nonnegative on the permitted head domain; controlled outlets use
capacity minus inflow as their rating. No differentiability is assumed. -/
theorem outlet_flow_bound_redundant (flow inflow opening rating : ℝ)
    (nonnegativeOpening : 0 ≤ opening) (nonnegativeRating : 0 ≤ rating) :
    (flow = inflow + opening * rating ∧ inflow ≤ flow) ↔
      flow = inflow + opening * rating := by
  constructor
  · exact And.left
  · intro law
    exact ⟨law, by nlinarith [mul_nonneg nonnegativeOpening nonnegativeRating]⟩

/-- At positive remaining capacity, a prescribed release interval is exactly
the corresponding gate interval. Capacity includes natural river inflow. -/
theorem controlled_outlet_gate_bounds (inflow capacity opening lower upper : ℝ)
    (remainingCapacity : inflow < capacity) :
    (lower ≤ inflow + (capacity-inflow)*opening ∧
      inflow + (capacity-inflow)*opening ≤ upper) ↔
    ((lower-inflow)/(capacity-inflow) ≤ opening ∧
      opening ≤ (upper-inflow)/(capacity-inflow)) := by
  have positive : 0 < capacity-inflow := sub_pos.mpr remainingCapacity
  rw [div_le_iff₀ positive, le_div_iff₀ positive]
  constructor <;> rintro ⟨h₁,h₂⟩ <;> constructor <;> nlinarith

/-- With no natural inflow and positive capacity, an off river state fixes
the controlled opening to zero. -/
theorem controlled_outlet_off (capacity opening : ℝ) (positive : 0 < capacity) :
    (0 ≤ capacity*opening ∧ capacity*opening ≤ 0) ↔ opening=0 := by
  constructor
  · rintro ⟨lower,upper⟩
    have zeroFlow : capacity*opening=0 := le_antisymm upper lower
    exact (mul_eq_zero.mp zeroFlow).resolve_left (ne_of_gt positive)
  · intro zeroOpening
    simp [zeroOpening]

end OpenSHOP
