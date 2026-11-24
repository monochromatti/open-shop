import Mathlib

namespace OpenSHOP

/-- The leading second-order errors give the coarse-grid Richardson factor.
The expansion is a hypothesis, not a bound on a computed trajectory. -/
theorem second_order_coarse_error (coarse fine reference coefficient step : ℝ)
    (hc : coarse = reference + coefficient * step^2)
    (hf : fine = reference + coefficient * (step/2)^2) :
    coarse - reference = (4/3 : ℝ) * (coarse - fine) := by
  rw [hc, hf]
  ring

/-- A child error budget must include the correction used for its retained mesh. -/
theorem corrected_error_budget (raw tolerance factor : ℝ) (positive : 0 < factor)
    (child : raw ≤ tolerance / factor) : factor * raw ≤ tolerance := by
  have h := (le_div_iff₀ positive).mp child
  simpa only [mul_comm] using h

/-- A flow-unit stopping allowance is exactly a volume-unit balance allowance.
This is a real-arithmetic unit conversion, not a floating-point error bound. -/
theorem storage_balance_rounding_units (before after rate step tolerance rounding : ℝ)
    (hs : 0 < step) :
    |(after-before)/step-rate| ≤ tolerance+rounding/step ↔
      |after-before-step*rate| ≤ step*tolerance+rounding := by
  have identity : (after-before)/step-rate = (after-before-step*rate)/step := by
    field_simp
  rw [identity, abs_div, abs_of_pos hs, div_le_iff₀ hs,
    add_mul, div_mul_cancel₀ _ (ne_of_gt hs)]
  rw [mul_comm tolerance step]

/-- The polynomial discrete gradient has no division by the volume change. -/
theorem quadratic_storage_chain_rule (a b z slope curvature : ℝ) :
    (b-a)*(z+slope*((a+b)/2)+curvature*(((a+b)/2)^2+(b-a)^2/12)) =
      (z*b+slope*b^2/2+curvature*b^3/3)-
      (z*a+slope*a^2/2+curvature*a^3/3) := by ring

theorem mean_head_reversal (a b z slope curvature : ℝ) :
    z+slope*((a+b)/2)+curvature*(((a+b)/2)^2+(b-a)^2/12) =
    z+slope*((b+a)/2)+curvature*(((b+a)/2)^2+(a-b)^2/12) := by ring

theorem mean_head_diagonal (a z slope curvature : ℝ) :
    z+slope*((a+a)/2)+curvature*(((a+a)/2)^2+(a-a)^2/12) =
    z+slope*a+curvature*a^2 := by ring

/-- Integrating a linear table cell gives its exact trapezoidal work. -/
theorem linear_cell_chain_rule (a b z slope : ℝ) :
    (b-a)*(z+slope*(a+b)/2) =
      (z*b+slope*b^2/2)-(z*a+slope*a^2/2) := by ring

/-- Exact potential from adjacent SOS2 weights, rather than a linear chord. -/
theorem potential_sos2_correction (a width z slope t : ℝ) :
    (1-t)*(z*a+slope*a^2/2)+
      t*(z*(a+width)+slope*(a+width)^2/2)-slope*width^2*t*(1-t)/2 =
    z*(a+width*t)+slope*(a+width*t)^2/2 := by ring

/-- The integral weights used for the first partial derivative are positive. -/
theorem mean_slope_weight (a b : ℝ) (_ha : 0 ≤ a) (hab : a ≤ b) (hb : b ≤ 1) :
    0 ≤ (b-a)*(1-(a+b)/2) ∧
    (b-a)*(1-(a+b)/2) = (b-b^2/2)-(a-a^2/2) := by
  constructor
  · exact mul_nonneg (sub_nonneg.mpr hab) (by linarith)
  · ring

/-- Given the table slopes and their half-mass integration weights, the
bounded divided difference retains every true mean head. -/
theorem mean_slope_bounds {I : Type*} [Fintype I]
    (weight slope : I → ℝ) (lo hi : ℝ)
    (nonnegative : ∀ i, 0 ≤ weight i)
    (mass : ∑ i, weight i = 1/2)
    (bounds : ∀ i, lo ≤ slope i ∧ slope i ≤ hi) :
    lo/2 ≤ ∑ i, weight i*slope i ∧ ∑ i, weight i*slope i ≤ hi/2 := by
  constructor
  · calc
      lo/2 = ∑ i, weight i*lo := by rw [← Finset.sum_mul,mass]; ring
      _ ≤ ∑ i, weight i*slope i := Finset.sum_le_sum fun i _ =>
        mul_le_mul_of_nonneg_left (bounds i).1 (nonnegative i)
  · calc
      (∑ i, weight i*slope i) ≤ ∑ i, weight i*hi := Finset.sum_le_sum fun i _ =>
        mul_le_mul_of_nonneg_left (bounds i).2 (nonnegative i)
      _ = hi/2 := by rw [← Finset.sum_mul,mass]; ring

/-- The additional divided-difference row determines the head on the diagonal;
the work row determines it off the diagonal. No arbitrary idle-storage head. -/
theorem mean_head_graph_exact (a b level h d potentialA potentialB : ℝ)
    (headRow : h=level+d*(b-a))
    (workRow : (b-a)*h=potentialB-potentialA) :
    (b=a → h=level) ∧
      (b≠a → h=(potentialB-potentialA)/(b-a)) := by
  constructor
  · intro same
    rw [same] at headRow
    simpa using headRow
  · intro different
    apply (eq_div_iff (sub_ne_zero.mpr different)).2
    nlinarith [workRow]

/-- Subdivision cannot change ideal work when each step obeys the chain rule. -/
theorem storage_work_telescopes (potential : ℕ → ℝ) (n : ℕ) :
    ∑ k ∈ Finset.range n, (potential (k+1)-potential k) = potential n-potential 0 := by
  induction n with
  | zero => simp
  | succ n ih => rw [Finset.sum_range_succ,ih]; ring

/-- Multiply each node's volume balance by its own mean head. Internal
transfers then carry the source-minus-target work, for any finite graph. -/
theorem weighted_incidence {N E : Type*} [Fintype N] [Fintype E] [DecidableEq N]
    (source target : E → N) (volume : E → ℝ) (head : N → ℝ) :
    (∑ n, head n * ∑ e, ((if target e=n then volume e else 0)-
      (if source e=n then volume e else 0))) =
    -(∑ e, volume e*(head (source e)-head (target e))) := by
  classical
  simp_rw [Finset.mul_sum]
  rw [Finset.sum_comm,← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro e _
  simp_rw [mul_sub,Finset.sum_sub_distrib,mul_ite,mul_zero]
  simp
  ring

/-- Passive closed networks dissipate potential. This statement covers direct
internal transfers; delayed reaches require separate energy-port accounting. -/
theorem closed_network_energy {N E : Type*}
    [Fintype N] [Fintype E] [DecidableEq N]
    (source target : E → N) (volume : E → ℝ) (head change : N → ℝ)
    (balance : ∀ n, change n=head n * ∑ e,
      ((if target e=n then volume e else 0)-(if source e=n then volume e else 0)))
    (passive : ∀ e, 0 ≤ volume e*(head (source e)-head (target e))) :
    ∑ n, change n ≤ 0 := by
  simp_rw [balance]
  rw [weighted_incidence]
  exact neg_nonpos.mpr (Finset.sum_nonneg fun e _ => passive e)

/-- Chain-rule coefficients for composing a two-input constitutive function
with affine storage paths. The supplied Hessian is a hypothesis; this does not
assert differentiability at a measured curve's knots or verify floating point. -/
theorem constitutive_hessian_pullback
    (aa ab bb aj ak bj bk u v : ℝ) :
    aa * (aj*u+ak*v)^2 + 2*ab*(aj*u+ak*v)*(bj*u+bk*v) + bb*(bj*u+bk*v)^2 =
      (aa*aj^2+2*ab*aj*bj+bb*bj^2)*u^2 +
      2*(aa*aj*ak+ab*(aj*bk+bj*ak)+bb*bj*bk)*u*v +
      (aa*ak^2+2*ab*ak*bk+bb*bk^2)*v^2 := by ring

/-- Differentiating the first integral-mean partial gives the implemented
same-endpoint curvature, assuming its first-partial identity and nonzero width. -/
theorem integral_mean_curvature_before (width mean head g slope : ℝ)
    (hw : width ≠ 0) (hg : mean-head = g*width) :
    ((g-slope)*width+(mean-head))/(width*width) = (2*g-slope)/width := by
  rw [hg]
  field_simp [hw]
  ring

/-- The analogous curvature at the other endpoint. -/
theorem integral_mean_curvature_after (width mean head g slope : ℝ)
    (hw : width ≠ 0) (hg : head-mean = g*width) :
    ((slope-g)*width-(head-mean))/(width*width) = (slope-2*g)/width := by
  rw [hg]
  field_simp [hw]
  ring

/-- The mixed partial follows from the two first-partial identities. -/
theorem integral_mean_curvature_mixed (width mean head ga gb : ℝ)
    (hw : width ≠ 0) (hg : mean-head = ga*width) :
    (gb*width-(mean-head))/(width*width) = (gb-ga)/width := by
  rw [hg]
  field_simp [hw]

/-- Duration-weighted operating averages preserve integrated flow or energy.
The cells partition a positive-duration interval; arithmetic here is exact. -/
theorem duration_weighted_integral {K : Type*} [Fintype K]
    (duration : ℝ) (dt rate : K → ℝ) (nonzero : duration ≠ 0) :
    duration * (∑ k, rate k * (dt k / duration)) = ∑ k, dt k * rate k := by
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro k _
  field_simp [nonzero]

/-- Holding a control on smaller cells and taking its original-window mean
returns the same control. Refining a state grid adds no discharge decision. -/
theorem held_control_mean {K : Type*} [Fintype K]
    (duration control : ℝ) (dt : K → ℝ) (nonzero : duration ≠ 0)
    (partition : ∑ k, dt k = duration) :
    (∑ k, control * (dt k / duration)) = control := by
  apply mul_left_cancel₀ nonzero
  rw [duration_weighted_integral duration dt (fun _ => control) nonzero]
  rw [← Finset.sum_mul, partition]

end OpenSHOP
