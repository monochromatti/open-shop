import Mathlib

namespace OpenSHOP

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

end OpenSHOP
