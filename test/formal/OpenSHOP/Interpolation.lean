import Mathlib

namespace OpenSHOP

/-- Cubic Hermite coefficients used in the original discharge cell. -/
theorem hermite_coefficients (a b u v t : ℝ) :
    (2*t^3-3*t^2+1)*a+(t^3-2*t^2+t)*u+
      (-2*t^3+3*t^2)*b+(t^3-t^2)*v =
    a+u*t+(3*(b-a)-2*u-v)*t^2+(2*(a-b)+u+v)*t^3 := by ring

/-- The quadratic-product SOS2 formulation reproduces the cubic, not its chord.
Here t is the right weight on the single selected adjacent pair. -/
theorem cubic_sos2_identity (a₀ a₁ a₂ a₃ t : ℝ) :
    a₀*(1-t)+(a₀+a₁+a₂+a₃)*t+
      (-a₂-a₃)*(1-t)^2*t+(-a₂-2*a₃)*(1-t)*t^2 =
    a₀+a₁*t+a₂*t^2+a₃*t^3 := by ring

/-- The second cubic basis needs no independent variable on a selected pair. -/
theorem cubic_second_basis (t : ℝ) :
    (1-t)*t-(1-t)^2*t = (1-t)*t^2 := by ring

theorem adjacent_weight_product_bound (t : ℝ) (h₀ : 0 ≤ t) (h₁ : t ≤ 1) :
    0 ≤ (1-t)*t ∧ (1-t)*t ≤ 1/4 := by
  constructor
  · exact mul_nonneg (sub_nonneg.mpr h₁) h₀
  · nlinarith [sq_nonneg (t-1/2)]

theorem cubic_basis_bound (t : ℝ) (h₀ : 0 ≤ t) (h₁ : t ≤ 1) :
    0 ≤ (1-t)^2*t ∧ (1-t)^2*t ≤ 4/27 := by
  constructor
  · positivity
  · have hp : 0 ≤ (t-1/3)^2*(4/3-t) :=
      mul_nonneg (sq_nonneg _) (by linarith)
    nlinarith

/-- Affine reparameterization retains the polynomial on each common cell. -/
theorem shifted_cubic (a₀ a₁ a₂ a₃ origin scale t : ℝ) :
    a₀+a₁*(origin+scale*t)+a₂*(origin+scale*t)^2+a₃*(origin+scale*t)^3 =
    (a₀+origin*(a₁+origin*(a₂+origin*a₃)))+
    (scale*(a₁+origin*(2*a₂+3*origin*a₃)))*t+
    (scale^2*(a₂+3*origin*a₃))*t^2+(scale^3*a₃)*t^3 := by ring

/-- Product weights reproduce a trilinear product for arbitrary signed axes.
No nonnegative-head assumption is needed. -/
theorem tensor_product_moment
    (w x v y r z : Fin 2 → ℝ) :
    ∑ i, ∑ j, ∑ k, (w i*v j*r k)*(x i*y j*z k) =
      (∑ i, w i*x i)*(∑ j, v j*y j)*(∑ k, r k*z k) := by
  simp only [Fin.sum_univ_two]
  ring

/-- Nonnegative partition weights enclose a polynomial by its coefficient
bounds. Applying this to the Bernstein basis gives the support inequality;
it does not certify the floating-point coefficient calculation. -/
theorem partition_upper_bound {I : Type*} [Fintype I]
    (weight coefficient : I → ℝ) (upper : ℝ)
    (nonnegative : ∀ i, 0 ≤ weight i) (partition : ∑ i, weight i = 1)
    (enclosure : ∀ i, coefficient i ≤ upper) :
    ∑ i, weight i*coefficient i ≤ upper := by
  calc
    ∑ i, weight i*coefficient i ≤ ∑ i, weight i*upper :=
      Finset.sum_le_sum (fun i _ => mul_le_mul_of_nonneg_left (enclosure i) (nonnegative i))
    _ = upper := by rw [← Finset.sum_mul, partition, one_mul]


/-- Two endpoint weights reproduce the physical coordinate for signed bounds. -/
theorem endpoint_moment (lo hi t : ℝ) :
    (1-t)*lo+t*hi = lo+(hi-lo)*t := by ring

/-- The rank-one weights for the eight corners sum to one. -/
theorem corner_partition (a b c : ℝ) :
    (1-a)*(1-b)*(1-c)+(1-a)*(1-b)*c+
    (1-a)*b*(1-c)+(1-a)*b*c+
    a*(1-b)*(1-c)+a*(1-b)*c+a*b*(1-c)+a*b*c=1 := by ring

/-- Electrical-efficiency bounds may multiply a nonnegative shaft product.
For negative head this requires conversion and P≥0, not just commitment. -/
theorem electrical_product_bounds (shaft electrical lo hi : ℝ)
    (hs : 0 ≤ shaft) (he : lo ≤ electrical ∧ electrical ≤ hi) :
    lo*shaft ≤ electrical*shaft ∧ electrical*shaft ≤ hi*shaft :=
  ⟨mul_le_mul_of_nonneg_right he.1 hs, mul_le_mul_of_nonneg_right he.2 hs⟩


/-- Positive intercepts make P/η(P) strictly increasing on an electrical
linear-efficiency segment with positive endpoint efficiencies. -/
theorem electrical_shaft_monotone (a b p₁ p₂ : ℝ)
    (intercept : 0 < a) (ordered : p₁ < p₂)
    (eta₁ : 0 < a+b*p₁) (eta₂ : 0 < a+b*p₂) :
    p₁/(a+b*p₁) < p₂/(a+b*p₂) := by
  apply (div_lt_div_iff₀ eta₁ eta₂).2
  nlinarith [mul_pos intercept (sub_pos.mpr ordered)]

/-- The implicit electrical equation has the stated unique segment root
whenever its denominator is nonzero. Physical extension checks require it
positive, as well as positive electrical efficiency. -/
theorem electrical_segment_root (a b shaft p : ℝ) (denominator : 1-shaft*b ≠ 0) :
    p=shaft*(a+b*p) ↔ p=shaft*a/(1-shaft*b) := by
  constructor
  · intro h
    apply (eq_div_iff denominator).2
    nlinarith
  · intro h
    have := (eq_div_iff denominator).1 h
    nlinarith


/-- Power-basis to Bernstein conversion at the discharge degree used by the
power supports. These are exactly the coefficients for degree four. -/
theorem bernstein_four_identity (a₀ a₁ a₂ a₃ a₄ t : ℝ) :
    a₀*(1-t)^4 + 4*(a₀+a₁/4)*t*(1-t)^3 +
    6*(a₀+a₁/2+a₂/6)*t^2*(1-t)^2 +
    4*(a₀+3*a₁/4+a₂/2+a₃/4)*t^3*(1-t) +
    (a₀+a₁+a₂+a₃+a₄)*t^4 = a₀+a₁*t+a₂*t^2+a₃*t^3+a₄*t^4 := by ring

/-- The corresponding degree-two conversion for the head coordinate. -/
theorem bernstein_two_identity (a₀ a₁ a₂ t : ℝ) :
    a₀*(1-t)^2+2*(a₀+a₁/2)*t*(1-t)+(a₀+a₁+a₂)*t^2 =
    a₀+a₁*t+a₂*t^2 := by ring

/-- The discharge Bernstein basis is a nonnegative partition on its cell. -/
theorem bernstein_four_partition (t : ℝ) (h₀ : 0 ≤ t) (h₁ : t ≤ 1) :
    0 ≤ (1-t)^4 ∧ 0 ≤ 4*t*(1-t)^3 ∧ 0 ≤ 6*t^2*(1-t)^2 ∧
    0 ≤ 4*t^3*(1-t) ∧ 0 ≤ t^4 ∧
    (1-t)^4+4*t*(1-t)^3+6*t^2*(1-t)^2+4*t^3*(1-t)+t^4=1 := by
  have : 0 ≤ 1-t := sub_nonneg.mpr h₁
  refine ⟨by positivity,by positivity,by positivity,by positivity,by positivity,?_⟩
  ring

/-- The head Bernstein basis is also a nonnegative partition. -/
theorem bernstein_two_partition (t : ℝ) (h₀ : 0 ≤ t) (h₁ : t ≤ 1) :
    0 ≤ (1-t)^2 ∧ 0 ≤ 2*t*(1-t) ∧ 0 ≤ t^2 ∧
    (1-t)^2+2*t*(1-t)+t^2=1 := by
  have : 0 ≤ 1-t := sub_nonneg.mpr h₁
  exact ⟨by positivity,by positivity,by positivity,by ring⟩

end OpenSHOP
