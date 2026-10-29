import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/- A prescribed arrival rate 12 + 24t m³/s, where t is hours in [0,1].
   Its cumulative arrival volume is 43200(t+t²) m³.
   The reservoir starts at 2000 m³, with limits 400 and 3600 m³.
   These statements concern this explicit mass-balance trajectory, not the
   complete OpenSHOP model, a solver, or any benchmark's infeasibility. -/
namespace OpenSHOP.StoragePaths

def volume (q t : ℝ) : ℝ :=
  2000 + 43200 * (t + t ^ 2) - 3600 * q * t

theorem hourly_endpoints_pass :
    volume 24 0 = 2000 ∧ volume 24 1 = 2000 := by
  norm_num [volume]

theorem midpoint_is_dry : volume 24 (1 / 2) = -8800 := by
  norm_num [volume]

theorem constant_release_impossible (q : ℝ) :
    ¬ (400 ≤ volume q (1 / 2) ∧ volume q 1 ≤ 3600) := by
  intro h
  norm_num [volume] at h
  linarith [h.1, h.2]

def following_volume (t : ℝ) : ℝ :=
  2000 + 43200 * (t + t ^ 2) - 3600 * (12 * t + 12 * t ^ 2)

theorem following_release_preserves_storage (t : ℝ) :
    following_volume t = 2000 := by
  unfold following_volume
  ring

theorem following_release_within_flow_limits (t : ℝ)
    (h0 : 0 ≤ t) (h1 : t ≤ 1) :
    2 ≤ 12 + 24 * t ∧ 12 + 24 * t ≤ 64 := by
  constructor <;> linarith

end OpenSHOP.StoragePaths
