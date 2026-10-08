# Proof-speed experiment protocol

Baseline: OpenSHOP 0.4.3, commit 7743702366ad96ce2ecfa24881c66a44a887c9a5.
Dependencies remain fixed by the same Manifest SHA256.

Screen 1, commit 57cf943a6f1ee0f6f251166faede10968208ebeb,
https://github.com/monochromatti/open-shop/actions/runs/37729035609:
Seven policies, six cases, 120 seconds, one repetition. Separate native OBBT
projection/filter ablations, shared head lift, affine discharge lift and both.

Screen 2, commit 3eb507390e6f14408715995d072679d222c13a0a,
https://github.com/monochromatti/open-shop/actions/runs/37730112283:
Five policies, six cases, 120 seconds, one repetition. Baseline, three native
parameter combinations with the shared lift, and bounded certified root cuts.

Refinement, commit c4397c6909516c9ebd26f368323296adcf357f4d,
https://github.com/monochromatti/open-shop/actions/runs/37732208118:
Five policies, six cases, 120 seconds, one repetition. Baseline; limited cuts
with/without bilinear projection LPs; wider cuts (128 cuts, eight rounds, 256
coordinate checks) with/without projection LPs. Both keep the same 1.5-second
polynomial-certificate allowance. This follows the second screen's accepted
cuts improving all four operating bounds while using 0.12–0.18 seconds and
examining only 64 of 672 available coordinates on the 24-hour cases.

Within each matrix job: identical frozen case, initial controls, thread counts,
model equations, table supports and on-state certificate domains. Preparation
and compilation warmups are excluded. Construction, optimizer transfer,
callbacks and optimization consume the allowance. Extraction/replay overrun
is retained separately. Profile order reverses in repeated confirmation.

Confirmation, commit c4397c6909516c9ebd26f368323296adcf357f4d,
https://github.com/monochromatti/open-shop/actions/runs/37733769105:
Baseline, limited cuts without bilinear projection LPs, wider cuts without
projection LPs. Six cases, 300 seconds, two repetitions with reversed order.
The refinement found the limited budget better at 2 and 24 hours, but the wider
budget better at 6 hours. The no-projection combinations improved the 2- and
6-hour bounds and small-case certificate time, with unchanged 24-hour bounds.
Retain ordinary OBBT and its generalized variable bounds in both finalists.

Primary speed evidence: first observed globally valid upper bound reaching
specified gaps relative to one common audited objective. Unreached thresholds
are right censored. Final bound may give a conservative attainment time in
absence of an event. Pooled-objective times are retrospective and assume the
best audited schedule in the paired job is already known; they do not measure
incumbent discovery. Small-case final certificate time includes that discovery.

Screening selects candidates for paired 300-second, two-repetition confirmation,
not a production winner. Compare proof speed at common attained thresholds,
all four operating final gaps and small-case certificate cost. Reject errors,
invalid starts, failed schedule/replay audits, unusable global bounds and any
bound excluding an independently audited matched schedule. Keep one production
policy only if confirmation supports the tradeoff; archive experimental modes.

Exact original nonlinear equations remain. Certificates concern the declared
discrete hydro model under numerical solver tolerances. These experiments do
not establish full SHOP parity or continuous-time global optimality.
