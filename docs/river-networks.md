# River networks

A reach carries water from its upstream release to its downstream arrival.
Reservoir outlets, plant discharges, directed tunnel outfalls, upstream reaches
and natural inflow can supply it. Water in transit belongs to the river inventory
until it arrives; it is not credited to downstream storage early.

## Connections

A plant or tunnel sets `discharge_river` to a receiving reach's name. Its
`target` remains the hydraulic head reference. With `discharge_river` set,
water leaves the source and enters that reach, instead of entering `target`
immediately. The reach's eventual target determines where the water arrives.
A reservoir can therefore determine turbine tailwater while receiving its
water later, through the river network.

A routed tunnel is an outfall in its declared `source` to `target` direction:
its discharge is nonnegative. Swap the hydraulic endpoints to describe the
opposite outfall direction. Other pressurized tunnels retain signed flow.
A reverse solution at an outfall fails simulation and validation; it is never
converted to a positive release. The global model uses its direct quadratic
loss equation without an additional direction binary.

A river can set `target` to another river's name. The receiving reach uses
`source=:auto` and `law=:junction`. Any number of rivers, plants and tunnel
outfalls can feed the same receiving reach. Set `inflow` on a reach to add
natural water at its top; an inflow-only reach also uses `source=:auto`.
`OperationalSeries(attribute=:inflow, ...)` supplies time-dependent natural
inflow. Input changes must lie on scheduling edges.

`normalize_river_connections(system)` compiles direct connections from Julia
objects; `readcase` and `case_from_dict` do this automatically. Inferred
confluences are internal zero-storage nodes and disappear from `case_dict`.
Users do not need to declare a RiverJunction. Explicit junctions remain
supported for existing inputs.

One source stream has one water destination. Multiple controlled reservoir
outlets provide independent branches that can merge farther downstream.
Allocation of a single stream among several river branches requires an explicit
physical diversion model and is not inferred by this profile.

## Conservation

Let `L[d,t]` be the local injection into reach `d`, in m³/s, and `A[e,t]`
be an upstream reach's arrival volume, in Mm³. At a confluence,

```math
Q^{release}_{d,t}=L_{d,t}+\frac{\sum_{e\to d} A_{e,t}}{0.0036\Delta t_t},
\qquad
L_{d,t}=I_{d,t}+\sum_{g\to d}q_{g,t}+\sum_{e\to d}Q^{tunnel}_{e,t}.
```

A reservoir-source reach instead receives its reservoir outlet flow plus
natural inflow. The outlet law applies to `Qrelease - inflow`; only that
withdrawal depletes the source reservoir. For a controlled outlet with total
reach capacity `C` and natural inflow `I`, withdrawal is `(C-I)*gate`.
Orifice, weir and supplied discharge laws determine withdrawal independently
of natural inflow. Capacity, release ramps and release observations refer to
the total upstream reach flow. Natural inflow must not exceed reach capacity.

Each reach satisfies

```math
W_{d,t+1}=W_{d,t}+0.0036\Delta t_t Q^{release}_{d,t}-A_{d,t}.
```

Summing reservoir and reach inventories cancels internal transfers, including
plant and tunnel outfalls. Only reservoir/river natural inflow and boundary
exchange change the combined inventory. The optimization and chronological
simulation retain the same destinations and source equations.

## Travel time

`deterministic_delay` shifts each original release pulse without dispersion.
An all-deterministic network preserves those pulses through every confluence,
including zero delays and arrivals within the current interval. Hydraulic
storage and current arrivals are solved together. Capacity and pointwise
arrival constraints use pulse breakpoints, not just interval averages.

Alternatively, `curves` supplies normalized delay-bin probabilities. One curve
is independent of flow. Multiple curves have strictly increasing reference
flows covering zero through capacity; intermediate references are allowed.
At release flow `q`, only the neighboring reference distributions are blended.
The chosen distribution stays with the original cohort, including on restart.
The supplied probabilities must sum to one; invalid inputs are rejected.

OpenSHOP selects the distribution using contemporaneous release flow. SHOP
selects it using flow from the previous optimization iteration, so identical
input curves do not establish identical numerical schedules.

For multiple references the global model uses one shared SOS2 flow coordinate
per release; arrival, pointwise observations and terminal inventory all use
that same distribution. The released transfer remains nonlinear. Single-curve
transfers are linear; two references give a quadratic expression. Local
nonlinear dispatch evaluates the exact neighboring-reference polynomial and
its analytic derivatives. This polynomial is continuous but its derivative
can change at an interior reference flow.

Numerical distributed routing accumulates only transfers whose delay support
overlaps the arrival grid; it does not allocate a dense tensor for replay.
Optimization still compiles transfer coefficients on its decision grid.

A network containing distributed reaches mixes upstream arrivals over the
scheduling intervals. This also applies to deterministic reaches within that
network. Transport refinement checks convergence of cumulative arrivals and
inventory; it is not a continuous-time global certificate. `transport_audit`
reports exhaustion of its refinement limit explicitly.

## History and replay

`history_grid` and `history_release` describe water entering each reach before
the horizon. Include historical plant, tunnel and natural contributions in that
reach's upstream history. Each history ends at the horizon start. Upstream
historical arrivals propagate into downstream reaches only after that start;
pre-horizon transfers are already represented by downstream history.

`restart_case` retains original executed cohorts, delayed inventory and the
previous control-window average used by release ramps. Simulation refines
hydraulics and transport while holding the original controls. `validate` checks
conservation, outfall direction, operating laws and the scheduling-grid routing;
`replay_audit` checks the finer trajectory separately.

The [mixed-source example](../examples/mixed_source_network.jl) shows direct
connections and a globally solved two-plant case. Tests also combine a
pressurized loop, reservoir-release branches, successive confluences, routed
tunnel discharge, operating changes and unequal delays on an irregular grid.
