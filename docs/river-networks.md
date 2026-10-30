# River networks

A delayed reach carries water from its upstream release to its downstream arrival. Instantaneous reservoir-to-reservoir channels instead use signed [delta-level flow curves](input.md).
Reservoir outlets, plant discharges, directed tunnel outfalls, upstream reaches
and natural inflow can supply it. Water in transit belongs to the river inventory
until it arrives; it is not credited to downstream storage early.

Delayed reaches assume physically valid downstream connections and calibrated
entry responses. Routing does not check gravitational elevation compatibility
or calculate backwater from changing downstream levels. A one-level entry
curve cannot describe independently varying downstream counterpressure. See
[model validity](model-validity.md) for these assumptions and time-grid limits.

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

`ScheduleCase` compiles direct connections from native Julia objects;
`readcase` and `case_from_dict` use the same construction. Inferred
confluences are internal zero-storage nodes and disappear from `case_dict`.
Users do not need to declare a RiverJunction. Explicit junctions remain
supported for existing inputs.

To split a stream, give its incoming reach a target label shared by the outgoing
reaches' sources. Each outgoing reach declares its prescribed `allocation`:

```julia
River(name=:Incoming, source=:Lake, target=:Fork, law=:controlled,
    capacity=10.0, curves=[], deterministic_delay=0.25, water_value=0.0)
River(name=:Left, source=:Fork, target=:Receiver, allocation=0.4,
    capacity=6.0, curves=[], deterministic_delay=0.5, water_value=0.0)
River(name=:Right, source=:Fork, target=:Receiver, allocation=0.6,
    capacity=8.0, curves=[], deterministic_delay=1.0, water_value=0.0)
```

`Fork` is a connection label, not an object to construct. It must differ from
every object name. The same form describes a merge or a node with several
incoming and outgoing reaches. `case_dict` retains these labels and fractions.
Existing explicit `RiverJunction` inputs also support prescribed allocations.

Fractions are finite, lie in `[0,1]`, and sum to one at each node. They default
to one, so a single outgoing reach needs no extra data. Splits must supply all
shares; OpenSHOP does not infer missing shares or normalize them. Fractions are
constant over the horizon, with a `1e-12` absolute and relative sum tolerance.
They describe a supplied diversion rule, not a decision about gate operation.
Controllable and time-dependent allocation are outside this profile. Reservoir
outlets remain independent decisions and use `allocation=1`.

Local natural inflow and plant/tunnel outfalls enter their named branch after
allocation. To distribute a common supply, feed a reach targeting the fork.
A river targeting a reach directly still requires that receiving reach to use
`source=:auto`; it cannot also have an explicit fork source. Use the node label
to add a tributary before the fork, or a downstream direct merge to add one
after it. Names beginning with `__river_merge__` or `__river_node__` are reserved.

## Conservation

Let `L[d,t]` be the local injection into reach `d`, in m³/s, and `A[e,t]`
be an upstream reach's arrival volume, in Mm³. At a zero-storage river node,

```math
Q^{release}_{d,t}=L_{d,t}+\alpha_d\frac{\sum_{e\to d} A_{e,t}}{0.0036\Delta t_t},
\qquad
L_{d,t}=I_{d,t}+\sum_{g\to d}q_{g,t}+\sum_{e\to d}Q^{tunnel}_{e,t}.
```

Here `α[d]` is the outgoing reach's `allocation`, and `e→d` denotes reaches
arriving at its source node. At every such node, the outgoing shares sum to one.
Thus a 40% branch with 5 m³/s arriving upstream and 1 m³/s local inflow releases
3 m³/s. A zero-share branch still carries its local water. Each branch's capacity
limits its total flow; spare capacity on another branch does not change the rule.

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
An all-deterministic network preserves those pulses through every merge and split,
including zero delays and arrivals within the current interval. Hydraulic
storage and current arrivals are solved together. Capacity and pointwise
arrival constraints use pulse breakpoints, not just interval averages.
Capacity checks in hydraulic replay and fixed-delay routing allow `1e-6` m³/s
of numerical roundoff. They retain the computed flow and report its excess;
they do not clip the water. Distributed reference-flow coverage remains a
separate constitutive-domain check.

Alternatively, `curves` supplies normalized delay-bin probabilities. One curve
is independent of flow. Multiple curves have strictly increasing reference
flows covering zero through capacity; intermediate references are allowed.
At release flow `q`, only the neighboring reference distributions are blended.
The chosen distribution stays with the original cohort, including on restart.
The supplied probabilities must sum to one; invalid inputs are rejected.

For multiple references the global model uses one shared SOS2 flow coordinate
per release; arrival, pointwise observations and terminal inventory all use
that same distribution. The released transfer remains nonlinear. Single-curve
transfers are linear; two references give a quadratic expression. Local
nonlinear dispatch evaluates the exact neighboring-reference polynomial and
its analytic derivatives. This polynomial is continuous but its derivative
can change at an interior reference flow.

Numerical distributed routing accumulates only transfers whose delay support
overlaps the arrival grid; it does not allocate a dense tensor for replay.
Optimization compiles transfer coefficients between each reach's entry-cohort
grid and the river transport grid.

A network containing distributed reaches mixes upstream arrivals over the
transport cells. This also applies to deterministic reaches within that
network. The default uses four transport cells per control interval when the
network contains distributed confluences. An optional `routing_grid` on
`ScheduleCase` overrides this mesh; it must span the horizon and contain every
control boundary. Source releases, unit decisions, gate positions, operating
ramps and observation windows retain their original control grid. Refining
transport adds no reservoir states or control decisions. Downstream cell flows
are determined by conservation, and their interval averages remain the reported
river releases. Restart history preserves the actual cell cohorts.

A confluence's incoming wave must respect the outgoing reach's capacity before
cell averaging. The optimizers enforce this at every upstream arrival knot and
both sides of jumps, where a piecewise affine wave reaches its extrema. A mean
flow below capacity can still hide a peak above it. Independent transport
acceptance checks these peaks again, as well as minimum arrival rates and
declared arrival windows, on the finer routed cohorts.

Transport refinement checks convergence of cumulative arrivals and
inventory; it is not a continuous-time global certificate. `transport_audit`
uses at most six successive bisections by default, retaining the same available
resolution per control interval as the horizon grows. An optional
`max_internal_intervals` caps the total cell count. The check reports
nonconvergence when either limit stops refinement before the volume tolerance
is met. A larger tolerance permits a larger estimated discretization error.

Inventory is evaluated directly from the survival probabilities of the original
entry cohorts, independently of arrival-volume conservation. A cohort contributes
only between its release and its last possible arrival. Fixed-delay volumes use
interval overlaps without allocating dense matrices during numerical replay.

## History and replay

`history_grid` and `history_release` describe water entering each reach before
the horizon. Include historical plant, tunnel and natural contributions in that
reach's upstream history. Each history ends at the horizon start. Upstream
historical arrivals propagate into downstream reaches only after that start;
pre-horizon transfers are already represented by downstream history.
Incoming historical water arriving at a fork after the horizon start receives
the prescribed shares. A branch's own historical cohorts are already allocated
and are never multiplied by its share again.

`restart_case` retains original executed cohorts, delayed inventory and the
previous control-window average used by release ramps. Simulation refines
hydraulics and transport while holding the original controls. `validate` checks
conservation, outfall direction, operating laws and the scheduling-grid routing;
`replay_audit` checks the finer trajectory separately.

The [mixed-source example](../examples/mixed_source_network.jl) shows direct
connections and a globally solved two-plant case. The [river-split example](../examples/river_splits.jl)
optimizes a unit feeding unequal branches with a binding branch capacity.
Tests also combine a
pressurized loop, reservoir-release branches, successive confluences, routed
tunnel discharge, operating changes and unequal delays on an irregular grid.
