# Object-level extensions support Newton and optimizer trial points. Strict
# public table APIs and schedule domain audits determine physical admissibility.
head(r::Reservoir, v; extrapolation = :linear) =
    r.level_curve===nothing ? r.z0+r.slope*v+r.curvature*v^2 :
    table_value(r.level_curve, v; extrapolation)
efficiency(g::Generator, q, h; extrapolation = :linear) =
    g.turbine_table===nothing ?
    g.efficiency-g.qcurvature*((q-g.qbest)/g.qbest)^2-g.hcurvature*((h-g.hbest)/g.hbest)^2 :
    turbine_efficiency(g.turbine_table, q, h; extrapolation)
generator_efficiency(g::Generator, p; extrapolation = :linear) =
    g.generator_efficiency_curve===nothing ? one(p) :
    table_value(g.generator_efficiency_curve, p; extrapolation)
function power(g::Generator, q, h; extrapolation = :linear)
    q==0 && return zero(q*h)
    shaft=0.00981*q*h*efficiency(g, q, h; extrapolation)
    g.generator_efficiency_curve===nothing ? shaft :
    electrical_power(g.generator_efficiency_curve, shaft; extrapolation)
end
tailwater(p::Plant, q; extrapolation = :linear) =
    p.tailwater_curve===nothing ? zero(q) : table_value(p.tailwater_curve, q; extrapolation)

"""Turbine outlet reference; discharged water still enters the plant's target."""
outlet_head(p::Plant, h) = p.outlet_head_floor===nothing ? h : max(h, p.outlet_head_floor)

"""Validate curve semantics and coverage without modifying supplied data."""
function validate_curve_data(s::HydroSystem)
    for r in s.reservoirs
        c=r.level_curve
        c===nothing && continue
        first(c.x)<=r.vmin<=r.vmax<=last(c.x) ||
            throw(ArgumentError("$(r.name): level table must cover storage bounds"))
        all(diff(c.y) .> 0) || throw(
            ArgumentError("$(r.name): reservoir level must strictly increase with volume"),
        )
    end
    for g in s.generators
        c=g.turbine_table
        if c!==nothing
            c.head_extrapolation==:linear ||
                first(c.heads)<=g.hmin<=g.hmax<=last(c.heads) ||
                throw(ArgumentError("$(g.name): turbine table must cover operating head"))
            first(c.discharge)<=g.qmin<=g.qmax<=last(c.discharge) ||
                throw(ArgumentError("$(g.name): turbine table must cover discharge bounds"))
            all(c.qmin .<= g.qmax) && all(c.qmax .>= g.qmin) || throw(
                ArgumentError(
                    "$(g.name): turbine and static discharge limits do not overlap",
                ),
            )
        end
        c=g.generator_efficiency_curve
        if c!==nothing
            first(c.x)==0 && last(c.x)>=g.pmax || throw(
                ArgumentError(
                    "$(g.name): electrical efficiency table must cover zero to pmax",
                ),
            )
            validate_electrical_curve(c)
        end
    end
    for p in s.plants
        c=p.tailwater_curve
        c===nothing && continue
        total=sum(g.qmax for g in s.generators if g.plant==p.name; init = 0.0)
        first(c.x)==0 && last(c.x)>=total || throw(
            ArgumentError(
                "$(p.name): tailwater table must cover zero to aggregate discharge",
            ),
        )
        all(c.y .>= 0) && all(diff(c.y) .>= 0) || throw(
            ArgumentError(
                "$(p.name): added tailwater head must be nonnegative and nondecreasing",
            ),
        )
    end
    for r in s.rivers
        c=r.discharge_curve
        c===nothing && continue
        all(c.y .>= 0) && all(diff(c.y) .>= 0) || throw(
            ArgumentError(
                "$(r.name): discharge must be nonnegative and nondecreasing with level",
            ),
        )
        source=findfirst(x->x.name==r.source, s.reservoirs)
        source===nothing &&
            throw(ArgumentError("$(r.name): discharge curve requires a reservoir source"))
        res=s.reservoirs[source]
        first(c.x)<=head(res, res.vmin)<=head(res, res.vmax)<=last(c.x) || throw(
            ArgumentError("$(r.name): discharge table must cover source reservoir levels"),
        )
    end
    true
end

"""Tabulate a synthetic analytic unit without claiming measured calibration.
The efficiency table includes q=0 as a numerical extension. Head-dependent
limits can optionally narrow the low-head envelope by a fraction of capacity.
"""
function synthetic_turbine_table(
    g::Generator;
    head_points = 5,
    discharge_points = 9,
    head_derating = 0.05,
    interpolation = :pchip_discharge,
)
    0<=head_derating<1 || throw(ArgumentError("head derating must lie in [0,1)"))
    g.hmax>g.hmin ||
        throw(ArgumentError("synthetic turbine sampling requires a nonzero head range"))
    hh=collect(range(g.hmin, g.hmax; length = head_points))
    qq=collect(range(0.0, g.qmax; length = discharge_points))
    ee=[
        g.efficiency-g.qcurvature*((q-g.qbest)/g.qbest)^2-g.hcurvature*(
            (h-g.hbest)/g.hbest
        )^2 for q in qq, h in hh
    ]
    lo=fill(g.qmin, length(hh))
    hi=[g.qmax*(1-head_derating*(g.hmax-h)/(g.hmax-g.hmin)) for h in hh]
    TurbineTable(hh, qq, ee, lo, hi; interpolation)
end
