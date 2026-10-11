# A delta-level reach is a zero-storage hydraulic connection. Its signed flow
# enters both reservoir balances immediately; it carries no delayed inventory.
is_channel(r::River) = r.law==:delta
channel_curve(r::River) = r.delta_level_ref_up_flow_curve===nothing ?
    r.delta_level_ref_down_flow_curve : r.delta_level_ref_up_flow_curve
channel_reference(r::River, upstream, downstream) =
    r.delta_level_ref_up_flow_curve===nothing ? downstream : upstream
channel_flow(r::River, upstream, downstream) = level_flow_value(
    channel_curve(r), upstream-downstream, channel_reference(r,upstream,downstream);
    extrapolation=:linear)

function validate_channel(r::River, s::HydroSystem)
    supplied=(r.delta_level_ref_up_flow_curve!==nothing)+
             (r.delta_level_ref_down_flow_curve!==nothing)
    if !is_channel(r)
        supplied==0 || throw(ArgumentError("$(r.name): delta-level curves require law=delta"))
        return
    end
    isfinite(r.capacity) && r.capacity>0 || throw(ArgumentError("invalid channel capacity"))
    supplied==1 || throw(ArgumentError("$(r.name): choose exactly one delta-level reference"))
    any(x.name==r.source for x in s.reservoirs) &&
    any(x.name==r.target for x in s.reservoirs) && r.source!=r.target ||
        throw(ArgumentError("$(r.name): delta-level channel needs two distinct reservoirs"))
    r.deterministic_delay==0.0 && isempty(r.curves) ||
        throw(ArgumentError("$(r.name): signed channels require zero delay and no transport curves"))
    r.inflow==0 && r.min_arrival==0 && r.gate_min==0 &&
    r.arrival_policy==:interval_average && isempty(r.arrival_window_grid) &&
    r.discharge_curve===nothing && r.coefficient==0 && r.crest==0 && !r.allow_dry ||
        throw(ArgumentError("$(r.name): delta-level law cannot be combined with outlet, inflow or arrival laws"))
    all(iszero,r.history_release) ||
        throw(ArgumentError("$(r.name): instantaneous channels have no historical inventory"))
end

# Sign-aware domains are shared by routing and independent physics. Ordinary
# downstream reaches retain nonnegative release and transport-cohort semantics.
river_flow_lower(r, capacity, inflow=0.0) = is_channel(r) ? -capacity : inflow
river_flow_valid(r, q; tolerance=1e-10) = isfinite(q) &&
    river_flow_lower(r,r.capacity)-tolerance<=q<=r.capacity+tolerance

function _check_river_releases(r, grid, q)
    RiverRouting.check_grid(grid)
    length(q)==length(grid)-1 || throw(DimensionMismatch("river release grid"))
    if is_channel(r)
        all(x->river_flow_valid(r,x),q) ||
            throw(ArgumentError("$(r.name): signed release exceeds channel capacity"))
    else
        RiverRouting.check_releases(grid,q)
    end
    nothing
end
