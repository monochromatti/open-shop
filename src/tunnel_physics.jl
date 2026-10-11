"""Head multiplier and signed-flow resistance at a prescribed tunnel opening.

With a supplied gate-loss table, its interpolated coefficient is added to the
fixed tunnel resistance. Without one, retain the existing opening-scaled law.
An opening of zero is a closed gate regardless of its tabulated coefficient.
"""
function tunnel_coefficients(t::Tunnel, opening=t.opening)
    c=t.gate_loss_factor_curve
    c===nothing ? (opening, t.resistance) :
        (1.0, t.resistance+table_value(c, opening))
end

function tunnel_loss_residual(t::Tunnel, dh, q, opening=t.opening)
    opening==0 && return q
    factor, resistance=tunnel_coefficients(t, opening)
    factor*dh-resistance*q*abs(q)
end

function tunnel_flow(t::Tunnel, dh, opening=t.opening)
    opening==0 && return zero(dh)
    factor, resistance=tunnel_coefficients(t, opening)
    sign(dh)*sqrt(factor*abs(dh)/resistance)
end

# Positive gate positions can all use the unscaled head difference. The legacy
# opening multiplier is equivalent to increasing resistance by its reciprocal.
function tunnel_effective_resistance(e::Tunnel,opening)
    opening>0 || throw(ArgumentError("effective resistance requires an open gate"))
    resistance=e.gate_loss_factor_curve===nothing ? e.resistance/opening :
        e.resistance+table_value(e.gate_loss_factor_curve,opening)
    isfinite(resistance) && resistance>0 ||
        throw(ArgumentError("$(e.name): gate resistance must be finite and positive"))
    resistance
end
