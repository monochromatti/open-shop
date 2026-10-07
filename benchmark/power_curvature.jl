module PowerCurvature
import OpenSHOP
export power_curvature_bounds

# Coefficients of exact shaft power in normalized discharge/head coordinates.
# Original PCHIP slopes and secant extensions are preserved.
function _power_coefficients(table, qlo, qhi, hlo, hhi, electrical_max)
    a=OpenSHOP._global_tensor_polynomial(table,qlo,qhi,hlo)
    z=OpenSHOP._global_tensor_polynomial(table,qlo,qhi,hhi)
    d=z.-a
    dq=qhi-qlo;dh=hhi-hlo
    C=zeros(5,3)
    for k in 1:4
        for (row,flow) in ((k,qlo),(k+1,dq))
            C[row,1]+=flow*hlo*a[k]
            C[row,2]+=flow*(hlo*d[k]+dh*a[k])
            C[row,3]+=flow*dh*d[k]
        end
    end
    C .*= 0.00981*electrical_max
    C
end

function _shift_axis!(target,source,origin,scale)
    degree=length(source)-1
    for k in 0:degree
        target[k+1]=scale^k*sum(
            binomial(j,k)*source[j+1]*origin^(j-k) for j in k:degree)
    end
    target
end

# Restrict the coefficient polynomial to a subbox and convert to Bernstein.
# Its smallest Bernstein coefficient is a lower enclosure in exact arithmetic.
function _box_lower(C,qlo=0.0,qhi=1.0,hlo=0.0,hhi=1.0)
    all(iszero,C) && return 0.0
    nq,nh=size(C)
    Q=similar(C);D=similar(C)
    for j in 1:nh
        _shift_axis!(view(Q,:,j),view(C,:,j),qlo,qhi-qlo)
    end
    for i in 1:nq
        _shift_axis!(view(D,i,:),view(Q,i,:),hlo,hhi-hlo)
    end
    lower=Inf
    for i in 0:(nq-1),j in 0:(nh-1)
        coefficient=sum(D[k+1,l+1]*binomial(i,k)/binomial(nq-1,k)*
            binomial(j,l)/binomial(nh-1,l) for k in 0:i,l in 0:j)
        lower=min(lower,coefficient)
    end
    # Account for cancellation while shifting, including off-state extensions.
    lower-1e-9*max(1.0,sum(abs,C),sum(abs,Q),sum(abs,D))
end

_cells(nodes)=length(nodes)==1 ? [(only(nodes),only(nodes))] :
    collect(zip(nodes[1:end-1],nodes[2:end]))
_portion(lo,hi,box)=(max(lo,box[1]),min(hi,box[2]))
_normalized(lo,hi,portion)=lo==hi ? (0.0,0.0) :
    ((portion[1]-lo)/(hi-lo),(portion[2]-lo)/(hi-lo))
function _interpolates(lo,hi,box,portion)
    portion[1]<portion[2] ||
        (box[1]==box[2] && lo<box[1]<hi)
end

"""Enclose negative curvature of exact composed turbine power.

The interpolation error along one axis depends on curvature over its full
coordinate cell, even when the on-state domain intersects only part of it.
Only the other axis may be restricted to its physical on-state domain.
Singleton coordinates need no interpolation-error correction on that axis.
"""
function power_curvature_bounds(table,qbox,hbox,qnodes,hnodes,electrical_max)
    all(isfinite,(qbox...,hbox...,electrical_max)) &&
        qbox[1]<=qbox[2] && hbox[1]<=hbox[2] && electrical_max>=0 ||
        throw(ArgumentError("finite ordered curvature domains required"))
    length(qnodes)>=1 && length(hnodes)>=1 &&
        all(isfinite,qnodes) && all(isfinite,hnodes) &&
        all(diff(qnodes).>0) && all(diff(hnodes).>0) &&
        first(qnodes)<=qbox[1]<=qbox[2]<=last(qnodes) &&
        first(hnodes)<=hbox[1]<=hbox[2]<=last(hnodes) ||
        throw(ArgumentError("coordinate cells must cover the on-state box"))
    qminimum=0.0;hminimum=0.0
    for (ql,qr) in _cells(qnodes),(hl,hr) in _cells(hnodes)
        qp=_portion(ql,qr,qbox);hp=_portion(hl,hr,hbox)
        qp[1]<=qp[2] && hp[1]<=hp[2] || continue
        C=_power_coefficients(table,ql,qr,hl,hr,electrical_max)
        all(isfinite,C) || throw(ArgumentError("power-curvature coefficients overflow"))
        if ql<qr && _interpolates(ql,qr,qbox,qp)
            Dq=[k*(k-1)*C[k+1,l+1]/(qr-ql)^2 for k in 2:4,l in 0:2]
            hpart=_normalized(hl,hr,hp)
            # Do not clip discharge to qp: chord endpoints are ql and qr.
            qminimum=min(qminimum,_box_lower(Dq,0.0,1.0,hpart...))
        end
        if hl<hr && _interpolates(hl,hr,hbox,hp)
            Dh=reshape(2*C[:,3]/(hr-hl)^2,5,1)
            qpart=_normalized(ql,qr,qp)
            # Head curvature is constant in head within each source cell.
            hminimum=min(hminimum,_box_lower(Dh,qpart...,0.0,1.0))
        end
    end
    result=(discharge=max(0.0,-0.5*qminimum),head=max(0.0,-0.5*hminimum))
    all(isfinite,result) || throw(ArgumentError("power curvature enclosure overflow"))
    result
end
end
