"""
    bdf_coefficients(order, dt, dt_prev) → Vector{Float64}

Return BDF coefficients `c` such that the BDF approximation is
    c[end]*y^k + c[end-1]*y^{k-1} + ... = dt * f(y^k)

order=1: c = [-1, 1]
order=2: c = [ck, bk, ak] with rk = dt/dt_prev,
         ak=(1+2rk)/(1+rk), bk=-(1+rk), ck=rk²/(1+rk)
"""
function bdf_coefficients(order::Int, dt::Float64, dt_prev::Float64)
    if order == 1
        return [-1.0, 1.0]
    elseif order == 2
        rk = dt / dt_prev
        ak = (1 + 2rk) / (1 + rk)
        bk = -(1 + rk)
        ck = rk^2 / (1 + rk)
        return [ck, bk, ak]
    else
        error("Only BDF orders 1 and 2 are supported.")
    end
end
