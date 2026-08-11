# Resolution of the Oh_eff coupling quadrature.
#
# The default `n_x` scales with `M` (see the STExactParams docstring). Nothing
# else in the suite would notice if it stopped: every other `Oh_eff` test runs
# at small `M`, where the old fixed 30-node rule was already accurate, and the
# quantity has no closed form to compare against at large `M`. These tests
# compare against a converged reference rule instead.

"""
Worst-mode relative difference in `Oh_eff` between an `(n_r, n_x)` rule and a
converged one. The reference must be far finer in `n_x` than anything being
measured: an under-resolved reference floors the reported error and hides the
very convergence being tested (a 300-node reference makes 182 and 270 nodes
look equally good at `M = 90`, when they differ by 20x against a 500-node one).
"""
function _q_err(M::Int, adot::Vector{Float64}, n_r::Int, n_x::Int;
                ref_r::Int = 40, ref_x::Int = 600)
    mk(nr, nx) = STExactParams(M, 0.5, 20.0, 0.75, 1.0; viscous = :lamb,
                               eta_inf_ratio = 0.004, n_r = nr, n_x = nx)
    ref = oh_eff_all_coupled(mk(ref_r, ref_x), 0.5, adot)
    got = oh_eff_all_coupled(mk(n_r, n_x), 0.5, adot)
    maximum(abs.(got .- ref) ./ ref)
end

# A spectrum that decays like l^-1: broad enough that high modes carry real
# weight, which is what makes the angular resolution matter. A steeply decaying
# state would hide the defect by putting all the shear in low modes. This is a
# deliberately harsher test than a real contact-phase state, whose spectrum
# falls off faster (worst-mode error there is ~8e-4 at M = 90, versus ~1e-2
# here) -- the point is to bound the bad case, not the typical one.
_q_adot(M::Int) = [1.0 / (l - 1) for l in 2:M]

@testset "Oh_eff coupling quadrature" begin

    @testset "Default n_x scales with M" begin
        # The rule itself, independent of any solve.
        @test STExactParams(6,  0.5, 20.0, 0.75, 1.0; viscous = :lamb) |>
              s -> length(s.x_nodes) == 30            # floor holds at small M
        @test length(STExactParams(60, 0.5, 20.0, 0.75, 1.0; viscous = :lamb).x_nodes) == 122
        @test length(STExactParams(90, 0.5, 20.0, 0.75, 1.0; viscous = :lamb).x_nodes) == 182
        # Radial resolution is deliberately NOT scaled.
        @test length(STExactParams(90, 0.5, 20.0, 0.75, 1.0; viscous = :lamb).r_nodes) == 20
    end

    @testset "Accuracy does not degrade as M grows" begin
        # This -- not any particular error value -- is what the scaling rule
        # buys, so it is what gets asserted. Asserting an absolute threshold
        # instead would pass just as happily on a rule that was uniformly
        # accurate for the wrong reason, and would need refitting whenever the
        # reference or the test spectrum changed.
        nx(M) = length(STExactParams(M, 0.5, 20.0, 0.75, 1.0; viscous = :lamb).x_nodes)
        scaled_lo = _q_err(16, _q_adot(16), 20, nx(16))
        scaled_hi = _q_err(90, _q_adot(90), 20, nx(90))
        fixed_lo  = _q_err(16, _q_adot(16), 20, 30)
        fixed_hi  = _q_err(90, _q_adot(90), 20, 30)

        @test scaled_hi / scaled_lo < 3.0        # measured ~1.8
        @test fixed_hi / fixed_lo   > 5.0        # measured ~8.4 -- the defect
        @test scaled_hi < fixed_hi / 10          # and it is far more accurate
    end

    @testset "The fixed 30-node rule is under-resolved at production M" begin
        # Guards the reasoning, not just the outcome: if this ever passes, the
        # integrand has changed and the scaling rule needs rederiving rather
        # than being carried along out of habit.
        @test _q_err(90, _q_adot(90), 20, 30) > 1e-2
    end

    @testset "Radial rule is already converged at n_r = 20" begin
        # Justifies not scaling n_r: at fixed fine angular resolution, going
        # from 20 to 120 radial nodes moves Oh_eff by well under the angular
        # error, because Gauss nodes cluster where r^(l-2) lives.
        M = 90
        adot = _q_adot(M)
        mk(nr) = STExactParams(M, 0.5, 20.0, 0.75, 1.0; viscous = :lamb,
                               n_r = nr, n_x = 300)
        coarse = oh_eff_all_coupled(mk(20), 0.5, adot)
        fine   = oh_eff_all_coupled(mk(120), 0.5, adot)
        @test maximum(abs.(coarse .- fine) ./ fine) < 1e-3
    end
end
