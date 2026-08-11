#!/usr/bin/env julia
# Render a shear-thinning (Carreau) drop-impact animation as an MP4.
# Reuses the same rasterizer as run_animation.jl; only rheology differs.
#
# Usage:
#   julia --project=.. scripts/run_animation_st.jl [output.mp4]
#
# Requires: ffmpeg on PATH.

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using DropSolver
using Printf
using Logging

# ---------------------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------------------
const M      = 90
const Oh     = 0.3038
const Bo     = 1/53.9
const We     = 0.5
const v0     = -sqrt(We)

# Carreau parameters.
# Only Γ_2 = 1783566/385 is known analytically (inviscid limit, l=2).
# Set Gamma[k]=0 for k>1 to use only the validated value.
# ε_ST=0.08 keeps the bracket (1 - ε_ST·Λ²·Γ₂·Ȧ₂²) positive up to Ȧ₂≈0.92
# even at the stronger We=0.5 impact (v₀≈−0.707).
const GAMMA2 = 1783566.0 / 385.0
const EPS_ST   = 0.08
const LAMBDA_C = 0.02

const T_END      = 8.0
const SAVE_EVERY = 0.025

# ---------------------------------------------------------------------------
# Video parameters
# ---------------------------------------------------------------------------
const FPS     = 40
const WIDTH   = 720
const HEIGHT  = 720
const MP4_OUT = length(ARGS) >= 1 ? ARGS[1] :
                joinpath(@__DIR__, "..", "outputs", "figures", "impact_st.mp4")

# ---------------------------------------------------------------------------
# Colour palette — orange tones to distinguish from the OB (blue) animation
# ---------------------------------------------------------------------------
const COL_BG      = (240, 244, 248)
const COL_DROP    = (210,  90,  30)   # burnt orange
const COL_OUTLINE = (150,  50,  10)   # darker orange
const COL_FLOOR   = (120, 120, 130)
const COL_AXIS    = (180, 180, 190)

# ---------------------------------------------------------------------------
# Rasterise one frame
# ---------------------------------------------------------------------------
function render_frame(state::DropState, cfg::SimConstants,
                      x_range::Tuple{Float64,Float64},
                      z_range::Tuple{Float64,Float64})
    W, H = WIDTH, HEIGHT
    buf  = fill(UInt8(COL_BG[1]), 3 * W * H)

    function px(x, z)
        col = round(Int, (x - x_range[1]) / (x_range[2] - x_range[1]) * (W - 1)) + 1
        row = round(Int, (z_range[2] - z) / (z_range[2] - z_range[1]) * (H - 1)) + 1
        return clamp(col, 1, W), clamp(row, 1, H)
    end

    set_pixel!(col, row, rgb) = begin
        1 ≤ col ≤ W && 1 ≤ row ≤ H || return
        idx = 3 * ((row - 1) * W + (col - 1)) + 1
        buf[idx]   = UInt8(rgb[1])
        buf[idx+1] = UInt8(rgb[2])
        buf[idx+2] = UInt8(rgb[3])
    end

    function draw_line!(c1, r1, c2, r2, rgb; thick=1)
        dc, dr = c2 - c1, r2 - r1
        steps  = max(abs(dc), abs(dr), 1)
        for k in 0:steps
            c = round(Int, c1 + dc * k / steps)
            r = round(Int, r1 + dr * k / steps)
            for dc2 in -thick:thick, dr2 in -thick:thick
                set_pixel!(c + dc2, r + dr2, rgb)
            end
        end
    end

    c1, r1 = px(x_range[1], 0.0); c2, r2 = px(x_range[2], 0.0)
    draw_line!(c1, r1, c2, r2, COL_FLOOR; thick=2)
    ca1, ra1 = px(0.0, z_range[1]); ca2, ra2 = px(0.0, z_range[2])
    draw_line!(ca1, ra1, ca2, ra2, COL_AXIS; thick=1)

    xs, zs = drop_profile(state, cfg; n_theta=300)
    n = length(xs)

    for row in 1:H
        z_world = z_range[2] - (row - 1) / (H - 1) * (z_range[2] - z_range[1])
        x_cross = Float64[]
        for i in 1:n-1
            z1, z2 = zs[i], zs[i+1]
            if (z1 - z_world) * (z2 - z_world) <= 0 && z1 != z2
                t = (z_world - z1) / (z2 - z1)
                push!(x_cross, xs[i] + t * (xs[i+1] - xs[i]))
            end
        end
        isempty(x_cross) && continue
        x_max = maximum(x_cross)
        if x_max > 0
            c_left,  _ = px(-x_max, z_world)
            c_right, _ = px( x_max, z_world)
            for c in max(1, c_left):min(W, c_right)
                set_pixel!(c, row, COL_DROP)
            end
        end
    end

    for i in 1:n-1
        c1, r1 = px( xs[i],   zs[i]);   c2, r2 = px( xs[i+1], zs[i+1])
        draw_line!(c1, r1, c2, r2, COL_OUTLINE; thick=1)
        c1m, _ = px(-xs[i],   zs[i]);   c2m, _ = px(-xs[i+1], zs[i+1])
        draw_line!(c1m, r1, c2m, r2, COL_OUTLINE; thick=1)
    end

    return buf
end

# ---------------------------------------------------------------------------
# Run simulation
# ---------------------------------------------------------------------------
println("Running Carreau simulation: Oh=$Oh, Bo=$Bo, v0=$(round(v0;digits=3)), " *
        "ε_ST=$EPS_ST, λ_c=$LAMBDA_C, M=$M ...")

dt_max    = make_dt_max(M)
theta_vec = make_theta_vec(M)
precomp   = precompute_integrals(NaN, M)[1]
cfg       = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, dt_max)
ob        = OBParams()

Gamma_vec = [GAMMA2; zeros(M - 2)]   # Γ_2 known; set higher modes to 0
st        = STParams(EPS_ST, LAMBDA_C, Gamma_vec)

init = DropState(M)
init.z = 1.1; init.v = v0; init.dt = dt_max; init.cp = 0

times, states = with_logger(NullLogger()) do
    solve_drop!(cfg, ob, deepcopy(init); st=st, t_end=T_END, save_every=SAVE_EVERY)
end
println("  $(length(states)) frames saved.")

# ---------------------------------------------------------------------------
# Bounding box
# ---------------------------------------------------------------------------
all_xs = Float64[]; all_zs = Float64[]
for s in states
    xs, zs = drop_profile(s, cfg; n_theta=100)
    append!(all_xs, xs); append!(all_xs, -xs)
    append!(all_zs, zs)
end
push!(all_zs, 0.0)
x_pad = 0.3 * (maximum(all_xs) - minimum(all_xs) + 0.1)
z_pad = 0.2 * (maximum(all_zs) - minimum(all_zs) + 0.1)
x_range = (minimum(all_xs) - x_pad, maximum(all_xs) + x_pad)
z_range = (max(-0.1, minimum(all_zs) - z_pad), maximum(all_zs) + z_pad)

# ---------------------------------------------------------------------------
# Encode MP4
# ---------------------------------------------------------------------------
mkpath(dirname(MP4_OUT))

ffmpeg_cmd = `ffmpeg -y -f rawvideo -pixel_format rgb24
    -video_size $(WIDTH)x$(HEIGHT) -framerate $FPS
    -i pipe:0
    -c:v libx264 -preset fast -crf 22 -pix_fmt yuv420p
    $MP4_OUT`

println("Encoding $(length(states)) frames at $(FPS) fps → $MP4_OUT")
proc = open(ffmpeg_cmd, "w")
for (i, state) in enumerate(states)
    buf = render_frame(state, cfg, x_range, z_range)
    write(proc, buf)
    i % 50 == 0 && @printf("  frame %d / %d  (t=%.2f)\n", i, length(states), times[i])
end
close(proc)
println("Done: $MP4_OUT")
