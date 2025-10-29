# ================================================
# Monte Carlo estimate of effective capture radius r_eff
# for a surface-tethered sensor of contour length l.
#
# Model:
# - Discrete WLC approximated as a correlated random chain (N segments)
# - Segment length b = 2*lp (Kuhn length)
# - One end grafted at z = h0 ≥ 0, chain confined to half-space z ≥ 0
# - A conformation is "capturable" if end-point 0 < z < zcap
# - r_eff = sqrt( ⟨ x^2 + y^2 ⟩ | 0<z<zcap )
#
# Also returns: capture_fraction = (#capturable conformations) / M
#
# Notes:
# - Correlations implemented via <cosθ> = a = exp(-b/ℓp), achieved by
#   blending the previous direction with an isotropic orthogonal component.
# - If you want a freely-jointed chain, set use_WLC = false (a = 0).
#
# Units:
# - Use any consistent length unit (e.g., nm). l, lp, zcap, h0, b all in same unit.
#
# Dependencies: Base only.
# ================================================

using Random, Printf

# -----------------------------
# Utility: random unit vector on S^2
# -----------------------------
function rand_unit_vector(rng::AbstractRNG=Random.GLOBAL_RNG)
    # Sample from 3D normal and normalize (Marsaglia method would also work)
    x = randn(rng, 3)
    n = sqrt(sum(abs2, x))
    return x ./ n
end

# -----------------------------
# Utility: correlated direction update (approximate WLC step)
# Inputs:
#   prev :: Vector{Float64} (unit)
#   a    :: Float64 in [0,1], desired <cosθ> correlation
# Returns:
#   new unit direction with ⟨prev·new⟩ ≈ a
# -----------------------------
function step_direction_correlated(prev::Vector{Float64}, a::Float64, rng::AbstractRNG)
    # a=0 → isotropic; a→1 → strongly forward-peaked
    if a <= 0
        return rand_unit_vector(rng)
    elseif a >= 1
        return copy(prev)
    end

    # Draw an isotropic direction and make it orthogonal to prev
    u = rand_unit_vector(rng)
    # Orthogonal component to prev
    dotpu = prev[1]*u[1] + prev[2]*u[2] + prev[3]*u[3]
    u_perp = [u[1] - dotpu*prev[1], u[2] - dotpu*prev[2], u[3] - dotpu*prev[3]]
    nperp = sqrt(sum(abs2, u_perp))
    if nperp < 1e-12
        # Rare degenerate case: fall back to slight random tilt
        return rand_unit_vector(rng)
    end
    u_perp ./= nperp

    # Blend prev with orthonormal component so that cosθ ≈ a
    s = sqrt(max(0.0, 1.0 - a^2))
    newdir = [a*prev[1] + s*u_perp[1],
              a*prev[2] + s*u_perp[2],
              a*prev[3] + s*u_perp[3]]
    nnew = sqrt(sum(abs2, newdir))
    return newdir ./ nnew
end

# -----------------------------
# Single conformation sampler (half-space constrained)
# Inputs:
#   N     :: Int         # number of segments
#   b     :: Float64     # segment length (same units as l, lp, zcap)
#   a     :: Float64     # correlation parameter, a = exp(-b/lp) if WLC
#   h0    :: Float64     # anchor height above plane (≥ 0)
#   rng   :: RNG
# Returns:
#   endpoint R::Vector{Float64} of size 3 (x,y,z) with z ≥ 0
# -----------------------------
function sample_endpoint_halfspace(N::Int, b::Float64, a::Float64, h0::Float64, rng::AbstractRNG)
    # Start at anchor
    R = [0.0, 0.0, max(h0, 0.0)]

    # First direction (biased to half-space z≥0 by rejection)
    dir = rand_unit_vector(rng)
    tries = 0
    max_tries_init = 1_000
    while R[3] + b*dir[3] < 0.0
        if tries < max_tries_init
            dir = rand_unit_vector(rng)
            tries += 1
        else
            # Fall back: ensure positive z-component
            dir[3] = abs(dir[3])
            ndir = sqrt(sum(abs2, dir))
            dir ./= ndir
            break
        end
    end
    # Take step
    R[1] += b*dir[1]; R[2] += b*dir[2]; R[3] += b*dir[3]

    # Remaining steps
    for i in 2:N
        # Propose correlated direction
        prop = step_direction_correlated(dir, a, rng)
        # Enforce half-space by resampling direction if crossing occurs
        tries = 0
        max_tries = 1_000
        while R[3] + b*prop[3] < 0.0
            if tries < max_tries
                # Try resampling with correlation
                prop = step_direction_correlated(dir, a, rng)
            else
                # Fall back to reflection: flip the z-component
                prop[3] = -prop[3]
                # Renormalize
                nprop = sqrt(sum(abs2, prop))
                prop ./= nprop
                break
            end
            tries += 1
        end
        dir = prop
        R[1] += b*dir[1]; R[2] += b*dir[2]; R[3] += b*dir[3]
    end
    return R
end

# -----------------------------
# Monte Carlo driver to estimate r_eff
# Inputs:
#   l        :: contour length
#   lp       :: persistence length
#   zcap     :: capture window thickness above plane (0<z<zcap)
#   M        :: number of MC samples
#   h0       :: anchor height (default 0)
#   use_WLC  :: Bool, if false use FJC (a=0)
#   seed     :: Int, RNG seed
# Returns:
#   (r_eff, capture_fraction)
# -----------------------------
function estimate_reff(l::Float64; lp::Float64=1.0, zcap::Float64=1.0, M::Int=100_000,
                       h0::Float64=0.0, use_WLC::Bool=true, seed::Int=1234)

    # Input validation
    l > 0.0 || throw(ArgumentError("Contour length l must be positive, got $l"))
    lp > 0.0 || throw(ArgumentError("Persistence length lp must be positive, got $lp"))
    zcap > 0.0 || throw(ArgumentError("Capture window zcap must be positive, got $zcap"))
    M > 0 || throw(ArgumentError("Number of samples M must be positive, got $M"))
    h0 >= 0.0 || throw(ArgumentError("Anchor height h0 must be non-negative, got $h0"))

    rng = MersenneTwister(seed)

    # Discretization: Kuhn length b = 2*lp; ensure at least 1 segment
    b = 2.0*lp
    N = max(1, Int(floor(l / b)))
    # If remainder > 0, you can slightly stretch the last step; we keep uniform steps for simplicity.
    # Correlation parameter: a = exp(-b/lp) for discrete WLC; a=0 for FJC
    a = use_WLC ? exp(-b / lp) : 0.0
    # Clamp numerical issues
    a = min(max(a, 0.0), 1.0)

    r2sum = 0.0
    hits  = 0

    for m in 1:M
        R = sample_endpoint_halfspace(N, b, a, h0, rng)
        z = R[3]
        if (z > 0.0) && (z < zcap)
            r2sum += R[1]^2 + R[2]^2
            hits  += 1
        end
    end

    capture_fraction = hits / M
    r_eff = (hits > 0) ? sqrt(r2sum / hits) : 0.0
    return r_eff, capture_fraction, (N=N, b=b, a=a)
end

# -----------------------------
# Example usage
# -----------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    # Example parameters (nm):
    # PEG-like tether at room temperature: lp ~ 0.38–0.5 nm (choose 0.5 nm here)
    # Contour lengths l ~ 2, 5, 10 nm; capture window zcap ~ 1 nm
    lp = 0.5
    zcap = 1.0
    for l in (2.0, 5.0, 10.0)
        r_eff_WLC, fcap_WLC, metaW = estimate_reff(l; lp=lp, zcap=zcap, M=50_000, use_WLC=true, seed=42)
        r_eff_FJC, fcap_FJC, metaF = estimate_reff(l; lp=lp, zcap=zcap, M=50_000, use_WLC=false, seed=42)
        @printf("\n--- l = %.1f nm, lp = %.2f nm, zcap = %.1f nm ---\n", l, lp, zcap)
        @printf("WLC approx: N=%d, b=%.2f, a=%.3f  →  r_eff = %.3f nm, capture_fraction = %.4f\n",
                metaW.N, metaW.b, metaW.a, r_eff_WLC, fcap_WLC)
        @printf("FJC (a=0):  N=%d, b=%.2f, a=%.3f  →  r_eff = %.3f nm, capture_fraction = %.4f\n",
                metaF.N, metaF.b, metaF.a, r_eff_FJC, fcap_FJC)
    end

    # Once r_eff is known and sensor surface density ρ_s (nm^-2) is set,
    # the single-pass capture probability is P = 1 - exp(-π r_eff^2 ρ_s).
    # Example: pitch p = 3 nm → ρ_s = 1/p^2 = 1/9 nm^-2
    rhos = 1.0 / 9.0
    l = 5.0
    r_eff, fcap, _ = estimate_reff(l; lp=lp, zcap=zcap, M=100_000, use_WLC=true, seed=7)
    P = 1.0 - exp(-π * r_eff^2 * rhos)
    @printf("\nFor l=%.1f nm, estimated r_eff=%.3f nm → P(ρ_s=1/9 nm^-2) = %.3f (single-pass)\n", l, r_eff, P)
end