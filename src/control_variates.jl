export cv_estimate

"""
    cv_estimate(X, zero_mean_cvs...)

Control-variates estimator for E[X] using zero_mean_cvs as controls.
ASSUMES means of control variates are 0.

Inputs:
  X :: AbstractVector       (length n)
  zero_mean_cvs :: AbstractVector...  


Returns NamedTuple:
  (mean_raw, mean_cv, alpha, var_raw, var_residual, variance_reduction)
"""
function cv_estimate(X::AbstractVector{T}, zero_mean_cvs::AbstractVector{T}...) where T
    
    mean_raw = mean(X)
    var_raw  = var(X; corrected=true)

    # Fallback so I can test effect of no CVs
    if length(zero_mean_cvs) == 0
        @warn "No control variates provided, returning raw estimate." maxlog=1
        return (mean_raw=mean_raw, mean_cv=mean_raw, alpha=zeros(T,0), var_raw=var_raw,
                var_cv=var_raw, variance_reduction=1.0)
    end

    W = hcat(zero_mean_cvs...) # n x p 
    n, p = size(W)
    @assert n > p "Need more samples (n=$(n)) than control variates (p=$(p))."
    @assert length(X) == n "X and W must have same number of rows."

    # C.V. should already have zero mean, but this ensures
    # the intercept is actually 0.
    μ_estimate = mean(W; dims=1)               # p
    Z  = W .- μ_estimate                       # n×p, zero-mean by sample

    X_centered = X .- mean_raw                 # n, zero-mean by sample ← ADD THIS
    
    α = Z \ X_centered                         # Regression on centered data ← CHANGED
    resid = X_centered .- Z * α               # Residuals from centered data ← CHANGED
    mean_cv = mean_raw + mean(resid)          # Adjust back to original scale ← CHANGED
    var_cv = var(resid; corrected=true)

    return (mean_raw=mean_raw, mean_cv=mean_cv, alpha=α, var_raw=var_raw,
            var_cv=var_cv, variance_reduction=var_raw/var_cv)

end

# Assumes harmonic reference
function get_cv_estimates(X, V2, T, n_atoms)

    # Build some control variates
    # These all have zero mean by construction
    f = 3*n_atoms - 3
    μ₂ = 0.5*f*kB*T  # ⟨V2⟩
    ∂μ₂_∂T = 0.5*f*kB
    σ₂² = 0.5*f*(kB*T)^2 # var(V2)

    C1 = V2 .- μ₂
    C1_sq = C1 .^ 2
    C2 = C1_sq .- σ₂²
    C3 = C1 .^ 3
    C4 = (C1 .^ 4) .- ((3*σ₂²).*C1_sq)

    # Estimate for <X>
    res1 = cv_estimate(X, C1, C2, C3, C4)
    μX = res1.mean_cv

    # Estimate for ∂<X>/∂T
    ∂X_∂T, red = derivative_with_opt_c(X, V2, T, n_atoms, K = 2)

    # Estimate for ∂²<X>/∂T²
    dXZ_1, red = derivative_with_opt_c(X .* C1, V2, T, n_atoms, K = 3)
    dXZ = dXZ_1 - (∂μ₂_∂T * μX)
    ∂²X_∂T² = (-2*∂X_∂T/T) + (dXZ/(kB*T*T))

    return μX, ∂X_∂T, ∂²X_∂T²

end


"""
    derivative_with_opt_c(X, V2, T, n_atoms; K::Int=2)


Single control variate estimator for the derivative d⟨X⟩/dT.

Estimate:
  cov_hat = ⟨(X - c*) * (V2 - μ2)⟩_0
  deriv_hat = cov_hat / (kB * T^2)

where μ2 = (f/2) * kB * T and Var(V2) = (f/2) * (kB*T)^2 with f = 3*n_atoms - 3.

Arguments
- `X::AbstractVector{<:Real}`: your `X` samples (e.g., V - V2 or any statistic).
- `V2::AbstractVector{<:Real}`: harmonic potential energy samples.
- `T::Real`: temperature (K).
- `n_atoms::Integer`: number of atoms (used for f = 3*n_atoms).
- `K`: folds for cross-fitting (default 2).

Returns
- `(; cov_hat, cov_se, deriv_hat, deriv_se, c_folds, mu2, varV2_analytic)`

"""
function derivative_with_opt_c(X::AbstractVector{<:Real},
                               V2::AbstractVector{<:Real},
                               T::Real,
                               n_atoms::Integer; K::Int=2)

    @assert length(X) == length(V2) "X and V2 must have the same length."
    n = length(X)
    n > 0 || error("Empty inputs.")

    # Degrees of freedom
    f = 3 * n_atoms - 3

    # True harmonic mean and variance of V2
    mu2 = 0.5 * f * kB * T
    varV2_analytic = 0.5 * f * (kB * T)^2  # E[Z^2] with Z = V2 - mu2

    Z = V2 .- mu2

    # Build K folds
    idx = collect(1:n)
    Random.shuffle!(idx)
    folds = [idx[round(Int, floor((k-1)*n/K))+1 : round(Int, floor(k*n/K))] for k in 1:K]
    # guard against rounding issues
    folds[end] = idx[sum(length.(folds[1:end-1]))+1:end]

    # Cross-fit c on complement of each fold; evaluate on the fold
    y = similar(float.(X))  # per-sample integrand values (with each fold's c)
    c_folds = Float64[]
    for fold in folds
        comp = setdiff(idx, fold)
        # Numerator: E[X Z^2] estimated on complement (use mean for scale invariance)
        num = mean(X[comp] .* (Z[comp]).^2)
        # Denominator: analytic E[Z^2] = Var(V2)
        c_hat = num / varV2_analytic
        push!(c_folds, c_hat)

        # Evaluate integrand on held-out fold
        @inbounds @views y[fold] .= (X[fold] .- c_hat) .* Z[fold]
    end

    cov_hat = mean(y)
    deriv_hat = cov_hat / (kB * T^2)

    raw_varY = var(X .* Z; corrected=true)
    red = raw_varY / var(y; corrected=true)

    return deriv_hat, red
end



# is = zeros(Int, length(X))
# μX_b = zeros(n_boot)
# ∂X_∂T_b = zeros(n_boot)
# ∂2X_∂T²_b = zeros(n_boot)
# for i in 1:n_boot
#     sample!(1:length(X), is; replace = true)
#     r1 = mean(X[is])
#     r2 = cov(X[is], V2[is]) / (kB * T * T)
#     μX_b[i] = r1
#     ∂X_∂T_b[i] = r2
#     dXZ = cov(X[is].*C1[is], V2[is]) / (kB * T * T) - (∂μ₂_∂T * μX_b[i])
#     # dXV2 = cov(X[is] .* V2[is], V2[is]) / (kB * T * T)
#     # dAB = (μX_b[i] * ∂μ₂_∂T) + (μ₂ * ∂X_∂T_b[i])
#     ∂2X_∂T²_b[i] = (-2*∂X_∂T_b[i]/T) + ((dXZ/(kB*T*T)))# * (dXV2 - dAB))
# end
# μX_SE_raw = std(μX_b)
# ∂X_∂T_SE_raw = std(∂X_∂T_b)
# ∂2X_∂T²_SE_raw = std(∂2X_∂T²_b)
# println("Bootstrap SE for ⟨X⟩ (raw): $(μX_SE_raw)")
# println("Bootstrap SE for ∂⟨X⟩/∂T (raw): $(∂X_∂T_SE_raw)")
# println("Bootstrap SE for ∂²⟨X⟩/∂T² (raw): $(∂2X_∂T²_SE_raw)")

    # mu_X_raw = mean(X)
# dx_dT_raw = cov(X, V2) / (kB * T * T)
# dXV2 = cov(X .* V2, V2) / (kB * T * T)
# dAB = (mu_X_raw * ∂μ₂_∂T) + (μ₂ * dx_dT_raw)
# dXZ = cov(X.*C1, V2) / (kB * T * T) - (∂μ₂_∂T * mu_X_raw)
# println("Raw <X> estimate: $(mean(X))")
# println("Raw ∂<X>/∂T estimate: $(cov(X, V2) / (kB * T * T))")
# println("Raw ∂²<X>/∂T² estimate: $((-2*dx_dT_raw/T) + ((1/(kB*T*T)) * (dXV2 - dAB)))")
# println("Raw ∂²<X>/∂T² estimate: $((-2*dx_dT_raw/T) + (dXZ/(kB*T*T)))")

# Estimate for ∂²<X>/∂T² = (-2/T)*∂⟨X⟩/∂T + (1/(kB T²)) (∂⟨X*V₂⟩/∂T - ∂⟨X⟩⟨V₂⟩/∂T)
# ∂XV₂_∂T, red =  derivative_with_opt_c(X .* V2, V2, T, n_atoms, K = 2)
# ∂AB_∂T = (μX * ∂μ₂_∂T) + (μ₂ * ∂X_∂T)
# ∂²X_∂T² = (-2*∂X_∂T/T) + ((1/(kB*T*T)) * (∂XV₂_∂T - ∂AB_∂T))

    # Bootstrap to get error estimates
# is = zeros(Int, length(X))
# μX_b = zeros(n_boot)
# ∂X_∂T_b = zeros(n_boot)
# ∂²X_∂T²_b = zeros(n_boot)

# # Storage
# X_b = similar(X); Y_b = similar(X); Z_b = similar(X)
# C1_b = similar(C1); C2_b = similar(C2)
# C3_b = similar(C3); C4_b = similar(C4)
# for i in 1:n_boot
#     sample!(1:length(X), is; replace = true)
#     X_b .= X[is]; V2_b = V2[is]
#     C1_b .= C1[is]; C2_b .= C2[is];
#     C3_b .= C3[is]; C4_b .= C4[is];
#     Y_b = X_b .* V2_b
#     Z_b = X_b .* C1_b

#     r1 = cv_estimate(X_b, C1_b, C2_b, C3_b, C4_b)
#     ∂X_∂T_b[i], _ =  derivative_with_opt_c(X_b, V2_b, T, n_atoms, K = 3)
#     # ∂XV₂_∂T_b, _  =  derivative_with_opt_c(Y_b, V2_b, T, n_atoms, K = 3)
#     dXZ_1, _ = derivative_with_opt_c(Z_b, V2_b, T, n_atoms, K = 3)
#     dXZ = dXZ_1 - (∂μ₂_∂T * μX_b[i])

#     μX_b[i] = r1.mean_cv
#     ∂AB_∂T_b = (μX_b[i] * ∂μ₂_∂T) + (μ₂ * ∂X_∂T_b[i])
#     # ∂²X_∂T²_b[i] = (-2*∂X_∂T_b[i]/T) + ((1/(kB*T*T)) * (∂XV₂_∂T_b - ∂AB_∂T_b))
#     ∂²X_∂T²_b[i] = (-2*∂X_∂T_b[i]/T) + (dXZ/(kB*T*T))
# end

# μX_SE = std(μX_b)
# ∂X_∂T_SE = std(∂X_∂T_b)
# ∂²X_∂T²_SE = std(∂²X_∂T²_b)

# println("⟨X⟩ = $(round(μX, digits=5)) ± $(round(μX_SE, digits=5))")
# println("∂⟨X⟩/∂T = $(round(∂X_∂T, digits=5)) ± $(round(∂X_∂T_SE, digits=5))")
# println("∂²⟨X⟩/∂T² = $(round(∂²X_∂T², digits=5)) ± $(round(∂²X_∂T²_SE, digits=5))")

# return μX, μX_SE, ∂X_∂T, ∂X_∂T_SE, ∂²X_∂T², ∂²X_∂T²_SE