export cv_estimate

#! NOT SURE TEHSE ARE RIGHT

# ⟨X⟩ - αᵀ (⟨W⟩ - μ_W)
function apply_cv(X::AbstractVector{T},
                  cvd::ControlVariateData{T},
                  zero_mean_cvs::AbstractVector{T}...
                ) where T

    W = hcat(zero_mean_cvs...)
    Z = W .- cvd.μW_estimate'
    return mean(X - Z * cvd.αs)
end

function apply_cv(X::AbstractVector{T},
                  cvd::ControlVariateData{T},
                  W::AbstractMatrix{T}
                ) where T
    Z = W .- cvd.μW_estimate'
    return mean(X - Z * cvd.αs) 
end

# Path when not using control variates
function apply_cv(X::AbstractVector{T},
                  cvd::ControlVariateData{T}
                ) where T
    return mean(X)
end


function cv_estimate(X::AbstractVector{T}) where T
    mean_raw = mean(X)
    @warn "No control variates provided, returning raw estimate." maxlog=1
    return mean_raw, ControlVariateData(zeros(T,0), zeros(T,0), T(1.0))
end

"""
    cv_estimate(X, zero_mean_cvs...)

Control-variates estimator for E[X] using zero_mean_cvs as controls.
ASSUMES means of control variates are 0. Does 2-fold cross-fitting
to remove a little bias.

Inputs:
  X :: AbstractVector       (length n)
  zero_mean_cvs :: AbstractVector...  


Returns NamedTuple:
  (mean_raw, mean_cv, alpha, var_raw, var_residual, variance_reduction)
"""
function cv_estimate(X::AbstractVector{T}, zero_mean_cvs::AbstractVector{T}...;
                     tol = 1e-4, ridge::Real = 1e-6) where T
    
    mean_raw = mean(X)
    var_raw  = var(X; corrected=true)

    W = hcat(zero_mean_cvs...) # n x p 
    n, p = size(W)
    @assert n > p "Need more samples (n=$(n)) than control variates (p=$(p))."
    @assert length(X) == n "X and W must have same number of rows."

    # C.V. should already have zero mean, but this ensures
    # numerical stability and corrects any small deviations
    μ_estimate = mean(W; dims=1)               # 1×p
    Z = W .- μ_estimate                        # n×p

    λ = T(ridge)
    A = Z'Z .+ λ*I    # α = (Z'Z) \ (Z' * X)
    α = A \ (Z' * X)

    cn = cond(A)
    if cn > 1e6
        @warn "Possibly ill-condition control variates, cond number: $(cn)"
    end

    # Estimator for E[X]
    mean_cv = mean(X .- Z * α)
    resid = X .- Z * α 
    var_cv = var(resid; corrected=true)
    
    # Variance reduction ratio
    vr = var_raw / var_cv
    
    # Sanity check (mean should be unchanged when CVs have zero mean)
    rel_err = abs(mean_cv - mean_raw) / max(abs(mean_raw), eps(T))
    if rel_err > tol
        @warn "Control variate mean differs from raw mean by $(round(rel_err*100, digits=2))%."
    end

    return mean_cv, ControlVariateData(α, vec(μ_estimate), vr)
end

function build_cvs(V2, T, n_atoms)
    f = 3*n_atoms - 3
    μ₂ = 0.5*f*kB*T  # ⟨V2⟩
    ∂μ₂_∂T = 0.5*f*kB
    σ₂² = 0.5*f*(kB*T)^2 # var(V2)

    # Build some control variates
    # These all have zero mean by construction
    C1 = V2 .- μ₂
    C1_sq = C1 .^ 2
    C2 = C1_sq .- σ₂²
    C3 = C1 .^ 3 .- (3*σ₂²*C1)

    return C1, C2, C3, ∂μ₂_∂T
end

# Fit from scratch
function get_cv_estimates(X, V2, T, n_atoms, use_cvs::Bool)
    
    cvs, ∂μ₂_∂T = build_cvs(V2, T, n_atoms)

    cvs = use_cvs ? cvs : ()

    # Estimate for <X>
    μX, cvd1 = cv_estimate(X, cvs...)

    # Estimate for ∂<X>/∂T
    Y = X .* C1
    cov_XZ, cvd2 = cv_estimate(Y, cvs...)
    ∂X_∂T = cov_XZ / (kB * T^2)

    # Estimate for ∂²<X>/∂T²
    cov_XZZ, cvd3 = cv_estimate(Y .* C1, cvs...)
    dXZ_1 = cov_XZZ / (kB * T^2)
    dXZ = dXZ_1 - (∂μ₂_∂T * μX)
    ∂²X_∂T² = (-2*∂X_∂T/T) + (dXZ/(kB*T*T))

    return μX, cvd1, ∂X_∂T, cvd2, ∂²X_∂T², cvd3

end

# Re-use alphas from before, 
function get_cv_estimates(X, V2, T, n_atoms, use_cvs::Bool, cvds...)

    cvs, ∂μ₂_∂T = build_cvs(V2, T, n_atoms)

    cvs = use_cvs ? cvs : ()

    # Estimate for <X>
    μX = apply_cv(X, cvds[1], cvs...)

    # Estimate for ∂<X>/∂T
    Y = X .* C1
    ∂X_∂T = apply_cv(Y, cvds[2], cvs...) / (kB * T^2)

    # Estimate for ∂²<X>/∂T²
    dXZ_1 = apply_cv(Y .* C1, cvds[3], cvs...) / (kB * T^2)
    dXZ = dXZ_1 - (∂μ₂_∂T * μX)
    ∂²X_∂T² = (-2*∂X_∂T/T) + (dXZ/(kB*T*T))

    return μX, ∂X_∂T, ∂²X_∂T²

end

# Assumes harmonic reference
# function get_cv_estimates(X, V2, T, n_atoms)

#     # Build some control variates
#     # These all have zero mean by construction
#     f = 3*n_atoms - 3
#     μ₂ = 0.5*f*kB*T  # ⟨V2⟩
#     ∂μ₂_∂T = 0.5*f*kB
#     σ₂² = 0.5*f*(kB*T)^2 # var(V2)

#     C1 = V2 .- μ₂
#     C1_sq = C1 .^ 2
#     C2 = C1_sq .- σ₂²
#     C3 = C1 .^ 3 .- (3*σ₂²*C1)

#     # Estimate for <X>
#     res1 = cv_estimate(X, C1, C2, C3)
#     μX = res1.mean_cv

#     # Estimate for ∂<X>/∂T
#     ∂X_∂T, red = derivative_with_opt_c(X, V2, T, n_atoms, K = 2)

#     # Estimate for ∂²<X>/∂T²
#     dXZ_1, red = derivative_with_opt_c(X .* C1, V2, T, n_atoms, K = 3)
#     dXZ = dXZ_1 - (∂μ₂_∂T * μX)
#     ∂²X_∂T² = (-2*∂X_∂T/T) + (dXZ/(kB*T*T))

#     return μX, ∂X_∂T, ∂²X_∂T²

# end


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