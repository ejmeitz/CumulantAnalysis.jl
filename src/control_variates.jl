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
    @warn "No control variates provided, returning raw estimate." maxlog=1
    return mean(X), ControlVariateData(zeros(T,0), zeros(T,0), T(1.0))
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
    R = diagm(0 => fill(λ, p))
    A = (Z' * Z) .+ R    # α = (Z'Z) \ (Z' * X)
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

    return C1, C2, C3, μ₂, ∂μ₂_∂T
end

function my_cov(X, V2, μ2, μX = mean(X))
    return sum((X .- μX) .* (V2 .- μ2)) / length(X)
end

# Fit from scratch
function get_cv_estimates(X, V2, T, n_atoms, use_cvs::Bool)
    
    cvs...,  μ₂, ∂μ₂_∂T = build_cvs(V2, T, n_atoms)
    Z = cvs[1]

    cvs = use_cvs ? cvs : ()

    # Estimate for <X>
    μX, cvd1 = cv_estimate(X, cvs...)

    # Estimate for ∂<X>/∂T ∝ cov(X, V2) = <XZ>
    Y1 = (X .- μX) .* Z
    cov_XZ, cvd2 = cv_estimate(Y1, cvs...)
    ∂X_∂T = cov_XZ / (kB * T^2)

    # Estimate for ∂²<X>/∂T² ∝ cov(XZ, V2) = <XZ^2>
    Y2 = (Y1 .- cov_XZ) .* Z
    cov_XZZ, cvd3 = cv_estimate(Y2, cvs...)
    dXZ_1 = cov_XZZ / (kB * T^2)
    dXZ = dXZ_1 - (∂μ₂_∂T * μX)
    ∂²X_∂T² = (-2*∂X_∂T/T) + (dXZ/(kB*T*T))

    return μX, cvd1, ∂X_∂T, cvd2, ∂²X_∂T², cvd3

end

# Re-use alphas from before, 
function get_cv_estimates(X, V2, T, n_atoms, use_cvs::Bool, cvds...)

    cvs..., μ₂, ∂μ₂_∂T = build_cvs(V2, T, n_atoms)
    Z = cvs[1]

    cvs = use_cvs ? cvs : ()

    # Estimate for <X>
    μX = apply_cv(X, cvds[1], cvs...)

    # Estimate for ∂<X>/∂T
    Y1 = (X .- μX) .* Z
    cov_XZ = apply_cv(Y1, cvds[2], cvs...) 
    ∂X_∂T = cov_XZ / (kB * T^2)

    # Estimate for ∂²<X>/∂T²
    Y2 = (Y1 .- cov_XZ) .* Z
    dXZ_1 = apply_cv(Y2, cvds[3], cvs...) / (kB * T^2)
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

