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

    # cn = cond(A)
    # if cn > 1e6
    #     @warn "Possibly ill-condition control variates, cond number: $(cn)"
    # end

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



"""
Cross-fitted (contiguous folds, constant ridge) control-variates estimator for E[X].

- Splits data into K contiguous blocks.
- On each fold: center & standardize CVs using *train* stats, fit ridge OLS,
  apply to test block using the same train means/scales, collect residuals.
- Returns (mean_cv, ControlVariateData(α_bar, μ_bar, vr)).

Notes:
- α_bar, μ_bar are fold-size–weighted diagnostics; reapplying them won’t
  exactly reproduce the cross-fitted mean (normal for cross-fitting).
"""
function cv_estimate_crossfit(X::AbstractVector{T}, cvs::AbstractVector{T}...;
                              K::Int=4,
                              ridge::Real=1e-6,
                              tol::Real=1e-4) where {T}

    n = length(X)

    # Stack CVs (n × p)
    W = hcat(cvs...)
    nW, p = size(W)
    @assert nW == n "X and W must have same number of rows."
    @assert n > p  "Need more samples (n=$(n)) than control variates (p=$(p))."

    mean_raw = mean(X)
    var_raw  = var(X; corrected=true)

    # ----- Contiguous K folds -----
    # boundaries: 0, ⌊n/K⌋, ⌊2n/K⌋, ..., n
    starts = round.(Int, range(0, n; length=K+1))
    folds  = [ (starts[k]+1):starts[k+1] for k in 1:K ]

    resid   = similar(X, Float64)
    α_accum = zeros(Float64, p)   # diagnostics (fold-size–weighted)
    μ_accum = zeros(Float64, p)
    n_accum = 0

    λ = T(ridge)
    R = diagm(0 => fill(λ, p))

    for k in 1:K
        test  = folds[k]
        train = vcat(folds[setdiff(1:K,(k,))]...)
        nt, ntr = length(test), length(train)
        @assert ntr > p "Train block (size=$(ntr)) must exceed #CVs (p=$(p))."

        # ---- Train-fold centering & scaling for conditioning ----
        μ_tr   = @views vec(mean(W[train, :]; dims=1))            # p
        Ztr_c  = @views W[train, :] .- μ_tr'                      # ntr×p
        σ_tr   = vec(std(Ztr_c; dims=1, corrected=true))         # p
        @inbounds for j in 1:p
            if σ_tr[j] == 0
                σ_tr[j] = one(Float64)  # guard constant column
            end
        end
        Ztr = Ztr_c ./ σ_tr'

        # Ridge solve (SPD ⇒ use Cholesky)
        A = (Ztr' * Ztr ).+ R
        α_std = @views cholesky(A) \ (Ztr' * X[train])         # p
        α = α_std ./ σ_tr   # un-normalize

        # ---- Apply on test block using *train* μ ----
        Zte = (@views W[test, :]) .- μ_tr'
        resid[test] = @views X[test] .- (Zte * α)

        # diagnostics
        α_accum .+= nt .* α
        μ_accum .+= nt .* μ_tr
        n_accum  += nt
    end

    mean_cv = mean(resid)
    var_cv  = var(resid; corrected=true)
    vr      = var_raw / var_cv

    # sanity: mean drift relative to raw
    rel_err = abs(mean_cv - mean_raw) / max(abs(mean_raw), eps(T))
    if rel_err > tol
        @warn "Cross-fitted CV mean differs from raw mean by $(round(rel_err*100, digits=2))%."
    end

    α_bar = α_accum / n_accum
    μ_bar = μ_accum / n_accum

    return T(mean_cv), ControlVariateData(T.(α_bar), T.(μ_bar), T(vr))
end


function build_zero_mean_cvs(V2, V3, T, n_atoms)
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
    C4 = V3 .* C1

    return C1, C2, C3, C4, V3, μ₂, ∂μ₂_∂T
end


function my_cov(X, V2, μ2, μX = mean(X))
    return sum((X .- μX) .* (V2 .- μ2)) / length(X)
end

# Fit from scratch
function get_cv_estimates(X, V2, V3, T, n_atoms, use_cvs::Bool)
    
    cvs...,  μ₂, ∂μ₂_∂T = build_zero_mean_cvs(V2, V3, T, n_atoms)
    Z = cvs[1]
    cvs_nz = (X, X .^ 2, X .* Z)

    cvs = use_cvs ? cvs : ()
    cvs_nz = use_cvs ? cvs_nz : ()

    # Estimate for <X>
    μX, cvd1 = cv_estimate_crossfit(X, cvs...)

    # Estimate for ∂<X>/∂T ∝ cov(X, V2) = <XZ>
    Y1 = (X .- μX) .* Z
    cov_XZ, cvd2 = cv_estimate_crossfit(Y1, cvs..., cvs_nz[1], cvs_nz[2])
    ∂X_∂T = cov_XZ / (kB * T^2)

    # Estimate for ∂²<X>/∂T² ∝ cov(XZ, V2) = <XZ^2>
    tmp = X .* V2
    Y2 = (tmp .- mean(tmp)) .* Z #! CAN I REUSE DATA FROM BEFORE TO IMPROVE ACCURACY?
    cov_XZZ, cvd3 = cv_estimate_crossfit(Y2, cvs..., cvs_nz...)
    dXZ_1 = cov_XZZ / (kB * T^2)
    dXZ = dXZ_1 - (∂μ₂_∂T * μX)
    ∂²X_∂T² = (-2*∂X_∂T/T) + (dXZ/(kB*T*T))

    return μX, cvd1, ∂X_∂T, cvd2, ∂²X_∂T², cvd3

end

# Re-use alphas from before, 
function get_cv_estimates(X, V2, V3, T, n_atoms, use_cvs::Bool, cvds...)

    cvs...,  μ₂, ∂μ₂_∂T = build_zero_mean_cvs(V2, V3, T, n_atoms)
    Z = cvs[1]
    cvs_nz = (X, X .^ 2, X .* Z)

    cvs = use_cvs ? cvs : ()
    cvs_nz = use_cvs ? cvs_nz : ()

    # Estimate for <X>
    μX = apply_cv(X, cvds[1], cvs...) # do not pass cvs that depend on X

    # Estimate for ∂<X>/∂T
    Y1 = (X .- μX) .* Z
    cov_XZ = apply_cv(Y1, cvds[2], cvs..., cvs_nz[1], cvs_nz[2]) 
    ∂X_∂T = cov_XZ / (kB * T^2)

    # Estimate for ∂²<X>/∂T²
    Y2 = (Y1 .- cov_XZ) .* Z
    dXZ_1 = apply_cv(Y2, cvds[3], cvs..., cvs_nz...) / (kB * T^2)
    dXZ = dXZ_1 - (∂μ₂_∂T * μX)
    ∂²X_∂T² = (-2*∂X_∂T/T) + (dXZ/(kB*T*T))

    return μX, ∂X_∂T, ∂²X_∂T²

end



# """
#     derivative_with_opt_c(X, V2, T, n_atoms; K::Int=2)


# Single control variate estimator for the derivative d⟨X⟩/dT.

# Estimate:
#   cov_hat = ⟨(X - c*) * (V2 - μ2)⟩_0
#   deriv_hat = cov_hat / (kB * T^2)

# where μ2 = (f/2) * kB * T and Var(V2) = (f/2) * (kB*T)^2 with f = 3*n_atoms - 3.

# Arguments
# - `X::AbstractVector{<:Real}`: your `X` samples (e.g., V - V2 or any statistic).
# - `V2::AbstractVector{<:Real}`: harmonic potential energy samples.
# - `T::Real`: temperature (K).
# - `n_atoms::Integer`: number of atoms (used for f = 3*n_atoms).
# - `K`: folds for cross-fitting (default 2).

# Returns
# - `(; cov_hat, cov_se, deriv_hat, deriv_se, c_folds, mu2, varV2_analytic)`

# """
# function derivative_with_opt_c(X::AbstractVector{<:Real},
#                                V2::AbstractVector{<:Real},
#                                T::Real,
#                                n_atoms::Integer; K::Int=2)

#     @assert length(X) == length(V2) "X and V2 must have the same length."
#     n = length(X)
#     n > 0 || error("Empty inputs.")

#     # Degrees of freedom
#     f = 3 * n_atoms - 3

#     # True harmonic mean and variance of V2
#     mu2 = 0.5 * f * kB * T
#     varV2_analytic = 0.5 * f * (kB * T)^2  # E[Z^2] with Z = V2 - mu2

#     Z = V2 .- mu2

#     # Build K folds
#     idx = collect(1:n)
#     Random.shuffle!(idx)
#     folds = [idx[round(Int, floor((k-1)*n/K))+1 : round(Int, floor(k*n/K))] for k in 1:K]
#     # guard against rounding issues
#     folds[end] = idx[sum(length.(folds[1:end-1]))+1:end]

#     # Cross-fit c on complement of each fold; evaluate on the fold
#     y = similar(float.(X))  # per-sample integrand values (with each fold's c)
#     c_folds = Float64[]
#     for fold in folds
#         comp = setdiff(idx, fold)
#         # Numerator: E[X Z^2] estimated on complement (use mean for scale invariance)
#         num = mean(X[comp] .* (Z[comp]).^2)
#         # Denominator: analytic E[Z^2] = Var(V2)
#         c_hat = num / varV2_analytic
#         push!(c_folds, c_hat)

#         # Evaluate integrand on held-out fold
#         @inbounds @views y[fold] .= (X[fold] .- c_hat) .* Z[fold]
#     end

#     cov_hat = mean(y)
#     deriv_hat = cov_hat / (kB * T^2)

#     raw_varY = var(X .* Z; corrected=true)
#     red = raw_varY / var(y; corrected=true)

#     return deriv_hat, red
# end

