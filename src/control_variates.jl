export cv_estimate

"""
    cv_estimate(X, W; crossfit::Bool=true, K::Int=5, verbose::Bool=true)

Control-variates estimator for E[X] using columns of W as controls.
Control means are estimated from data. With `crossfit=true` (default),
uses K-fold cross-fitting for an unbiased estimator.

Inputs:
  X :: AbstractVector       (length n)
  W :: AbstractMatrix       (n × p), each column a control variate

Keywords:
  crossfit :: Bool = true   # unbiased via K-fold cross-fit
  K        :: Int  = 5
  verbose  :: Bool = true   # print diagnostics

Returns NamedTuple:
  (mean_raw, mean_cv, alpha, var_raw, var_residual, variance_reduction)
"""
function cv_estimate(X::AbstractVector, cvs::AbstractVector...;
                     crossfit::Bool=true, K::Int=5, verbose::Bool=true)

    W = hcat(cvs...)

    n, p = size(W)
    @assert length(X) == n "X and W must have same number of rows."

    mean_raw = mean(X)
    var_raw  = var(X)

    if !crossfit
        # Single-fit (small O(1/n) bias): center by full-sample means
        μ_estimate = mean(W; dims=1)               # p
        Z  = W .- μ_estimate                       # n×p, zero-mean by sample
        α  = (Z' * Z) \ (Z' * X)                   # p
        resid = X .- Z * α
        mean_cv = mean(resid)                   # ≈ unbiased as n→∞
        var_resid = var(resid)
        if verbose
            println("Control Variates (single fit; small O(1/n) bias)")
            println("  raw mean             = $(mean_raw)")
            println("  CV mean              = $(mean_cv)")
            println("  variance(raw)        = $(var_raw)")
            println("  variance(residual)   = $(var_resid)")
            println("  variance reduction   = $(var_raw/var_resid)×")
        end
        return (mean_raw=mean_raw, mean_cv=mean_cv, alpha=α,
                var_raw=var_raw, var_residual=var_resid,
                variance_reduction=var_raw/var_resid)
    end

    # K-fold cross-fit (unbiased)
    idx = collect(1:n)
    shuffle!(idx)
    cuts = round.(Int, range(0, n, length=K+1))
    folds = [idx[cuts[k]+1 : cuts[k+1]] for k in 1:K]

    fold_means = Float64[]
    resid_all  = similar(X)   # to report variance of residual integrand
    alphas     = Vector{Float64}[]

    for F in folds
        Tset = setdiff(1:n, F)
        # Means from training split
        μ_estimate = vec(mean(W[Tset, :]; dims=1))      # p
        Z_tr  = W[Tset, :] .- ones(length(Tset)) * μ_estimate'
        α     = (Z_tr' * Z_tr) \ (Z_tr' * X[Tset])    # p
        push!(alphas, α)

        # Unbiased fold estimate: mean(X_te - W_te*α) + μ_estimate' α
        mean_fold = @views mean(X[F] .- (W[F, :] * α)) + sum(μ_estimate .* α)
        push!(fold_means, mean_fold)

        # Residuals on test (for variance reporting)
        Z_te = @views W[F, :] .- ones(length(F)) * μ_estimate'
        resid_all[F] = X[F] .- Z_te * α
    end

    mean_cv   = mean(fold_means) # unbiased
    var_resid = var(resid_all)
    if verbose
        println("Control Variates (K-fold cross-fit, unbiased)")
        println("  raw mean             = $(mean_raw)")
        println("  CV mean              = $(mean_cv)")
        println("  variance(raw)        = $(var_raw)")
        println("  variance(residual)   = $(var_resid)")
        println("  variance reduction   = $(var_raw/var_resid)×")
    end
    # Return average α just for diagnostics
    α_avg = reduce(+, alphas) ./ length(alphas)

    return (mean_raw=mean_raw, mean_cv=mean_cv, alpha=α_avg,
            var_raw=var_raw, var_residual=var_resid,
            variance_reduction=var_raw/var_resid)
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

    # Degrees of freedom (adjust if you have constraints/rigid modes to remove)
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

    # Final estimates
    cov_hat = mean(y)
    cov_se  = sqrt(var(y; corrected=true) / n)  # standard error for the mean
    deriv_hat = cov_hat / (kB * T^2)
    deriv_se  = cov_se  / (kB * T^2)

    raw_varY = var(X .* Z; corrected=true)
    red = raw_varY / var(y; corrected=true)

    return deriv_hat, deriv_se, red
end