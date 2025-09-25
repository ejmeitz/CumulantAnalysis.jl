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
    var_raw  = var(X; corrected=true)

    if !crossfit
        # Single-fit (small O(1/n) bias): center by full-sample means
        μ_estimate = mean(W; dims=1)               # p
        Z  = W .- μ_estimate                       # n×p, zero-mean by sample
        α  = (Z' * Z) \ (Z' * X)                   # p
        resid = X .- Z * α
        mean_cv = mean(resid)                   # ≈ unbiased as n→∞
        var_resid = var(resid; corrected=true)
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
        mean_fold = @views mean(X[F] .- (W[F, :] * α)) + dot(μ_estimate, α)
        push!(fold_means, mean_fold)

        # Residuals on test (for variance reporting)
        Z_te = @views W[F, :] .- ones(length(F)) * μ_estimate'
        resid_all[F] = X[F] .- Z_te * α
    end

    mean_cv   = mean(fold_means)                  # unbiased
    var_resid = var(resid_all; corrected=true)
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

