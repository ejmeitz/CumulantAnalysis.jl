
# Custom implementation that takes advantage of the
# structure of the covariance to improve accuracy
function estimate_V0(ce::CumulantEstimator{O}, V, V₂, V₃, V₄, n_atoms, T)


    X = V₀_rv(ce, V, V₂, V₃, V₄)

    μ₂ = 0.5*(3*n_atoms - 3)*kB*T # ⟨V2⟩
    Z = V₂ .- μ₂
    Y = X .* Z
    
    # use zero mean control variates for now
    return cv_estimate(Y, Z, V₃, V₃ .* V₂; K=10, verbose=true)

end

