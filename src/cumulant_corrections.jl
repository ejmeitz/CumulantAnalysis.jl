# Various derivatives of ⟨O⟩
∂A_∂T(A, V, T) = cov(A, V) / (kB * T * T)
∂AB_∂T(A, B, V, T, dA = ∂A_∂T(A, V, T), dB =  ∂A_∂T(B, V, T)) = (mean(A) * dB) + (mean(B) * dA)
∂²A_∂T²(A, V, T, dA = ∂A_∂T(A, V, T)) = (-2*dA/T) + ((1/(kB*T*T)) * (∂A_∂T(A.*V, V, T) - ∂AB_∂T(A, V, V, T, dA)))

# probably biased
function central_moment(X, n::Int)
    return mean( (X .- mean(X)) .^ n )
end

skew(X) = central_moment(X, 3)

## CONSTANT CORRECTION (Order Zero) ##

function CumulantData(V, V₂, V₃, V₄, V_ref, T, ::Val{0}, ce::AnalyticalEstimator, V4_analytical)

    X = V₀_rv(ce, V, V₂, V₃, V₄)

    # Uses V3/V4 as control variates for estimating <V - V2 - V3 - V4>_0
    # Allocates each iteration...
    # Xcv, α = cv_analytical_estimator(X, V₃, V₄, V4_analytical)

    V₀ = mean(X)
    ∂V₀ = ∂A_∂T(X, V_ref, T)
    ∂²V₀ = ∂²A_∂T²(X, V_ref, T, ∂V₀)

    return CumulantData{0, typeof(V₀), typeof(∂V₀), typeof(∂²V₀)}(V₀, ∂V₀, ∂²V₀)

end

function constant_corrections(c0::CumulantData{0}, T)

    F_corr = c0.κ
    S_corr = -c0.∂κ_∂T
    U_corr = c0.κ - T*c0.∂κ_∂T
    Cv_corr = -T * c0.∂²κ_∂T²

    return F_corr, S_corr, U_corr, Cv_corr
    
end

function first_order_corrections(c1::CumulantData{1}, T)
   
    F_correction = c1.κ
    S_correction = -c1.∂κ_∂T
    U_correction = c1.κ - T*c1.∂κ_∂T
    Cv_correction = -T*c1.∂²κ_∂T²

    return F_correction, S_correction, U_correction, Cv_correction
end

function second_order_corrections(c2::CumulantData{2}, T, stochastic::Bool)

    pref = stochastic ? -1.0 : 1.0
    β = 1 / (kB*T)

    F_correction = pref * c2.κ / (2*kB*T)
    S_correction =  pref * (c2.κ - T*c2.∂κ_∂T) / (2*kB*T*T)
    U_correction =  pref * β * (c2.κ - 0.5*T*c2.∂κ_∂T)
    # Cv_correction =  pref * ((-U_correction/T) + β*(0.5*c2.∂κ_∂T - 0.5*T*c2.∂²κ_∂T²))
    Cv_correction = pref * β * ((-c2.κ/T) + c2.∂κ_∂T - ((T/2)*(c2.∂²κ_∂T²)))

    return F_correction, S_correction, U_correction, Cv_correction
end


## CONSTANT CORRECTION ##
# function CumulantData(V, V₂, V₃, V₄, T, n_atoms, ::Val{0}, ce::SamplingCumulantEstimator)

#     X = V₀_rv(ce, V, V₂, V₃, V₄)

#     V₀ = mean(X)
#     ∂V₀ = ∂A_∂T(X, V₂, T)
#     ∂²V₀ = ∂²A_∂T²(X, V₂, T, ∂V₀)

#     # # This estimator uses a user provided V0
#     # if ce isa MixedEstimator
#     #     V₀ = get_V₀(ce, V, V₂, V₃, V₄)
#     # end

#     return CumulantData{0, typeof(V₀), typeof(∂V₀), typeof(∂²V₀)}(V₀, ∂V₀, ∂²V₀)

# end

# ## FIRST CUMULANTS ##
# function CumulantData(V, V₂, V₃, V₄, T, n_atoms, ::Val{1}, ce::SamplingCumulantEstimator)

#     X = X1(ce, V, V₂, V₃, V₄)

#     κ₁ = mean(X)
#     ∂κ₁_∂T = ∂A_∂T(X, V₂, T)
#     ∂²κ₁_∂T² = ∂²A_∂T²(X, V₂, T)

#     return CumulantData{1, typeof(κ₁), typeof(∂κ₁_∂T), typeof(∂²κ₁_∂T²)}(κ₁, ∂κ₁_∂T, ∂²κ₁_∂T²)
# end


# ## SECOND CUMULANTS ##
# function CumulantData(V, V₂, V₃, V₄, T, n_atoms, c1::CumulantData{1}, ::Val{2}, ce::SamplingCumulantEstimator)

#     X = X2(ce, V, V₂, V₃, V₄)
#     X² = X .^ 2

#     # μX², cvd1, ∂X²_∂T, cvd2, ∂²X²_∂T², cvd3 = get_cv_estimates(X², V₂, V₃, T, n_atoms, use_cvs)
#     # κ₂ =  μX² - c1.κ^2

#     κ₂ = var(X; corrected = true)
#     ∂X²_∂T = ∂A_∂T(X², V₂, T)
#     ∂²X²_∂T² = ∂²A_∂T²(X², V₂, T)

#     ∂κ₂_∂T = ∂X²_∂T - (2*c1.κ*c1.∂κ_∂T)
#     ∂²κ₂_∂T² = ∂²X²_∂T² - 2*(((c1.∂κ_∂T)^2) + (c1.κ*c1.∂²κ_∂T²))

#     return CumulantData{2, typeof(κ₂), typeof(∂κ₂_∂T), typeof(∂²κ₂_∂T²)}(κ₂, ∂κ₂_∂T, ∂²κ₂_∂T²)
# end


## THRID CUMULANTS ##

# function CumulantData(V, V₂, V₃, V₄, T, c1::CumulantData{1}, ::Val{3}, ce::SamplingCumulantEstimator)

#     X = X3(ce, V, V₂, V₃, V₄)
#     X² = X .^ 2
#     X³ = X .^ 3

#     t1 = Threads.@spawn skew(X)
#     t2 = Threads.@spawn ∂A_∂T(X³, V₂, T)
#     t3 = Threads.@spawn ∂²A_∂T²(X³, V₂, T)
#     t4 = Threads.@spawn mean(X²)

#     # These are both re calculated from second order
#     t5 = Threads.@spawn ∂A_∂T(X², V₂, T)
#     t6 = Threads.@spawn ∂²A_∂T²(X², V₂, T)

#     κ₃ = fetch(t1)
#     ∂X³_∂T = fetch(t2)
#     ∂²X³_∂T² = fetch(t3)
#     μ_X² = fetch(t4)
#     ∂X²_∂T = fetch(t5)
#     ∂²X²_∂T² = fetch(t6)

#     ∂κ₃_∂T = ∂X³_∂T - 3*c1.κ*∂X²_∂T + 3*μ_X²*c1.∂κ_∂T
#     ∂²κ₃_∂T² = ∂²X³_∂T² - 3*(c1.∂κ_∂T*∂X²_∂T + c1.κ*∂²X²_∂T²) + 3*(∂X²_∂T*c1.∂κ_∂T + μ_X²*c1.∂²κ_∂T²)

#     return CumulantData{3, typeof(κ₃), typeof(∂κ₃_∂T), typeof(∂²κ₃_∂T²)}(
#                             κ₃, ∂κ₃_∂T, ∂²κ₃_∂T²)
# end

# function third_order_corrections(c3::CumulantData{3}, T)

#     β = 1 / (kB*T)
#     β² = β^2; β³ = β^3

#     # no clue which of these are correct
#     F_correction = c3.κ * β² / 6
#     S_correction = 0.0 #(c3.κ*kB*β³/3) - (β²*c3.∂κ_∂T/6)
#     U_correction = 0.0 #T*((0.5*kB*β³*c3.κ) - (β²*c3.∂κ_∂T/6))
#     Cv_correction = 0.0 #(U_correction/T) - (3*β²*c3.κ/2) + (5*β²*c3.∂κ_∂T/6) + (T*β²*c3.∂²κ_∂T²/6)

#     return F_correction, S_correction, U_correction, Cv_correction
# end
