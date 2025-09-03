# Various derivatives of ⟨O⟩
∂A_∂T(A, V, T) = cov(A, V; corrected = false) / (kB * T * T)
∂AB_∂T(A, B, V, T, dA = ∂A_∂T(A, V, T), dB =  ∂A_∂T(B, V, T)) = (mean(A) * dB) + (mean(B) * dA)
∂²A_∂T²(A, V, T, dA = ∂A_∂T(A, V, T)) = (-2*dA/T) + ((1/(kB*T*T)) * (∂A_∂T(A.*V, V, T) - ∂AB_∂T(A, V, V, T, dA)))

function central_moment(X, n::Int)
    return mean( (X .- mean(X)) .^ n )
end

skew(X) = central_moment(X, 3)

## FIRST CUMULANTS ##
function CumulantData(V, V₂, V₃, V₄, T, ::Val{1}, ce::CumulantEstimator)

    X = X1(ce, V, V₂, V₃, V₄)
    t1 = Threads.@spawn mean(X)
    t2 = Threads.@spawn ∂A_∂T(X, V₂, T)
    t3 = Threads.@spawn ∂²A_∂T²(X, V₂, T)

    κ₁ = fetch(t1)
    ∂κ₁_∂T = fetch(t2)
    ∂²κ₁_∂T² = fetch(t3)

    return CumulantData{1, typeof(κ₁), typeof(∂κ₁_∂T), typeof(∂²κ₁_∂T²)}(
                            κ₁, ∂κ₁_∂T, ∂²κ₁_∂T²)
end

## SECOND CUMULANTS ##

function CumulantData(V, V₂, V₃, V₄, T, c1::CumulantData{1}, ::Val{2}, ce::CumulantEstimator)

    X = X2(ce, V, V₂, V₃, V₄)
    X² = X .^ 2

    t1 = Threads.@spawn var(X; corrected = false)
    t2 = Threads.@spawn ∂A_∂T(X², V₂, T)
    t3 = Threads.@spawn ∂²A_∂T²(X², V₂, T)

    κ₂ = fetch(t1)
    ∂X²_∂T = fetch(t2)
    ∂²X²_∂T² = fetch(t3)

    ∂κ₂_∂T = ∂X²_∂T - (2*c1.κ*c1.∂κ_∂T)
    ∂²κ₂_∂T² = ∂²X²_∂T² - 2*(((c1.∂κ_∂T)^2) + (c1.κ*c1.∂²κ_∂T²))

    return CumulantData{2, typeof(κ₂), typeof(∂κ₂_∂T), typeof(∂²κ₂_∂T²)}(
                        κ₂, ∂κ₂_∂T, ∂²κ₂_∂T²)
end

## THRID CUMULANTS ##

function CumulantData(V, V₂, V₃, V₄, T, c1::CumulantData{1}, ::Val{3}, ce::CumulantEstimator)

    X = X3(ce, V, V₂, V₃, V₄)
    X² = X .^ 2
    X³ = X .^ 3

    t1 = Threads.@spawn skew(X)
    t2 = Threads.@spawn ∂A_∂T(X³, V_ref, T)
    t3 = Threads.@spawn ∂²A_∂T²(X³, V_ref, T)
    t4 = Threads.@spawn mean(X²)

    # These are both re calculated from second order
    t5 = Threads.@spawn ∂A_∂T(X², V_ref, T)
    t6 = Threads.@spawn ∂²A_∂T²(X², V_ref, T)

    κ₃ = fetch(t1)
    ∂X³_∂T = fetch(t2)
    ∂²X³_∂T² = fetch(t3)
    μ_X² = fetch(t4)
    ∂X²_∂T = fetch(t5)
    ∂²X²_∂T² = fetch(t6)

    ∂κ₃_∂T = ∂X³_∂T - 3*c1.κ*∂X²_∂T + 3*μ_X²*c1.∂κ_∂T
    ∂²κ₃_∂T² = ∂²X³_∂T² - 3*(c1.∂κ_∂T*∂X²_∂T + c1.κ*∂²X²_∂T²) + 3*(∂X²_∂T*c1.∂κ_∂T + μ_X²*c1.∂²κ_∂T²)

    return CumulantData{3, typeof(κ₃), typeof(∂κ₃_∂T), typeof(∂²κ₃_∂T²)}(
                            κ₃, ∂κ₃_∂T, ∂²κ₃_∂T²)
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

function third_order_corrections(c3::CumulantData{3}, T)

    β = 1 / (kB*T)
    β² = β^2; β³ = β^3

    # no clue which of these are correct
    F_correction = c3.κ * β² / 6
    S_correction = 0.0 #(c3.κ*kB*β³/3) - (β²*c3.∂κ_∂T/6)
    U_correction = 0.0 #T*((0.5*kB*β³*c3.κ) - (β²*c3.∂κ_∂T/6))
    Cv_correction = 0.0 #(U_correction/T) - (3*β²*c3.κ/2) + (5*β²*c3.∂κ_∂T/6) + (T*β²*c3.∂²κ_∂T²/6)

    return F_correction, S_correction, U_correction, Cv_correction
end
