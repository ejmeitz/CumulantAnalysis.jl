
abstract type Limit end
struct Quantum <: Limit end
struct Classical <: Limit end

struct CumulantData{A,B,C,D,E,F}
    κ₁::A
    ∂κ₁_∂T::B
    ∂²κ₁_∂T²::C
    κ₂::D
    ∂κ₂_∂T::E
    ∂²κ₂_∂T²::F
end

#* definitely check the second derivatives!!
function CumulantData(V, ΔV, kB, T)

    ΔV² = ΔV .^ 2

    κ₁ = mean(ΔV)
    ∂κ₁_∂T = cov(V, ΔV)/(kB*T*T)
    ∂²κ₁_∂T² = (2.0/(kB*T*T*T)) * (cov(ΔV, V.^2) - mean(V)*cov(ΔV, V)) - ∂κ₁_∂T

    κ₂ = var(ΔV)
    ∂κ₂_∂T = (1/(kB*T*T)) * (cov(ΔV², V) - 2*U0*cov(ΔV, V))
    ∂²κ₂_∂T² = ((-2/kB*T*T*T)*∂κ₂_∂T) + ((1/(kB*kB*T*T*T*T)) * (cov((ΔV² .* V, V) - mean(V)*cov(ΔV², V) - mean(ΔV²)*var(V))))

    return CumulantData(κ₁, ∂κ₁_∂T, ∂²κ₁_∂T², κ₂, ∂κ₂_∂T, ∂²κ₂_∂T²)
end

function V₂()

end


#* ARE ALL OF THESE MISSING A FACTOR OF 1/Nq ???

function U₀(ω, ħ, kB, T, limit::Quantum)
    mode_internal_energies = @. (ħ*ω) * ((1 / (exp(ħ*ω/(kB*T)) - 1)) + 0.5)
    return sum(mode_internal_energies)
end

function U₀(ω, ħ, kB, T, limit::Classical)
    return length(ω)*kB*T
end

function F₀(ω, ħ, kB, T, limit::Quantum)
    kBT = kB * T
    mode_free_energies = @. (0.5*ħ*ω) + kBT * log(1 - exp(-ħ*ω/kBT))
    return sum(mode_free_energies)
end

function F₀(ω, kB, T, limit::Classical)

end

function S₀(ω, ħ, kB, T, limit::Quantum)
    ħω_kBT = ħ .* ω ./ (kB * T)
    mode_entropies = @. ħω_kBT / (exp(ħω_kBT) - 1) - log(1 - exp(-ħω_kBT))
    return kB * sum(mode_entropies)
end

function S₀(ω, kB, T, limit::Classical)

end

function Cᵥ₀(ω, kB, T, limit::Quantum)
    ħω_2kBT = ħ * ω / (2 * kB * T)
    mode_heat_capacities = @. ((ħω_2kBT)^2) * (csch(ħω_2kBT)^2)
    return kB * sum(mode_heat_capacities)
end

function Cᵥ₀(ω, kB, T, limit::Classical)
    return length(ω)*kB
end


function first_order(cd::CumulantData, T)
   
    F_correction = cd.κ₁
    S_correction = -cd.∂κ₁_∂T
    U_correction = cd.κ₁ - T*cd.∂κ₁_∂T
    Cv_correction = -T*cd.∂²κ₁_∂T²

    return F_correction, S_correction, U_correction, Cv_correction
end

function second_order(cd::CumulantData, kB, T, stochastic::Bool)

    pref = stochastic ? -1.0 : 1.0

    F_correction = pref * cd.κ₂ / (2*kB*T)
    S_correction =  pref * (cd.κ₂ - T*cd.∂κ₂_∂T) / (2*kB*T*T)
    U_correction =  pref * (cd.κ₂ - 0.5*T*cd.∂κ₂_∂T) / (kB*T)
    #* Need to double check this
    Cv_correction =  (pref / kB) * ((cd.∂κ₂_∂T /T) - (cd.κ₂/(T*T)) - (0.5 * cd.∂²κ₂_∂T²))

    return F_correction, S_correction, U_correction, Cv_correction
end