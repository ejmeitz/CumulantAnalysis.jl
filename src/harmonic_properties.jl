function V_harmonic(ifc2::Matrix{T}, u::AbstractVector{U})
    return 0.5 * ((transpose(u) * ifc2) * u)
end


#* ARE ALL OF THESE MISSING A FACTOR OF 1/Nq ???

function U_harmonic(ω, ħ, kB, T, limit::Quantum)
    mode_internal_energies = @. (ħ*ω) * ((1 / (exp(ħ*ω/(kB*T)) - 1)) + 0.5)
    return sum(mode_internal_energies)
end

function U_harmonic(ω, ħ, kB, T, limit::Classical)
    return length(ω)*kB*T
end

function F_harmonic(ω, ħ, kB, T, limit::Quantum)
    kBT = kB * T
    mode_free_energies = @. (0.5*ħ*ω) + kBT * log(1 - exp(-ħ*ω/kBT))
    return sum(mode_free_energies)
end

function F_harmonic(ω, kB, T, limit::Classical)

end

function S_harmonic(ω, ħ, kB, T, limit::Quantum)
    ħω_kBT = ħ .* ω ./ (kB * T)
    mode_entropies = @. ħω_kBT / (exp(ħω_kBT) - 1) - log(1 - exp(-ħω_kBT))
    return kB * sum(mode_entropies)
end

function S_harmonic(ω, kB, T, limit::Classical)

end

function Cᵥ_harmonic(ω, kB, T, limit::Quantum)
    ħω_2kBT = ħ * ω / (2 * kB * T)
    mode_heat_capacities = @. ((ħω_2kBT)^2) * (csch(ħω_2kBT)^2)
    return kB * sum(mode_heat_capacities)
end

function Cᵥ_harmonic(ω, kB, T, limit::Classical)
    return length(ω)*kB
end