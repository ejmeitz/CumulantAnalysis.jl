export Quantum, Classical

abstract type Limit end
struct Quantum <: Limit end
struct Classical <: Limit end

#* ACCOUNT FOR MISSING DOFS IN CLASSICAL LIMITS?

function V_harmonic(ifc2::AbstractMatrix, u::AbstractVector)
    return 0.5 * ((transpose(u) * ifc2) * u)
end

function sum_over_freqs!(freqs, f::Function, freq_tol = 1e-6)
    res = 0.0
    for freq in freqs
        if freq > freq_tol
            res += f(freq)
        end
    end
    return res
end

function U_harmonic(ω, ħ, kB, T, limit::Quantum)
    f = (freq) -> (ħ*freq) * ((1 / (exp(ħ*freq/(kB*T)) - 1)) + 0.5)
    return sum_over_freqs!(ω, f)
end

function U_harmonic(ω, ħ, kB, T, limit::Classical)
    return length(ω)*kB*T
end

function F_harmonic(ω, ħ, kB, T, limit::Quantum)
    kBT = kB * T
    f = (freq) -> @. (0.5*ħ*freq) + kBT * log(1 - exp(-ħ*freq/kBT))
    return sum_over_freqs!(ω, f)
end

function F_harmonic(ω, ħ, kB, T, limit::Classical)
    kBT = kB * T
    f = (freq) -> log(ħ*freq/kBT)
    return kBT * sum_over_freqs!(ω, f)
end

function S_harmonic(ω, ħ, kB, T, limit::Quantum)
    kBT = kB * T
    f = (freq) -> ((ħ*freq/kBT) / (exp(ħ*freq/kBT) - 1)) - log(1 - exp(-ħ*freq/kBT))
    return kB * sum_over_freqs!(ω, f)
end

function S_harmonic(ω, ħ, kB, T, limit::Classical)
    kBT = kB * T
    f = (freq) -> (1 - log(ħ*freq/kBT))
    return kB * sum_over_freqs!(ω, f)
end

function Cᵥ_harmonic(ω, kB, T, limit::Quantum)
    tkBT =  2 * kB * T
    f = (freq) -> ((ħ*freq/tkBT)^2) * (csch(ħ*freq/tkBT)^2)
    return kB * sum_over_freqs!(ω, f)
end

function Cᵥ_harmonic(ω, kB, T, limit::Classical)
    return length(ω)*kB
end