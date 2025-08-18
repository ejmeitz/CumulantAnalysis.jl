export Quantum, Classical

abstract type Limit end
struct Quantum <: Limit end
struct Classical <: Limit end


function harmonic_properties(
    estim::Estimator,
    ω::AbstractVector,
    kB, ħ
)
    F₀ = F_harmonic(ω, ħ, kB, estim.temperature, limit(estim))
    S₀ = S_harmonic(ω, ħ, kB, estim.temperature, limit(estim))
    U₀ = U_harmonic(ω, ħ, kB, estim.temperature, limit(estim))
    Cᵥ₀ = Cᵥ_harmonic(ω, kB, estim.temperature, limit(estim))

    return F₀, S₀, U₀, Cᵥ₀
end

#* FIX CONTRIBUTION FROM Zero-Point MOTION ON RIGID TRANSLATION MODES??
#* SEE HOW TDEP IMPLEMENTS THINGS LIKE FREE ENERGY

function V_harmonic(ifc2::AbstractMatrix, u::AbstractVector)
    return 0.5 * ((transpose(u) * ifc2) * u)
end

function sum_over_freqs(freqs, f::Function)
    res = 0.0
    for freq in freqs
        if freq > FREQ_TOL
            res += f(freq)
        end
    end
    return res
end

function U_harmonic(ω, ħ, kB, T, limit::Quantum)
    f = (freq) -> (ħ*freq) * ((1 / (exp(ħ*freq/(kB*T)) - 1)) + 0.5)
    return sum_over_freqs(ω, f)
end

function U_harmonic(ω, ħ, kB, T, limit::Classical)
    n_nonzero = count(freq -> freq > FREQ_TOL, ω)
    return n_nonzero*kB*T
end

function F_harmonic(ω, ħ, kB, T, limit::Quantum)
    kBT = kB * T
    f = (freq) -> (0.5*ħ*freq) + kBT * log(1 - exp(-ħ*freq/kBT))
    return sum_over_freqs(ω, f)
end

function F_harmonic(ω, ħ, kB, T, limit::Classical)
    kBT = kB * T
    f = (freq) -> log(ħ*freq/kBT)
    return kBT * sum_over_freqs(ω, f)
end

function S_harmonic(ω, ħ, kB, T, limit::Quantum)
    kBT = kB * T
    f = (freq) -> ((ħ*freq/kBT) / (exp(ħ*freq/kBT) - 1)) - log(1 - exp(-ħ*freq/kBT))
    return kB * sum_over_freqs(ω, f)
end

function S_harmonic(ω, ħ, kB, T, limit::Classical)
    f = (freq) -> (1 - log(ħ*freq/(kB * T)))
    return kB * sum_over_freqs(ω, f)
end

function Cᵥ_harmonic(ω, kB, T, limit::Quantum)
    tkBT =  2 * kB * T
    f = (freq) -> ((ħ*freq/tkBT)^2) * (csch(ħ*freq/tkBT)^2)
    return kB * sum_over_freqs(ω, f)
end

function Cᵥ_harmonic(ω, kB, T, limit::Classical)
    n_nonzero = count(freq -> freq > FREQ_TOL, ω)
    return n_nonzero*kB
end