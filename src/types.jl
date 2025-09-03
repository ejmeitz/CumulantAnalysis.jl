export Quantum, Classical

abstract type Limit end
struct Quantum <: Limit end
struct Classical <: Limit end

########################################

struct CumulantData{O,A,B,C}
    κ::A
    ∂κ_∂T::B
    ∂²κ_∂T²::C
end

order(::CumulantData{O}) where O = O

########################################

struct BootstrapCumualantEstimate{O,H}
    harmonic::H
    corrections::SVector{O, H}
    correction_SEs::SVector{O, H}
    total::H
    total_SE::H
    property::String
    unit_str::String
end

########################################

abstract type CumulantEstimator{O, L} where {L <: Limit} end

order(::CumulantEstimator{O}) where O = O

struct EffectiveHamiltonianEstimator{O,L,T} <: CumulantEstimator{O,L}
    lim::L
    ifc2_path::String
    ifc3_path::String
    ifc4_path::String
    U0::T
end

########################################

# Random variable used in nth cumulant
X1(::EffectiveHamiltonianEstimator, V, V₂, V₃, V₄) = V₄
X2(::EffectiveHamiltonianEstimator, V, V₂, V₃, V₄) = V₃ .+ V₄
X3(e::EffectiveHamiltonianEstimator, V, V₂, V₃, V₄) = X2(e, V, V₂, V₃, V₄)

struct HarmonicEstimator{O,L,T} <: CumulantEstimator{O,L}
    lim::L
    ifc2_path::String
    # U0::T #! CAN BE CALCULATED FROM DATA
end

########################################

# Random variable used in nth cumulant
X1(::HarmonicEstimator, V, V₂, V₃, V₄) = zero.(V₄)
X2(::HarmonicEstimator, V, V₂, V₃, V₄) = V .- V₂
X3(e::HarmonicEstimator, V, V₂, V₃, V₄) = X2(e, V, V₂, V₃, V₄)


struct ResidualEstimator{O,L} <: CumulantEstimator{O,L}
    lim::L
    ifc2_path::String
    ifc3_path::String
    ifc4_path::String
end

########################################

# Random variable used in nth cumulant
X1(::ResidualEstimator, V, V₂, V₃, V₄) = V₄
function X2(::ResidualEstimator, V, V₂, V₃, V₄)
    U0 = mean(V - V₂ - V₃ - V₄)
    return V .- V₂ .- U0 # R + V₃ + V₄
end
X3(e::ResidualEstimator, V, V₂, V₃, V₄) = X2(e, V, V₂, V₃, V₄)
