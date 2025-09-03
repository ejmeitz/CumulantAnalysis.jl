export Quantum, Classical

abstract type Limit end
struct Quantum <: Limit end
struct Classical <: Limit end

abstract type CumulantEstimator end

struct EffectiveHamiltonian{T} <: CumulantEstimator
    ifc2_path::String
    ifc3_path::String
    ifc4_path::String
    U0::T
end

struct HarmonicEstimator{T} <: CumulantEstimator
    ifc2_path::String
    U0::T
end

struct ResidualEstimator <: CumulantEstimator
    ifc2_path::String
    ifc3_path::String
    ifc4_path::String
end