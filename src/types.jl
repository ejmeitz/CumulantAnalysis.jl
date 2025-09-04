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

abstract type CumulantEstimator{O, L <: Limit} end

order(::CumulantEstimator{O}) where O = O

function check_paths(ce::CumulantEstimator)
    for p in ifc_paths(ce)
        isfile(p) || throw(ArgumentError("Force constant path is not a file: $(ifc_path)"))
    end
end

########################################

struct EffectiveHamiltonianEstimator{O,L,T} <: CumulantEstimator{O,L}
    lim::L
    ifc2_path::String
    ifc3_path::String
    ifc4_path::String
    U0::T
    nconf::Int
    n_boot::Int
    boot_size::Int
end

rv(::EffectiveHamiltonianEstimator, V, V₂, V₃, V₄) = V₃ .+ V₄

# Random variable used in nth cumulant
X1(::EffectiveHamiltonianEstimator, V, V₂, V₃, V₄) = V₄
X2(ehe::EffectiveHamiltonianEstimator, V, V₂, V₃, V₄) = rv(ehe, V, V₂, V₃, V₄)
X3(ehe::EffectiveHamiltonianEstimator, V, V₂, V₃, V₄) = rv(ehe, V, V₂, V₃, V₄)

ifc_paths(ehe::EffectiveHamiltonianEstimator) = [ehe.ifc2_path, ehe.ifc3_path, ehe.ifc4_path]
needs_true_V(::EffectiveHamiltonianEstimator) = false
U0(e::EffectiveHamiltonianEstimator, V, V₂, V₃, V₄) = e.U0

########################################


struct HarmonicEstimator{O,L,C} <: CumulantEstimator{O,L}
    lim::L
    force_calculator::C
    ifc2_path::String
    nconf::Int
    n_boot::Int
    boot_size::Int
end

rv(::HarmonicEstimator, V, V₂, V₃, V₄) = V .- V₂

# Random variable used in nth cumulant
X1(::HarmonicEstimator, V, V₂, V₃, V₄) = zero.(V₄)
X2(he::HarmonicEstimator, V, V₂, V₃, V₄) = rv(he, V, V₂, V₃, V₄)
X3(he::HarmonicEstimator, V, V₂, V₃, V₄) = rv(he, V, V₂, V₃, V₄)

ifc_paths(ehe::HarmonicEstimator) = [ehe.ifc2_path]
needs_true_V(::HarmonicEstimator) = true
U0(::HarmonicEstimator, V, V₂, V₃, V₄) = mean(V .- V₂)

########################################

struct ResidualEstimator{O,L,C} <: CumulantEstimator{O,L}
    lim::L
    force_calculator::C
    ifc2_path::String
    ifc3_path::String
    ifc4_path::String
    nconf::Int
    n_boot::Int
    boot_size::Int
end

rv(::ResidualEstimator, V, V₂, V₃, V₄) = V .- V₂ .- U0(re, V, V₂, V₃, V₄) # R + V₃ + V₄

# Random variable used in nth cumulant
X1(::ResidualEstimator, V, V₂, V₃, V₄) = V₄
X2(re::ResidualEstimator, V, V₂, V₃, V₄) = rv(re, V, V₂, V₃, V₄)
X3(re::ResidualEstimator, V, V₂, V₃, V₄) = rv(re, V, V₂, V₃, V₄)

ifc_paths(ehe::ResidualEstimator) = [ehe.ifc2_path, ehe.ifc3_path, ehe.ifc4_path]
needs_true_V(::ResidualEstimator) = true
U0(::ResidualEstimator, V, V₂, V₃, V₄) = mean(V .- V₂ .- V₃ .- V₄)