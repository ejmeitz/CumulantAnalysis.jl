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

function check_ifc_paths(ce::CumulantEstimator)
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
    V₀::T # expects eV
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
get_V₀(e::EffectiveHamiltonianEstimator, V, V₂, V₃, V₄) = e.V₀

function move_ifcs(ehe::EffectiveHamiltonianEstimator, outpath::String)

    check_ifc_paths(ehe)

    new_ifc2_path = joinpath(outpath, "infile.forceconstant")
    new_ifc3_path = joinpath(outpath, "infile.forceconstant_thirdorder")
    new_ifc4_path = joinpath(outpath, "infile.forceconstant_fourthorder")

    isfile(new_ifc2_path) || cp(ehe.ifc2_path, new_ifc2_path; force = true)
    isfile(new_ifc3_path) || cp(ehe.ifc3_path, new_ifc3_path; force = true)
    isfile(new_ifc4_path) || cp(ehe.ifc4_path, new_ifc4_path; force = true)

end

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
get_V₀(::HarmonicEstimator, V, V₂, V₃, V₄) = mean(V .- V₂)

function move_ifcs(ehe::HarmonicEstimator, outpath::String)
    check_ifc_paths(ehe)
    new_ifc2_path = joinpath(outpath, "infile.forceconstant")
    isfile(new_ifc2_path) || cp(ehe.ifc2_path, new_ifc2_path; force = true)
end

########################################


struct FourthOrderEstimator{O,L,C} <: CumulantEstimator{O,L}
    lim::L
    force_calculator::C
    ifc2_path::String
    ifc3_path::String
    ifc4_path::String
    nconf::Int
    n_boot::Int
    boot_size::Int
end

rv(::FourthOrderEstimator, V, V₂, V₃, V₄) = V .- V₂ .- V₃ .- V₄

# Random variable used in nth cumulant
X1(::FourthOrderEstimator, V, V₂, V₃, V₄) = zero.(V₄)
X2(he::FourthOrderEstimator, V, V₂, V₃, V₄) = rv(he, V, V₂, V₃, V₄)
X3(he::FourthOrderEstimator, V, V₂, V₃, V₄) = rv(he, V, V₂, V₃, V₄)

ifc_paths(foe::FourthOrderEstimator) = [foe.ifc2_path, foe.ifc3_path, foe.ifc4_path]
needs_true_V(::FourthOrderEstimator) = true
get_V₀(::FourthOrderEstimator, V, V₂, V₃, V₄) = mean(V .- V₂ .- V₃ .- V₄)

function move_ifcs(foe::FourthOrderEstimator, outpath::String)

    check_ifc_paths(foe)

    new_ifc2_path = joinpath(outpath, "infile.forceconstant")
    new_ifc3_path = joinpath(outpath, "infile.forceconstant_thirdorder")
    new_ifc4_path = joinpath(outpath, "infile.forceconstant_fourthorder")

    isfile(new_ifc2_path) || cp(foe.ifc2_path, new_ifc2_path; force = true)
    isfile(new_ifc3_path) || cp(foe.ifc3_path, new_ifc3_path; force = true)
    isfile(new_ifc4_path) || cp(foe.ifc4_path, new_ifc4_path; force = true)

end

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

rv(::ResidualEstimator, V, V₂, V₃, V₄) = V .- V₂ .- get_V₀(re, V, V₂, V₃, V₄) # R + V₃ + V₄

# Random variable used in nth cumulant
X1(::ResidualEstimator, V, V₂, V₃, V₄) = V₄
X2(re::ResidualEstimator, V, V₂, V₃, V₄) = rv(re, V, V₂, V₃, V₄)
X3(re::ResidualEstimator, V, V₂, V₃, V₄) = rv(re, V, V₂, V₃, V₄)

ifc_paths(re::ResidualEstimator) = [re.ifc2_path, re.ifc3_path, re.ifc4_path]
needs_true_V(::ResidualEstimator) = true
get_V₀(::ResidualEstimator, V, V₂, V₃, V₄) = mean(V .- V₂ .- V₃ .- V₄)

function move_ifcs(re::ResidualEstimator, outpath::String)

    check_ifc_paths(re)

    new_ifc2_path = joinpath(outpath, "infile.forceconstant")
    new_ifc3_path = joinpath(outpath, "infile.forceconstant_thirdorder")
    new_ifc4_path = joinpath(outpath, "infile.forceconstant_fourthorder")

    isfile(new_ifc2_path) || cp(re.ifc2_path, new_ifc2_path; force = true)
    isfile(new_ifc3_path) || cp(re.ifc3_path, new_ifc3_path; force = true)
    isfile(new_ifc4_path) || cp(re.ifc4_path, new_ifc4_path; force = true)

end