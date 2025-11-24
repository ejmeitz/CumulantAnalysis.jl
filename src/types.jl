export 
    HarmonicEstimator, 
    FourthOrderEstimator,
    MixedEstimator,
    AnalyticalEstimator


########################################

struct ControlVariateData{T}
    αs::Vector{T}
    μW_estimate::Vector{T} # sample mean of control variates
    var_red::T
end

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

abstract type CumulantEstimator{O} end
abstract type SamplingCumulantEstimator{O} <: CumulantEstimator{O} end

order(::CumulantEstimator{O}) where O = O

function check_ifc_paths(ce::CumulantEstimator)
    for p in ifc_paths(ce)
        isfile(p) || throw(ArgumentError("Force constant path is not a file: $(ifc_path)"))
    end
end

########################################

struct HarmonicEstimator{O} <: SamplingCumulantEstimator{O}
    ifc2_path::String
    nconf::Int
    n_boot::Int
end

function HarmonicEstimator(order::Int, ifc2_path, nconf, n_boot)
    return HarmonicEstimator{order}(ifc2_path, nconf, n_boot)
end

rv(::HarmonicEstimator, V, V₂, V₃, V₄) = V .- V₂
V₀_rv(he::HarmonicEstimator, V, V₂, V₃, V₄) = rv(he, V, V₂, V₃, V₄)


# Random variable used in nth cumulant
X1(::HarmonicEstimator, V, V₂, V₃, V₄) = zero.(V₄)
X2(he::HarmonicEstimator, V, V₂, V₃, V₄) = rv(he, V, V₂, V₃, V₄)
X3(he::HarmonicEstimator, V, V₂, V₃, V₄) = rv(he, V, V₂, V₃, V₄)

ifc_paths(ehe::HarmonicEstimator) = [ehe.ifc2_path]
needs_true_V(::HarmonicEstimator) = true
get_V₀(he::HarmonicEstimator, V, V₂, V₃, V₄) = mean(V₀_rv(he, V, V₂, V₃, V₄))

function move_ifcs(ehe::HarmonicEstimator, outpath::String)
    check_ifc_paths(ehe)
    new_ifc2_path = joinpath(outpath, "infile.forceconstant")
    isfile(new_ifc2_path) || cp(ehe.ifc2_path, new_ifc2_path; force = true)
end

function load_ifcs(::HarmonicEstimator, ucposcar_path::String, basepath::String)

    ifc2_path = joinpath(basepath, "infile.forceconstant")
    ifc2 = read_ifc2(ifc2_path, ucposcar_path)

    return (ifc2, )
end

########################################


struct FourthOrderEstimator{O} <: SamplingCumulantEstimator{O}
    ifc2_path::String
    ifc3_path::String
    ifc4_path::String
    nconf::Int
    n_boot::Int
end

function FourthOrderEstimator(order::Int, ifc2_path, ifc3_path, 
        ifc4_path, nconf, n_boot)
    return FourthOrderEstimator{order}(ifc2_path,
        ifc3_path, ifc4_path, nconf, n_boot)
end


rv(::FourthOrderEstimator, V, V₂, V₃, V₄) = V₃ .+ V₄
V₀_rv(::FourthOrderEstimator, V, V₂, V₃, V₄) = V .- V₂ .- V₃ .- V₄

# Random variable used in nth cumulant
X1(::FourthOrderEstimator, V, V₂, V₃, V₄) = V₄
X2(foe::FourthOrderEstimator, V, V₂, V₃, V₄) = rv(foe, V, V₂, V₃, V₄)
X3(foe::FourthOrderEstimator, V, V₂, V₃, V₄) = rv(foe, V, V₂, V₃, V₄)

ifc_paths(foe::FourthOrderEstimator) = [foe.ifc2_path, foe.ifc3_path, foe.ifc4_path]
needs_true_V(::FourthOrderEstimator) = true
get_V₀(foe::FourthOrderEstimator, V, V₂, V₃, V₄) = mean(V₀_rv(foe, V, V₂, V₃, V₄))

function move_ifcs(foe::FourthOrderEstimator, outpath::String)

    check_ifc_paths(foe)

    new_ifc2_path = joinpath(outpath, "infile.forceconstant")
    new_ifc3_path = joinpath(outpath, "infile.forceconstant_thirdorder")
    new_ifc4_path = joinpath(outpath, "infile.forceconstant_fourthorder")

    isfile(new_ifc2_path) || cp(foe.ifc2_path, new_ifc2_path; force = true)
    isfile(new_ifc3_path) || cp(foe.ifc3_path, new_ifc3_path; force = true)
    isfile(new_ifc4_path) || cp(foe.ifc4_path, new_ifc4_path; force = true)

end

function load_ifcs(::FourthOrderEstimator, ucposcar_path::String, basepath::String)

    ifc2_path = joinpath(basepath, "infile.forceconstant")
    ifc3_path = joinpath(basepath, "infile.forceconstant_thirdorder")
    ifc4_path = joinpath(basepath, "infile.forceconstant_fourthorder")

    ifc2 = read_ifc2(ifc2_path, ucposcar_path)
    ifc3 = read_ifc3(ifc3_path, ucposcar_path)
    ifc4 = read_ifc4(ifc4_path, ucposcar_path)

    return ifc2, ifc3, ifc4
end


###########################

# Uses V0 from MD for Free energy
# Estimates V0 as <V - V2 - V3 - V4>_0 for derivatives
# Approximates V as (V0 + V2 + V3 + V4)
struct MixedEstimator{O,T} <: SamplingCumulantEstimator{O}
    ifc2_path::String
    ifc3_path::String
    ifc4_path::String
    V0::T
    nconf::Int
    n_boot::Int
end

function MixedEstimator(order::Int, ifc2_path, ifc3_path, 
        ifc4_path, V0::T, nconf, n_boot) where T
    return MixedEstimator{order, T}(ifc2_path,
        ifc3_path, ifc4_path, V0, nconf, n_boot)
end


rv(::MixedEstimator, V, V₂, V₃, V₄) = V₃ .+ V₄
V₀_rv(::MixedEstimator, V, V₂, V₃, V₄) = V .- V₂ .- V₃ .- V₄

# Random variable used in nth cumulant
X1(::MixedEstimator, V, V₂, V₃, V₄) = V₄
X2(me::MixedEstimator, V, V₂, V₃, V₄) = rv(me, V, V₂, V₃, V₄)
X3(me::MixedEstimator, V, V₂, V₃, V₄) = rv(me, V, V₂, V₃, V₄)

ifc_paths(me::MixedEstimator) = [me.ifc2_path, me.ifc3_path, me.ifc4_path]
needs_true_V(::MixedEstimator) = true
get_V₀(me::MixedEstimator, V, V₂, V₃, V₄) = me.V0

function move_ifcs(me::MixedEstimator, outpath::String)

    check_ifc_paths(me)

    new_ifc2_path = joinpath(outpath, "infile.forceconstant")
    new_ifc3_path = joinpath(outpath, "infile.forceconstant_thirdorder")
    new_ifc4_path = joinpath(outpath, "infile.forceconstant_fourthorder")

    isfile(new_ifc2_path) || cp(me.ifc2_path, new_ifc2_path; force = true)
    isfile(new_ifc3_path) || cp(me.ifc3_path, new_ifc3_path; force = true)
    isfile(new_ifc4_path) || cp(me.ifc4_path, new_ifc4_path; force = true)

end

function load_ifcs(::MixedEstimator, ucposcar_path::String, basepath::String)

    ifc2_path = joinpath(basepath, "infile.forceconstant")
    ifc3_path = joinpath(basepath, "infile.forceconstant_thirdorder")
    ifc4_path = joinpath(basepath, "infile.forceconstant_fourthorder")

    ifc2 = read_ifc2(ifc2_path, ucposcar_path)
    ifc3 = read_ifc3(ifc3_path, ucposcar_path)
    ifc4 = read_ifc4(ifc4_path, ucposcar_path)

    return ifc2, ifc3, ifc4
end


###########################

# Uses <V4> and <V3*V3> analytically
# Estimates derivatives of V0 numerically
struct AnalyticalEstimator <: CumulantEstimator{2}
    ifc2_path::String
    ifc3_path::String
    ifc4_path::String
    nconf::Int
    n_boot::Int
end

V₀_rv(::AnalyticalEstimator, V, V₂, V₃, V₄) = V .- V₂ .- V₃ .- V₄

ifc_paths(me::AnalyticalEstimator) = [me.ifc2_path, me.ifc3_path, me.ifc4_path]
needs_true_V(::AnalyticalEstimator) = true
get_V₀(me::AnalyticalEstimator, V, V₂, V₃, V₄) = mean(V₀_rv(foe, V, V₂, V₃, V₄))

function move_ifcs(me::AnalyticalEstimator, outpath::String)

    check_ifc_paths(me)

    new_ifc2_path = joinpath(outpath, "infile.forceconstant")
    new_ifc3_path = joinpath(outpath, "infile.forceconstant_thirdorder")
    new_ifc4_path = joinpath(outpath, "infile.forceconstant_fourthorder")

    isfile(new_ifc2_path) || cp(me.ifc2_path, new_ifc2_path; force = true)
    isfile(new_ifc3_path) || cp(me.ifc3_path, new_ifc3_path; force = true)
    isfile(new_ifc4_path) || cp(me.ifc4_path, new_ifc4_path; force = true)

end

function load_ifcs(::AnalyticalEstimator, ucposcar_path::String, basepath::String)

    ifc2_path = joinpath(basepath, "infile.forceconstant")
    ifc3_path = joinpath(basepath, "infile.forceconstant_thirdorder")
    ifc4_path = joinpath(basepath, "infile.forceconstant_fourthorder")

    ifc2 = read_ifc2(ifc2_path, ucposcar_path)
    ifc3 = read_ifc3(ifc3_path, ucposcar_path)
    ifc4 = read_ifc4(ifc4_path, ucposcar_path)

    return ifc2, ifc3, ifc4
end


########################################

# struct EffectiveHamiltonianEstimator{O,L,T} <: SamplingCumulantEstimator{O,L}
#     lim::L
#     ifc2_path::String
#     ifc3_path::String
#     ifc4_path::String
#     V₀::T # expects eV
#     nconf::Int
#     n_boot::Int
# end

# function EffectiveHamiltonianEstimator(order::Int, lim::L, ifc2_path, ifc3_path, 
#         ifc4_path, V₀_eV, nconf, n_boot) where L
#     return EffectiveHamiltonianEstimator{order, L, typeof(V₀_eV)}(lim, ifc2_path,
#         ifc3_path, ifc4_path, V₀_eV, nconf, n_boot)
# end

# rv(::EffectiveHamiltonianEstimator, V, V₂, V₃, V₄) = V₃ .+ V₄

# # Random variable used in nth cumulant
# X1(::EffectiveHamiltonianEstimator, V, V₂, V₃, V₄) = V₄
# X2(ehe::EffectiveHamiltonianEstimator, V, V₂, V₃, V₄) = rv(ehe, V, V₂, V₃, V₄)
# X3(ehe::EffectiveHamiltonianEstimator, V, V₂, V₃, V₄) = rv(ehe, V, V₂, V₃, V₄)

# ifc_paths(ehe::EffectiveHamiltonianEstimator) = [ehe.ifc2_path, ehe.ifc3_path, ehe.ifc4_path]
# needs_true_V(::EffectiveHamiltonianEstimator) = false
# get_V₀(e::EffectiveHamiltonianEstimator, V, V₂, V₃, V₄) = e.V₀

# function move_ifcs(ehe::EffectiveHamiltonianEstimator, outpath::String)

#     check_ifc_paths(ehe)

#     new_ifc2_path = joinpath(outpath, "infile.forceconstant")
#     new_ifc3_path = joinpath(outpath, "infile.forceconstant_thirdorder")
#     new_ifc4_path = joinpath(outpath, "infile.forceconstant_fourthorder")

#     isfile(new_ifc2_path) || cp(ehe.ifc2_path, new_ifc2_path; force = true)
#     isfile(new_ifc3_path) || cp(ehe.ifc3_path, new_ifc3_path; force = true)
#     isfile(new_ifc4_path) || cp(ehe.ifc4_path, new_ifc4_path; force = true)

# end


# function load_ifcs(::EffectiveHamiltonianEstimator, ucposcar_path::String, basepath::String)

#     ifc2_path = joinpath(basepath, "infile.forceconstant")
#     ifc3_path = joinpath(basepath, "infile.forceconstant_thirdorder")
#     ifc4_path = joinpath(basepath, "infile.forceconstant_fourthorder")

#     ifc2 = read_ifc2(ifc2_path, ucposcar_path)
#     ifc3 = read_ifc3(ifc3_path, ucposcar_path)
#     ifc4 = read_ifc4(ifc4_path, ucposcar_path)

#     return ifc2, ifc3, ifc4
# end