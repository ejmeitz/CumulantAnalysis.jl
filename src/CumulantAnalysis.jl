module CumulantAnalysis

using DataFrames
using AtomsCalculators
using StaticArrays
using Statistics
using StatsBase
using Unitful
using DelimitedFiles
using Measurements
using ProgressMeter
using HDF5
using TDEP
using FileIO
using Printf
using OrderedCollections

@static if !CumulantAnalysis.TDEP.TDEP_jll.is_available()
    @warn "Could not load TDEP on this platform. Try Linux or MacOs"
end

const FREQ_TOL = 1e-6
const kB = ustrip(u"eV / K", Unitful.k)
const ħ = ustrip(u"eV * s", Unitful.ħ)

include("types.jl")
include("harmonic_properties.jl")
include("cumulant_corrections.jl")
include("estimate.jl")


end # CumulantAnalysis