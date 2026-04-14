module CumulantAnalysis

using DataFrames
using AtomsCalculators
using StaticArrays
using Statistics
using StatsBase
using Unitful
using DelimitedFiles
using ProgressMeter
using HDF5
using FileIO
using LinearAlgebra
using Printf
using OrderedCollections
using Random
using LAMMPS
using SavitzkyGolay
using LatticeDynamicsToolkit
import LatticeDynamicsToolkit: Hartree_to_eV


const kB = ustrip(u"eV / K", Unitful.k)
const ħ = ustrip(u"eV * s", Unitful.ħ)

include("types.jl")
include("cumulant_corrections.jl")
include("bootstrap.jl")
include("estimate.jl")
include("run.jl")

end # CumulantAnalysis
