module CumulantAnalysis

using DataFrames
using StaticArrays
using Statistics
using Unitful
using DelimitedFiles
using Measurements
using ProgressMeter
using HDF5
using FileIO
using Printf
using OrderedCollections

const FREQ_TOL = 1e-6

include("lammps_dump_parser.jl")
include("harmonic_properties.jl")
include("block_averaging.jl")
include("error_analysis.jl")
include("cumulant_corrections.jl")
include("phonon.jl")


end # CumulantAnalysis