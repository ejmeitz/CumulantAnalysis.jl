module CumulantAnalysis

using DataFrames
using Statistics
using Unitful
using DelimitedFiles
using Measurements
using ProgressMeter
using Plots

#Supresses Plot Displays
ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"

include("lammps_dump_parser.jl")
include("harmonic_properties.jl")
include("block_averaging.jl")
include("cumulant_corrections.jl")
include("phonon.jl")


end # CumulantAnalysis