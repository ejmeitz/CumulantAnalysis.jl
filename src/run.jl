export crystal_thermodynamic_properties, make_stdep_ifcs

parse_path(path_func::Function, T) = path_func(T)
parse_path(path::String, T) = path 

"""
"""
function crystal_thermodynamic_properties(
    temperatures::AbstractVector{<:Real},
    outpath::Union{Function, String},
    ucposcar_path::Union{Function, String},
    ssposcar_path::Union{Function, String},
    ifc2_path::Union{Function, String},
    ifc3_path::Union{Function, String},
    ifc4_path::Union{Function, String},
    pot_cmds::Union{String, Vector{String}};
    quantum::Bool = false,
    nconf::Int = 100_000,
    nboot::Int = 2500,
    size_study::Bool = false,
    harmonic_q_mesh::AbstractVector{<:Integer} = [30,30,30],
    free_energy_q_mesh::AbstractVector{<:Integer} = [25,25,25],
    n_threads::Integer = Threads.nthreads(),
    kwargs...
)
    all_ifcs = Vector{IFC2}(undef, length(temperatures))
    all_ucs = Vector{CrystalStructure}(undef, length(temperatures))

    LIMIT = quantum ? Quantum : Classical

    for (i,T) in enumerate(temperatures)

        outpath_T = parse_path(outpath, T)
        ifc2_path_T = parse_path(ifc2_path, T)
        ifc3_path_T = parse_path(ifc3_path, T)
        ifc4_path_T = parse_path(ifc4_path, T)
        ucposcar_path_T = parse_path(ucposcar_path, T)
        ssposcar_path_T = parse_path(ssposcar_path, T)

        mkpath(outpath_T)

        estim = AnalyticalEstimator(
                ifc2_path_T, ifc3_path_T, ifc4_path_T, nconf, nboot
            )

        res, all_ifcs[i] = estimate(
            estim,
            Float64(T),
            outpath_T,
            pot_cmds;
            ucposcar_path = ucposcar_path_T,
            ssposcar_path = ssposcar_path_T,
            n_threads = n_threads,
            size_study = size_study,
            quantum = quantum,
            harmonic_q_mesh = harmonic_q_mesh,
            free_energy_q_mesh = free_energy_q_mesh,
        )

        all_ucs[i] = CrystalStructure(ucposcar_path_T)

    end
end

"""
"""
function make_stdep_ifcs(
    ucposcar_path::String,
    ssposcar_path::String,
    outdir::String,
    pot_cmds::Vector{String},
    n_iter::Int,
    r_cut::Real,
    T::Real,
    maximum_frequency::Real,
    quantum::Bool,
    kwargs...
)

    sys = CrystalStructure(ssposcar_path)
    lc = LAMMPSCalculator(sys, pot_cmds)

    mkpath(outdir)
    cp(ucposcar_path, joinpath(outdir, "infile.ucposcar"); force = true)
    cp(ssposcar_path, joinpath(outdir, "infile.ssposcar"); force = true)

    sTDEP(
        sys,
        lc,
        outdir,
        n_iter,
        r_cut,
        T,
        maximum_frequency;
        quantum = quantum,
        kwargs...
    )

end
