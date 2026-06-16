export crystal_thermodynamic_properties, make_stdep_ifcs

parse_path(path_func::Function, T) = path_func(T)
parse_path(path::String, T) = path 

"""
Calculate thermodynamic properties of a crystal across a
range of temperatures using the free energy cumulant expansion.

All path-like parameters may be passed as a single string, or as a function 
which takes on parameters (temperature) and returns a string. For example,
'''julia
make_ucposcar_path = (T) -> joinpath("/home/user", "infile.ucposcar$(T)")
ucposcar_path = "/home/user/infile.ucposcar"
'''

Parameters:
-----------
- `temperatures::AbstractVector{<:Real}`: The temperatures at which to calculate the thermodynamic properties.
- `outpath::Union{Function, String}`: The output path for the results.
- `ucposcar_path::Union{Function, String}`: The path/function to the POSCAR file for the unit cell.
- `ssposcar_path::Union{Function, String}`: The path/function to the POSCAR file for the supercell.
- `ifc2_path::Union{Function, String}`: The path/function to the infile.forceconstant file (TDEP format).
- `ifc3_path::Union{Function, String}`: The path/function to the infile.forceconstant file (TDEP format).
- `ifc4_path::Union{Function, String}`: The path/function to the infile.forceconstant file (TDEP format).
- `pot_cmds::Union{String, Vector{String}}`: The LAMMPS potential commands to calculate the potential energy. For example,
    `pot_cmds = pot_cmds = ["pair_style lj/cut 6.955", "pair_coeff * * 0.0032135 2.782", "pair_modify shift yes"]`
- `nconf::Int = 100_000`: The number of configurations to sample for the constant correction.
- `nboot::Int = 2500`: The number of bootstrap samples to use for error estimation of the constant correction.
- `size_study:Bool = false`: Whether to calculate the constant correction as a function of the number of samples (log-spaced).
- `harmonic_q_mesh::AbstractVector{<:Integer} = [30, 30, 30]`: The q-mesh to use for the harmonic contribution.
- `free_energy_q_mesh::AbstractVector{<:Integer} = [25, 25, 25]`: The q-mesh to use for the cumulant corrections.
- `n_threads::Integer = Threads.nthreads()`: The number of threads to use for parallelization.
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
)

    @info "Calculating Thermodynamic Properties with $n_threads threads."

    if n_threads > Threads.nthreads()
        @warn "You asked for $n_threads threads, but only $(Threads.nthreads()) are available. Be sure to set JULIA_NUM_THREADS in your environment."
    end

    for (i,T) in enumerate(temperatures)
        @info "Running T = $(T) K"

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

        estimate(
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
            use_hot = false
        )

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

    #! TODO FIT 3rd/4th ORDER IFCS

    #! RETURN F_H FOR EACH ITERATION TO CHECK CONVERGENCE

end
