module TDEPExt

using TDEP
using CumulantAnalysis
using AtomsCalculators
using OhMyThreads
using StaticArrays
using ProgressMeter
using DelimitedFiles
using StatsBase

function CumulantAnalysis.sTDEPEstimator(order::Int, nsamples::Int, temperature_K, quantum::Bool)
    limit = (quantum == true) ? Quantum : Classical
    T = ustrip(temperature_K) * u"K"
    cc = CanonicalConfiguration(
        temperature = Float64(ustrip(temperature_K)),
        nconf = nsamples,
        quantum = quantum, 
    )
    return sTDEPEstimator{order, limit, typeof(T), typeof(cc)}(cc, T)
end

function get_V(cc, calc, ssposcar_path, basedir, verbose, n_threads)

    sys_ss = TDEPSystem(ssposcar_path)
    n_atoms = length(sys_ss)

    get_filepath = (i) -> joinpath(basedir, "contcar_conf$(lpad(i, 4, '0'))")
    energy_unit = (AtomsCalculators.energy_unit(calc) == NoUnits) ? NoUnits : u"eV"

    if energy_unit == NoUnits
        @warn "Your calculator did not have units. Assuming eV."
    end

    @info "Generating Configurations"
    execute(cc, basedir, 1, verbose)

    V = zeros(typeof(1.0 * energy_unit), cc.nconf)
    V2 = zeros(typeof(1.0 * energy_unit), cc.nconf)

    posns = Matrix{Float64}(undef, 3, n_atoms)

    p = Progress(cc.nconf, desc = "Calculating Energies")
    L = typeof(1.0u"Å")
    for i in 1:cc.nconf   
        filepath = get_filepath(i)
        TDEP.read_poscar_positions!(posns, filepath; n_atoms = n_atoms)
        
        new_sys = TDEPSystem(sys_ss, vec(collect(reinterpret(SVector{3, L}, posns))))

        V[i] = uconvert(energy_unit, AtomsCalculators.potential_energy(new_sys, calc))

        next!(p)
    end
    finish!(p)

    # Parse V2 from logfile
    open(joinpath(basedir, "$(TDEP.cmd_name(cc)).log")) do f

        while strip(readline(f)) != "... remapped fc" end
        readline(f)

        for i in 1:cc.nconf
            data = split(strip(readline(f)))
            V2[i] = (1.5 * ustrip(CumulantAnalysis.kB) * parse(Float64, data[4])) * energy_unit
        end
    end

    return ustrip.(V), ustrip.(V2)
end

function CumulantAnalysis.estimate(
    e::sTDEPEstimator,
    calc,
    ω,
    basedir::String;
    ucposcar_path::String = joinpath(basedir, "infile.ucposcar"),
    ssposcar_path::String = joinpath(basedir, "infile.ssposcar"),
    ifc_path::String = joinpath(basedir, "infile.forceconstant"),
    n_threads = Threads.nthreads(),
    verbose::Bool = true,
    n_boot::Int = 100,
    boot_size::Int = 10_000
)
    if e.cc.nconf < boot_size
        error("Number of configurations $(e.cc.nconf) is less than the number of bootstrap samples $(boot_size).")
    end

    config_dir = joinpath(basedir, "configs")
    mkpath(config_dir)

    ω = CumulantAnalysis.convert_freq_units(ω)

    if !isfile(ifc_path)
        error(ArgumentError("Could not find infile.forceconstant in basedir: $(basedir)"))
    end

    cp(ifc_path, joinpath(config_dir, "infile.forceconstant"); force = true)
    cp(ucposcar_path, joinpath(config_dir, "infile.ucposcar"); force = true)
    cp(ssposcar_path, joinpath(config_dir, "infile.ssposcar"); force = true)

    V, V2 = get_V(e.cc, calc, ssposcar_path, config_dir, verbose, n_threads)
    ΔV = V .- V2

    header = ["V [eV]" "V2 [eV]" "ΔV [eV]"]
    open(joinpath(basedir, "potential_energies.txt"), "w") do f
        writedlm(f, [header; V V2 ΔV])
    end

    return CumulantAnalysis.bootstrap_corrections(e, ω, V, ΔV, n_boot, boot_size; normalize = true)

end

end