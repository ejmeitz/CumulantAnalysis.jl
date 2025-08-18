module TDEPExt

using TDEP
using CumulantAnalysis
using AtomsCalculators
using OhMyThreads
using StaticArrays

function CumulantAnalysis.sTDEPEstimator(order, nsamples, temperature_K, quamtum::Bool)
    limit = (quantum == true) ? Quantum : Classical
    T = ustrip(temperature_K) * u"K"
    cc = CanonicalConfiguration(
        temperature = temperature,
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

    p = Progress(cc.nconf, desc = "Calculating Energies")
    @tasks for i in 1:cc.nconf
        @set ntasks = n_threads

        filepath = get_filepath(i)
        read_poscar_positions(S.positions, filepath; n_atoms = n_atoms)
       
        V[i] = uconvert(energy_unit, AtomsCalculators.potential_energy(S, calc))

        next!(p)
    end
    finish!(p)

    # Parse V2 from logfile
    open(joinpath(basedir, "$(TDEP.cmd_name(cmd)).log")) do f
        readline(f)
        readline(f)
        readline(f)
        for i in 1:cc.nconf
            data = split(strip(readline(f)))
            V2[i] = (n_atoms * parse(Float64, data[4])) * energy_unit
        end
    end

    return V, V2
end

function CumulantAnalysis.estimate(
    e::sTDEPEstimator,
    calc,
    ω::AbstractVector{W},
    ssposcar_path::String;
    basedir = dirname(ssposcar_path),
    n_threads = Threads.nthreads(),
    true_F = missing,
    verbose::Bool = true
)

    ω = TDEP.convert_freq_units(ω)

    if !isfile(joinpath(basedir, "infile.forceconstant"))
        raise(ArgumentError("Could not find infile.forceconstant in basedir: $(basedir)"))
    end
    
    V, V2 = get_V(e.cc, calc, ssposcar_path, basedir, verbose, n_threads)
    ΔV = V .- V2

    F₀, ΔF, S₀, ΔS, U₀, ΔU, Cᵥ₀, ΔCᵥ = calculate_corrections(e, ω, V, ΔV)

    # we should be able to get elastic moduli and thermal expansion too
    F_corrections = CumulantCorrections(F₀, SVector(ΔF...), true_F, "F")
    S_corrections = CumulantCorrections(S₀,  SVector(ΔS...), missing, "S")
    U_corrections = CumulantCorrections(U₀,  SVector(ΔU...), missing, "U")
    Cv_corrections = CumulantCorrections(Cᵥ₀,  SVector(ΔCᵥ...), missing, "Cv")

    return F_corrections, S_corrections, U_corrections, Cv_corrections

end

end