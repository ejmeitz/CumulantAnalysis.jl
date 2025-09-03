export estimate

function calculate_corrections(T, O, V, ΔV, is_stochastic::Bool)

    ΔF = zeros(O); ΔS = zeros(O)
    ΔU = zeros(O); ΔCᵥ = zeros(O)

    c1 = CumulantData(V, ΔV, kB, T, Val{1}())
    ΔF[1], ΔS[1], ΔU[1], ΔCᵥ[1] = first_order_corrections(c1, T) 

    if O >= 2
        c2 = CumulantData(V, ΔV, kB, T, c1, Val{2}())
        ΔF[2], ΔS[2], ΔU[2], ΔCᵥ[2] = second_order_corrections(c2, kB, T, is_stochastic)
    end

    if O >= 3
        c3 = CumulantData(V, ΔV, kB, T, c1, Val{3}())
        ΔF[3], ΔS[3], ΔU[3], ΔCᵥ[3] = third_order_corrections(c3, kB, T)
    end

    return ΔF, ΔS, ΔU, ΔCᵥ

end

struct BootstrapCumualantEstimate{O,H}
    harmonic::H
    corrections::SVector{O, H}
    correction_SEs::SVector{O, H}
    total::H
    total_SE::H
    property::String
    unit_str::String
end

function bootstrap_corrections(T, V, ΔV, n_boot, boot_size, ifc_dir::String,
                                 Nat::Int, O::Int, is_stochastic::Bool, ::Type{L}) where {L <: Limit}

    # these are returned per-atom
    F₀, S₀, U₀, Cᵥ₀ = harmonic_properties(T, L, ifc_dir)

    idx_storage = zeros(Int, boot_size)

    ΔFs = zeros(O, n_boot); ΔSs = zeros(O, n_boot)
    ΔUs = zeros(O, n_boot); ΔCᵥs = zeros(O, n_boot)

    p = Progress(n_boot, "Bootstrapping Corrections")
    for i in 1:n_boot
        sample!(1:length(V), idx_storage; replace = true)
        ΔFs[:,i], ΔSs[:,i], ΔUs[:,i], ΔCᵥs[:,i] = calculate_corrections(T, O, V[idx_storage], ΔV[idx_storage], is_stochastic)
        next!(p)
    end
    finish!(p)

    F_totals = sum(ΔFs, dims = 1) .+ (F₀*Nat)
    S_totals = sum(ΔSs, dims = 1) .+ (S₀*Nat)
    U_totals = sum(ΔUs, dims = 1) .+ (U₀*Nat)
    Cᵥ_totals = sum(ΔCᵥs, dims = 1) .+ (Cᵥ₀*Nat)

    ΔF = mean(ΔFs, dims = 2); ΔS = mean(ΔSs, dims = 2)
    ΔU = mean(ΔUs, dims = 2); ΔCᵥ = mean(ΔCᵥs, dims = 2)
    
    F_SEs = std(ΔFs, dims = 2); S_SEs = std(ΔSs, dims = 2)
    U_SEs = std(ΔUs, dims = 2); Cᵥ_SEs = std(ΔCᵥs, dims = 2)

    kBNat = ustrip(CumulantAnalysis.kB) * Nat

    F = BootstrapCumualantEstimate(
        F₀, SVector(ΔF...) ./ Nat, SVector(F_SEs...) ./ Nat,
        mean(F_totals) / Nat, std(F_totals) / Nat, "F", "[eV/atom]"
    )

    S = BootstrapCumualantEstimate(
        S₀ / kB, SVector(ΔS...) ./ kBNat, SVector(S_SEs...) ./ kBNat,
        mean(S_totals) / kBNat, std(S_totals) / kBNat, "S", "[kB / atom]"
    )

    U = BootstrapCumualantEstimate(
        U₀, SVector(ΔU...) ./ Nat, SVector(U_SEs...) ./ Nat,
        mean(U_totals) / Nat, std(U_totals) / Nat, "U", "[eV/atom]"
    )

    Cᵥ = BootstrapCumualantEstimate(
        Cᵥ₀ / kB, SVector(ΔCᵥ...) ./ kBNat, SVector(Cᵥ_SEs...) ./ kBNat,
        mean(Cᵥ_totals) / kBNat, std(Cᵥ_totals) / kBNat, "Cv", "[kB / atom]"
    )

    return F, S, U, Cᵥ

end

#### sTDEP ESTIMATOR ####

function get_V(cc, calc, ssposcar_path, basedir, verbose)

    sys_ss = TDEPSystem(ssposcar_path)
    n_atoms = length(sys_ss)

    get_filepath = (i) -> joinpath(basedir, "contcar_conf$(lpad(i, 4, '0'))")
    energy_unit = (AtomsCalculators.energy_unit(calc) == NoUnits) ? NoUnits : u"eV"

    if energy_unit == NoUnits
        @warn "Your energy calculator did not have units. Assuming eV."
    end

    @info "Generating Configurations"
    execute(cc, basedir, 1, verbose)

    V = zeros(typeof(1.0 * energy_unit), cc.nconf)
    V2 = zeros(typeof(1.0 * energy_unit), cc.nconf)

    posns = zeros(Float64, 3, n_atoms)

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
            V2[i] = (1.5 * ustrip(kB) * parse(Float64, data[4])) * energy_unit
        end
    end

    return ustrip.(V), ustrip.(V2), n_atoms
end

function estimate(
    e::sTDEPEstimator,
    calc,
    basedir::String;
    ucposcar_path::String = joinpath(basedir, "infile.ucposcar"),
    ssposcar_path::String = joinpath(basedir, "infile.ssposcar"),
    ifc_path::String = joinpath(basedir, "infile.forceconstant"),
    verbose::Bool = true,
    n_boot::Int = 100,
    boot_size::Int = 10_000
)

    if e.cc.nconf < boot_size
        error("Number of configurations $(e.cc.nconf) is less than the number of bootstrap samples $(boot_size).")
    end

    config_dir = joinpath(basedir, "configs")
    mkpath(config_dir)

    isfile(ifc_path) || error(ArgumentError("Could not find infile.forceconstant at $(basedir)"))
    isfile(ucposcar_path) || error(ArgumentError("Could not find infile.ucposcar at $(ucposcar_path)"))
    isfile(ssposcar_path) || error(ArgumentError("Could not find infile.ssposcar at $(ssposcar_path)"))

    cp(ifc_path, joinpath(config_dir, "infile.forceconstant"); force = true)
    cp(ifc_path, joinpath(basedir, "infile.forceconstant"); force = true)
    cp(ucposcar_path, joinpath(config_dir, "infile.ucposcar"); force = true)
    cp(ucposcar_path, joinpath(basedir, "infile.ucposcar"); force = true)
    cp(ssposcar_path, joinpath(config_dir, "infile.ssposcar"); force = true)

    V, V2, n_atoms = get_V(e.cc, calc, ssposcar_path, config_dir, verbose)
    ΔV = V .- V2

    # If lots of configs are made this can take a bit
    t = Threads.@spawn rm(config_dir, recursive = true)

    header = ["V [eV]" "V2 [eV]" "ΔV [eV]"]
    open(joinpath(basedir, "potential_energies.txt"), "w") do f
        writedlm(f, [header; V V2 ΔV])
    end

    #! Use V2 in place of V for derivatives??
    T = Float64(ustrip(e.temperature))
    res = bootstrap_corrections(T, V2, ΔV, n_boot, boot_size, basedir, n_atoms, order(e), stochastic(e), limit(e))

    wait(t)

    return res

end



