export estimate

function calculate_corrections(e::ThermoEstimator, V, ΔV)

    T = ustrip(e.temperature)

    O = order(e)

    ΔF = zeros(O); ΔS = zeros(O)
    ΔU = zeros(O); ΔCᵥ = zeros(O)

    c1 = CumulantData(V, ΔV, kB, T, Val{1}())
    ΔF[1], ΔS[1], ΔU[1], ΔCᵥ[1] = first_order_corrections(c1, T) 

    if O >= 2
        c2 = CumulantData(V, ΔV, kB, T, c1, Val{2}())
        ΔF[2], ΔS[2], ΔU[2], ΔCᵥ[2] = second_order_corrections(c2, kB, T, stochastic(e))
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

function bootstrap_corrections(e::ThermoEstimator, V, ΔV, n_boot, boot_size, ifc_dir::String)
    
    F₀, S₀, U₀, Cᵥ₀ = harmonic_properties(e, ifc_dir)
    O = order(e)

    idx_storage = zeros(Int, boot_size)

    ΔFs = zeros(O, n_boot); ΔSs = zeros(O, n_boot)
    ΔUs = zeros(O, n_boot); ΔCᵥs = zeros(O, n_boot)

    F_totals = zeros(n_boot); S_totals = zeros(n_boot)
    U_totals = zeros(n_boot); Cᵥ_totals = zeros(n_boot)

    # Could parallelize here
    p = Progress(n_boot, "Bootstrapping Corrections")
    for i in 1:n_boot
        sample!(1:length(V), idx_storage; replace = true)
        ΔFs[:,i], ΔSs[:,i], ΔUs[:,i], ΔCᵥs[:,i] = calculate_corrections(e, V[idx_storage], ΔV[idx_storage])
        F_totals[i] = sum(ΔFs[:,i]) + F₀; S_totals[i] = sum(ΔSs[:,i]) + S₀
        U_totals[i] = sum(ΔUs[:,i]) + U₀; Cᵥ_totals[i] = sum(ΔCᵥs[:,i]) + Cᵥ₀
        next!(p)
    end
    finish!(p)

    ΔF = mean(ΔFs, dims = 2); ΔS = mean(ΔSs, dims = 2)
    ΔU = mean(ΔUs, dims = 2); ΔCᵥ = mean(ΔCᵥs, dims = 2)
    
    F_SEs = std(ΔFs, dims = 2); S_SEs = std(ΔSs, dims = 2)
    U_SEs = std(ΔUs, dims = 2); Cᵥ_SEs = std(ΔCᵥs, dims = 2)

    Nat = Int(length(ω) / 3)
    kBNat = ustrip(CumulantAnalysis.kB * Nat)

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

function convert_freq_units(freqs)
    
    freq_unit = unit(first(freqs))
    freq_dim = dimension(freq_unit)

    if freq_dim != Unitful.𝐓^-1
        error(ArgumentError("Frequency units should be inverse time. Got dimensions of $(freq_dim)"))
    end

    return  ustrip.(u"rad/s", freqs)

end

# Ensures that:
# Length Units are Angstrom
# Energy Units are eV
# Mass Units can either be u or g/mol
function convert_units(freqs, ifc2)

    ifc_unit = unit(first(ifc2))
    ifc_dim = dimension(ifc_unit)

    if ifc_unit == NoUnits && freq_unit == NoUnits
        @warn "You did not use Unitful quantities for freuqency or force constants. Assuming frequencies are rad/s and force constants are eV/Å^2"
        return freqs, ifc2
    end

    if ifc_dim != Unitful.𝐌*Unitful.𝐍^-1*Unitful.𝐓^-2 || ifc_dim != Unitful.𝐌*Unitful.𝐓^-2
        error(ArgumentError("Force constant dimensions should be energy / length^2, but got $(ifc_dim)"))
    end

    ifc_units_molar = ifc_dim == Unitful.𝐌*Unitful.𝐍^-1*Unitful.𝐓^-2

    ifc2_ev_ang = ifc_units_molar ? ustrip.(u"eV/Å^2", ifc2 ./ Unitful.Na) : ustrip.(u"eV/Å^2", ifc2)
    freqs_rad_s =  convert_freq_units(freqs)

    return freqs_rad_s, ifc2_ev_ang
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
            V2[i] = (1.5 * ustrip(kB) * parse(Float64, data[4])) * energy_unit
        end
    end

    return ustrip.(V), ustrip.(V2)
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

    V, V2 = get_V(e.cc, calc, ssposcar_path, config_dir, verbose)
    ΔV = V .- V2

    # If lots of configs are made this can take a bit
    t = Threads.@spawn rm(config_dir, recursive = true)

    header = ["V [eV]" "V2 [eV]" "ΔV [eV]"]
    open(joinpath(basedir, "potential_energies.txt"), "w") do f
        writedlm(f, [header; V V2 ΔV])
    end

    res = bootstrap_corrections(e, V, ΔV, n_boot, boot_size, basedir)

    wait(t)

    return res

end




#### LAMMPS ESTIMATOR ####

# # name in lammps dump ==> name used in program
# const header_dict = Dict(
#     "DisplacementX" => "xu",
#     "DisplacementY" => "yu",
#     "DisplacementZ" => "zu",
# )

# function thermo_prop_checks(lammps_dump_path, order, dump_fields)
#     if order ∉ [1,2,3]
#         @error "Can only calculate first and second order cumulant corrections"
#     end

#     @info "Parsing LAMMPS dump file at $(lammps_dump_path)"
#     ld = LammpsDump(lammps_dump_path)
#     @info "Found $(n_samples(ld)) samples and $(n_atoms(ld)) atoms"

#     actual_fields = fields(ld)
#     @assert issubset(dump_fields, actual_fields) "Dump file needs $(dump_fields) fields, got $(actual_fields). Can re-name with dump_fields kwarg"

#     return ld
# end

# function parse_files(lammps_eq_dump_path, stat_file_path, T)

#     @info "Parsing Equilibrium Positions"
#     eq = LammpsDump(lammps_eq_dump_path)

#     parse_timestep!(eq, 1)
#     atom_masses = get_col(eq, "mass")
#     #* HARDCODED COL NAMES!!
#     initial_positions = Matrix(eq.data_storage[!, ["xu", "yu", "zu"]])

#     bulk_properties = readdlm(stat_file_path, comments = true)

#     E_total = bulk_properties[:, 3]
#     V = bulk_properties[:, 4]
#     T_sim = bulk_properties[:, 6]

#     if mean(T_sim) - T > 1e-1
#         @warn "Simulation temperature $(mean(T_sim)) different from target temperature $(T). Results may be inaccurate."
#     end

#     return initial_positions, atom_masses, E_total, V
# end

# # ASSUME LAMMPS UNIT SYSTEM IS METAL
# function estimate(
#     e::LAMMPSEstimator,
#     lammps_dump_path::String,
#     lammps_eq_dump_path::String,
#     stat_file_path::String,
#     ifc2::AbstractMatrix{I},
#     ω::AbstractVector{W};
#     dump_x_unrolled_names::AbstractVector{String} = ["xu", "yu", "zu"],
#     true_F = missing
# ) where {I,W}

#     ifc2, ω = convert_units(ω, ifc2)

#     D = length(dump_x_unrolled_names)
#     N_atoms = Int(length(ω) / 3)

#     ld = thermo_prop_checks(lammps_dump_path, order, dump_x_unrolled_names)

#     initial_positions, _, E_total, V =
#         parse_files(lammps_eq_dump_path, stat_file_path, e.temperature)

#     u = load_displacements(ld, initial_positions, 
#                             dump_x_unrolled_names = dump_x_unrolled_names, D = D)

#     V₂ = V_harmonic.(Ref(ifc2), eachcol(u))
#     ΔV = V .- V₂

#! THIS DOESNT RETURN F0 S0 U0 Cv0 anymore..
#     F₀, ΔF, S₀, ΔS, U₀, ΔU, Cᵥ₀, ΔCᵥ = calculate_corrections(e, ω, V, ΔV)

#     # Estimate true internal energy and heat capacity
#     U_MD = mean(E_total)
#     Cᵥ_MD = var(E_total) / (e.kB * e.temperature^2)

#     true_S = (U_MD - true_F) / e.temperature
#     # we should be able to get elastic moduli and thermal expansion too
#     F_corrections = CumulantCorrections(F₀,
#                                         SVector(ΔF...) ./  N_atoms,
#                                         true_F,
#                                         "F", "[eV/atom]")
#     S_corrections = CumulantCorrections(S₀ / (ustrip(kB)), 
#                                         SVector(ΔS...) ./ (ustrip(kB) * N_atoms),
#                                         true_S,
#                                         "S", "[kB / atom]")
#     U_corrections = CumulantCorrections(U₀,
#                                         SVector(ΔU...) ./  N_atoms,
#                                         U_MD,
#                                         "U", "[eV/atom]")
#     Cv_corrections = CumulantCorrections(Cᵥ₀ / (ustrip(kB)),
#                                          SVector(ΔCᵥ...) ./ (ustrip(kB) * N_atoms),
#                                          Cᵥ_MD,
#                                          "Cv", "[kB / atom]")

#     return F_corrections, S_corrections, U_corrections, Cv_corrections

# end