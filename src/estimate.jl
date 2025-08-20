export estimate

function calculate_corrections(e::ThermoEstimator, ω::AbstractVector, V, ΔV)

    T = ustrip(e.temperature)

    F₀, S₀, U₀, Cᵥ₀ = harmonic_properties(e, ω, ustrip(kB), ustrip(ħ))
    O = order(e)

    ΔF = zeros(O); ΔS = zeros(O)
    ΔU = zeros(O); ΔCᵥ = zeros(O)

    @info "Calculating First Order Corrections"
    c1 = CumulantData(V, ΔV, kB, T, Val{1}())
    ΔF[1], ΔS[1], ΔU[1], ΔCᵥ[1] = first_order_corrections(c1, T) 

    if O >= 2
        @info "Calculating Secoind Order Corrections"
        c2 = CumulantData(V, ΔV, kB, T, c1, Val{2}())
        ΔF[2], ΔS[2], ΔU[2], ΔCᵥ[2] = second_order_corrections(c2, kB, T, stochastic(e))
    end

    if O >= 3
        @info "Calculating Third Order Corrections"
        c3 = CumulantData(V, ΔV, kB, T, c1, Val{3}())
        ΔF[3], ΔS[3], ΔU[3], ΔCᵥ[3] = third_order_corrections(c3, kB, T)
    end

    return F₀, ΔF, S₀, ΔS, U₀, ΔU, Cᵥ₀, ΔCᵥ

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

#### LAMMPS ESTIMATOR ####

# name in lammps dump ==> name used in program
const header_dict = Dict(
    "DisplacementX" => "xu",
    "DisplacementY" => "yu",
    "DisplacementZ" => "zu",
)

function thermo_prop_checks(lammps_dump_path, order, dump_fields)
    if order ∉ [1,2,3]
        @error "Can only calculate first and second order cumulant corrections"
    end

    @info "Parsing LAMMPS dump file at $(lammps_dump_path)"
    ld = LammpsDump(lammps_dump_path)
    @info "Found $(n_samples(ld)) samples and $(n_atoms(ld)) atoms"

    actual_fields = fields(ld)
    @assert issubset(dump_fields, actual_fields) "Dump file needs $(dump_fields) fields, got $(actual_fields). Can re-name with dump_fields kwarg"

    return ld
end

function parse_files(lammps_eq_dump_path, stat_file_path, T)

    @info "Parsing Equilibrium Positions"
    eq = LammpsDump(lammps_eq_dump_path)

    parse_timestep!(eq, 1)
    atom_masses = get_col(eq, "mass")
    #* HARDCODED COL NAMES!!
    initial_positions = Matrix(eq.data_storage[!, ["xu", "yu", "zu"]])

    bulk_properties = readdlm(stat_file_path, comments = true)

    E_total = bulk_properties[:, 3]
    V = bulk_properties[:, 4]
    T_sim = bulk_properties[:, 6]

    if mean(T_sim) - T > 1e-1
        @warn "Simulation temperature $(mean(T_sim)) different from target temperature $(T). Results may be inaccurate."
    end

    return initial_positions, atom_masses, E_total, V
end

# ASSUME LAMMPS UNIT SYSTEM IS METAL
function estimate(
    e::LAMMPSEstimator,
    lammps_dump_path::String,
    lammps_eq_dump_path::String,
    stat_file_path::String,
    ifc2::AbstractMatrix{I},
    ω::AbstractVector{W};
    dump_x_unrolled_names::AbstractVector{String} = ["xu", "yu", "zu"],
    true_F = missing
) where {I,W}

    ifc2, ω = convert_units(ω, ifc2)

    D = length(dump_x_unrolled_names)

    ld = thermo_prop_checks(lammps_dump_path, order, dump_x_unrolled_names)

    initial_positions, _, E_total, V =
        parse_files(lammps_eq_dump_path, stat_file_path, e.temperature)

    u = load_displacements(ld, initial_positions, 
                            dump_x_unrolled_names = dump_x_unrolled_names, D = D)

    V₂ = V_harmonic.(Ref(ifc2), eachcol(u))
    ΔV = V .- V₂

    F₀, ΔF, S₀, ΔS, U₀, ΔU, Cᵥ₀, ΔCᵥ = calculate_corrections(e, ω, V, ΔV)

    # Estimate true internal energy and heat capacity
    U_MD = mean(E_total)
    Cᵥ_MD = var(E_total) / (e.kB * e.temperature^2)

    true_S = (U_MD - true_F) / e.temperature
    # we should be able to get elastic moduli and thermal expansion too
    F_corrections = CumulantCorrections(F₀, SVector(ΔF...), true_F, "F")
    S_corrections = CumulantCorrections(S₀,  SVector(ΔS...), true_S, "S")
    U_corrections = CumulantCorrections(U₀,  SVector(ΔU...), U_MD, "U")
    Cv_corrections = CumulantCorrections(Cᵥ₀,  SVector(ΔCᵥ...), Cᵥ_MD, "Cv")

    return F_corrections, S_corrections, U_corrections, Cv_corrections

end