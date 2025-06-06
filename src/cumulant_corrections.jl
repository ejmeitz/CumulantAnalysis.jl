export estimate_thermo_properties, Cumulants

struct Cumulants{A,B,C,D,E,F,G,H,I}
    κ₁::A
    ∂κ₁_∂T::B
    ∂²κ₁_∂T²::C
    κ₂::D
    ∂κ₂_∂T::E
    ∂²κ₂_∂T²::F
    κ₃::G
    ∂κ₃_∂T::H
    ∂²κ₃_∂T²::I
end


function Cumulants(V, ΔV, kB, T)

    ΔV² = ΔV .^ 2
    ΔV³ = ΔV² .* ΔV
    
    # Various derivatives of ⟨O⟩
    ∂A_∂T(A) = cov(A, V) / (kB * T * T)
    ∂AB_∂T(A, B; dA = ∂A_∂T(A), dB =  ∂A_∂T(B)) = (mean(A) * dB) + (mean(B) * dA)
    ∂²A_∂T²(A; dA = ∂A_∂T(A)) = (-2*dA/T) + ((1/(kB*T*T)) * (∂A_∂T(A.*V) - ∂AB_∂T(A, V; dA = dA)))

    κ₁ = mean(ΔV)
    ∂κ₁_∂T = ∂A_∂T(ΔV)
    ∂²κ₁_∂T² = ∂²A_∂T²(ΔV; dA = ∂κ₁_∂T)

    κ₂ = var(ΔV)
    ∂ΔV²_∂T = ∂A_∂T(ΔV²)
    ∂²ΔV²_∂T² = ∂²A_∂T²(ΔV²; dA = ∂ΔV²_∂T)
    ∂κ₂_∂T = ∂ΔV²_∂T - (2*κ₁*∂κ₁_∂T)
    ∂²κ₂_∂T² = ∂²ΔV²_∂T² - 2*((∂κ₁_∂T^2) + (κ₁*∂²κ₁_∂T²))

    κ₃ = mean((ΔV .- κ₁).^3)
    ∂ΔV³_∂T = ∂A_∂T(ΔV³)
    μ_ΔV² = mean(ΔV²)
    ∂κ₃_∂T = ∂ΔV³_∂T - 3*κ₁*∂ΔV²_∂T + 3*μ_ΔV²*∂κ₁_∂T
    ∂²κ₃_∂T² = ∂²A_∂T²(ΔV³, dA = ∂ΔV³_∂T) - 3*(∂κ₁_∂T*∂ΔV²_∂T + κ₁*∂²ΔV²_∂T²) + 3*(∂ΔV²_∂T*∂κ₁_∂T + μ_ΔV²*∂²κ₁_∂T²)

    return Cumulants(κ₁, ∂κ₁_∂T, ∂²κ₁_∂T², κ₂, ∂κ₂_∂T, ∂²κ₂_∂T², κ₃, ∂κ₃_∂T, ∂²κ₃_∂T²)
end

function first_order(c::Cumulants, T)
   
    F_correction = c.κ₁
    S_correction = -c.∂κ₁_∂T
    U_correction = c.κ₁ - T*c.∂κ₁_∂T
    Cv_correction = -T*c.∂²κ₁_∂T²

    return F_correction, S_correction, U_correction, Cv_correction
end

function second_order(c::Cumulants, kB, T, stochastic::Bool)

    pref = stochastic ? -1.0 : 1.0
    β = 1 / (kB*T)

    F_correction = pref * c.κ₂ / (2*kB*T)
    S_correction =  pref * (c.κ₂ - T*c.∂κ₂_∂T) / (2*kB*T*T)
    U_correction =  pref * (c.κ₂ - 0.5*T*c.∂κ₂_∂T) / (kB*T)
    Cv_correction =  pref * ((-U_correction*β/T) + β*(0.5*c.∂κ₂_∂T - 0.5*T*c.∂²κ₂_∂T²))

    return F_correction, S_correction, U_correction, Cv_correction
end

function third_order(c::Cumulants, kB, T)

    β = 1 / (kB*T)
    β² = β^2; β³ = β^3

    F_correction = c.κ₃ * β² / 6
    S_correction = (c.κ₃*kB*β³/3) - (β²*c.∂κ₃_∂T/6)
    U_correction = T*((0.5*kB*β³*c.κ₃) - (β²*c.∂κ₃_∂T/6))
    Cv_correction = (U_correction/T) - (3*β²*c.κ₃/2) + (5*β²*c.∂κ₃_∂T/6) + (T*β²*c.∂²κ₃_∂T²/6)

    return F_correction, S_correction, U_correction, Cv_correction
end

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

"""
estimate_thermo_properties(lammps_dump_path::String, stat_file_path::String,
                            kB, temperature::T; limit::Limit = Quantum(), order::Int = 2, 
                            dump_fields::Dict{String, String} = header_dict)

Estimates thermodynamic properties via a cumulant expansion from a LAMMPS dump file and a TDEP infile.stat.

Parameters:
-----------
- `lammps_dump_path::String`: Path to the LAMMPS dump file.
- `stat_file_path::String`: Path to the TDEP infile.stat file.
- `kB`: Boltzmann constant.
- `temperature::T`: Temperature.
- `limit::Limit` = Classical(): Either `Quantum()` or `Classical()`.
- `order::Int` = 2: Order of the cumulant expansion (1 or 2).
- `dump_x_unrolled_names = ["xu", "yu", "zu"]` = header_dict: Name of dump columns corresponding to unrolled displacements. In x,y,z order!!

Returns:
-----------
- A tuple containing the free energy, entropy, internal energy, and heat capacity corrections.
"""
function estimate_thermo_properties(
    lammps_dump_path::String,
    lammps_eq_dump_path::String,
    stat_file_path::String,
    ifc2::AbstractMatrix,
    ω, kB, ħ, temperature; 
    limit::Limit = Classical(), order::Int = 3, 
    dump_x_unrolled_names::AbstractVector{String} = ["xu", "yu", "zu"],
    stochastic::Bool = false
)

    D = length(dump_x_unrolled_names)

    ld = thermo_prop_checks(lammps_dump_path, order, dump_x_unrolled_names)

    initial_positions, atom_masses, E_total, V =
        parse_files(lammps_eq_dump_path, stat_file_path, temperature)

    @info "Parsing Displacements"
    u = load_displacements(ld, initial_positions, 
                            dump_x_unrolled_names = dump_x_unrolled_names, D = D)
    V₂ = V_harmonic.(Ref(ifc2), eachcol(u))
    ΔV = V .- V₂

    F₀ = F_harmonic(ω, ħ, kB, temperature, limit)
    S₀ = S_harmonic(ω, ħ, kB, temperature, limit)
    U₀ = U_harmonic(ω, ħ, kB, temperature, limit)
    Cᵥ₀ = Cᵥ_harmonic(ω, kB, temperature, limit)

    #*TODO ADD BOOTSRAP ERROR ESTIMATES
    c = Cumulants(V, ΔV, kB, temperature)
    ΔF₁, ΔS₁, ΔU₁, ΔCᵥ₁ = first_order(c, temperature) 
    ΔF₂, ΔS₂, ΔU₂, ΔCᵥ₂ = 0.0, 0.0, 0.0, 0.0
    ΔF₃, ΔS₃, ΔU₃, ΔCᵥ₃ = 0.0, 0.0, 0.0, 0.0

    if order >= 2
        ΔF₂, ΔS₂, ΔU₂, ΔCᵥ₂ = second_order(c, kB, temperature, stochastic)
    end

    if order >= 3
        ΔF₃, ΔS₃, ΔU₃, ΔCᵥ₃ = third_order(c, kB, temperature)
    end

    df = DataFrame(
        F = [F₀, ΔF₁, ΔF₂, ΔF₃],
        S = [S₀, ΔS₁, ΔS₂, ΔS₃],
        U = [U₀, ΔU₁, ΔU₂, ΔU₃],
        Cv = [Cᵥ₀, ΔCᵥ₁, ΔCᵥ₂, ΔCᵥ₃]
    )

    # Estimate true internal energy and heat capacity
    U_MD = mean(E_total)
    Cᵥ_MD = var(E_total) / (kB * temperature^2)

    return df, U_MD, Cᵥ_MD
end