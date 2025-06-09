export estimate_thermo_properties, CumulantData

# Various derivatives of ⟨O⟩
∂A_∂T(A, V, kB, T) = cov(A, V) / (kB * T * T)
∂AB_∂T(A, B, V, kB, T, dA = ∂A_∂T(A, V, kB, T), dB =  ∂A_∂T(B, V, kB, T)) = (mean(A) * dB) + (mean(B) * dA)
∂²A_∂T²(A, V, kB, T, dA = ∂A_∂T(A, V, kB, T)) = (-2*dA/T) + ((1/(kB*T*T)) * (∂A_∂T(A.*V, V, kB, T) - ∂AB_∂T(A, V, V, kB, T, dA)))

function skew(X)
    return mean((X .- mean(X)).^3)
end 

struct CumulantData{O,A,B,C}
    κ::Measurement{A}
    ∂κ_∂T::Measurement{B}
    ∂²κ_∂T²::Measurement{C}
end

order(::CumulantData{O}) where O = O

function CumulantData(V, ΔV, kB, T, ::Val{1})

    ba_ΔV_V = BlockAveragable(ΔV, V)

    @sync begin
        Threads.@spawn κ₁ = @ba mean(ΔV)
        Threads.@spawn ∂κ₁_∂T = @ba nargs=2 ∂A_∂T(ΔV, V, kB, T)
        Threads.@spawn ∂²κ₁_∂T² = @ba nargs=2 ∂²A_∂T²(ΔV, V, kB, T)
    end

    return CumulantData{1}(κ₁, ∂κ₁_∂T, ∂²κ₁_∂T²)
end

function CumulantData(V, ΔV, kB, T, ::Val{2})

    ΔV² = ΔV .^ 2

    @sync begin
        Threads.@spawn κ₂ = @ba var(ΔV)
        Threads.@spawn ∂ΔV²_∂T = @ba nargs=2 ∂A_∂T(ΔV², V, kB, T)
        Threads.@spawn ∂²ΔV²_∂T² = @ba nargs=2 ∂²A_∂T²(ΔV², V, kB, T)
    end

    ∂κ₂_∂T = ∂ΔV²_∂T - (2*κ₁*∂κ₁_∂T)
    ∂²κ₂_∂T² = ∂²ΔV²_∂T² - 2*((∂κ₁_∂T^2) + (κ₁*∂²κ₁_∂T²))

    return CumulantData{2}(κ₂, ∂κ₂_∂T, ∂²κ₂_∂T²)
end

function CumulantData(V, ΔV, kB, T, ::Val{3})

    ΔV³ = ΔV² .* ΔV

    @sync begin
        Threads.@spawn κ₃ = @ba skew(ΔV)
        Threads.@spawn ∂ΔV³_∂T = @ba nargs=2 ∂A_∂T(ΔV³, V, kB, T)
        Threads.@spawn ∂²ΔV³_∂T² = @ba nargs=2 ∂²A_∂T²(ΔV³, V, kB, T)
        Threads.@spawn μ_ΔV² = @ba nargs=2 mean(ΔV², V, kB, T)
    end

    ∂κ₃_∂T = ∂ΔV³_∂T - 3*κ₁*∂ΔV²_∂T + 3*μ_ΔV²*∂κ₁_∂T
    ∂²κ₃_∂T² = ∂²ΔV³_∂T² - 3*(∂κ₁_∂T*∂ΔV²_∂T + κ₁*∂²ΔV²_∂T²) + 3*(∂ΔV²_∂T*∂κ₁_∂T + μ_ΔV²*∂²κ₁_∂T²)

    return CumulantData{3}(κ₃, ∂κ₃_∂T, ∂²κ₃_∂T²)
end

function first_order_corrections(c1::CumulantData{1}, T)
   
    F_correction = c1.κ
    S_correction = -c1.∂κ_∂T
    U_correction = c1.κ - T*c.∂κ_∂T
    Cv_correction = -T*c1.∂²κ_∂T²

    return F_correction, S_correction, U_correction, Cv_correction
end

function second_order_corrections(c2::CumulantData{2}, kB, T, stochastic::Bool)

    pref = stochastic ? -1.0 : 1.0
    β = 1 / (kB*T)

    F_correction = pref * c2.κ / (2*kB*T)
    S_correction =  pref * (c2.κ - T*c2.∂κ_∂T) / (2*kB*T*T)
    U_correction =  pref * (c2.κ - 0.5*T*c2.∂κ_∂T) / (kB*T)
    Cv_correction =  pref * ((-U_correction*β/T) + β*(0.5*c2.∂κ_∂T - 0.5*T*c2.∂²κ_∂T²))

    return F_correction, S_correction, U_correction, Cv_correction
end

function third_order_corrections(c3::CumulantData{3}, kB, T)

    β = 1 / (kB*T)
    β² = β^2; β³ = β^3

    F_correction = c3.κ * β² / 6
    S_correction = (c3.κ*kB*β³/3) - (β²*c3.∂κ_∂T/6)
    U_correction = T*((0.5*kB*β³*c3.κ) - (β²*c3.∂κ_∂T/6))
    Cv_correction = (U_correction/T) - (3*β²*c3.κ/2) + (5*β²*c3.∂κ_∂T/6) + (T*β²*c3.∂²κ_∂T²/6)

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

    u = load_displacements(ld, initial_positions, 
                            dump_x_unrolled_names = dump_x_unrolled_names, D = D)
    V₂ = V_harmonic.(Ref(ifc2), eachcol(u))
    ΔV = V .- V₂

    F₀ = F_harmonic(ω, ħ, kB, temperature, limit)
    S₀ = S_harmonic(ω, ħ, kB, temperature, limit)
    U₀ = U_harmonic(ω, ħ, kB, temperature, limit)
    Cᵥ₀ = Cᵥ_harmonic(ω, kB, temperature, limit)

    c1 = Cumulants(V, ΔV, kB, temperature, Val{1}())
    ΔF₁, ΔS₁, ΔU₁, ΔCᵥ₁ = first_order(c1, temperature) 

    ΔF₂, ΔS₂, ΔU₂, ΔCᵥ₂ = 0.0, 0.0, 0.0, 0.0
    if order >= 2
        c2 = Cumulants(V, ΔV, kB, temperature, Val{2}())
        ΔF₂, ΔS₂, ΔU₂, ΔCᵥ₂ = second_order(c2, kB, temperature, stochastic)
    end

    ΔF₃, ΔS₃, ΔU₃, ΔCᵥ₃ = 0.0, 0.0, 0.0, 0.0
    if order >= 3
        c3 = Cumulants(V, ΔV, kB, temperature, Val{3}())
        ΔF₃, ΔS₃, ΔU₃, ΔCᵥ₃ = third_order(c3, kB, temperature)
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