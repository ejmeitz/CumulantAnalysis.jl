export estimate_thermo_properties, Quantum, Classical, Cumulants

struct Cumulants{A,B,C,D,E,F}
    κ₁::A
    ∂κ₁_∂T::B
    ∂²κ₁_∂T²::C
    κ₂::D
    ∂κ₂_∂T::E
    ∂²κ₂_∂T²::F
end

function Cumulants(V, ΔV, kB, T)

    ΔV² = ΔV .^ 2
      
    # ∂⟨A⟩/∂T = cov(A, V)/(kB*T*T)
    ∂A_∂T(A) = cov(A, V) / (kB * T * T)
    # ∂(⟨A⟩⟨B⟩)/∂T = <A>∂<B>∂T + <B>∂<A>∂T
    ∂AB_∂T(A, B; dA = ∂A_∂T(A), dB =  ∂A_∂T(B)) = (mean(A) * dB) + (mean(B) * dA)
    # ∂²<A>/∂T²
    ∂²A_∂T²(A; dA = ∂A_∂T(A)) = (-2*dA/(kB*T*T*T)) + ((1/(kB*T*T)) * (∂A_∂T(A.*B) - ∂AB_∂T(A, V; dA = dA)))

    κ₁ = mean(ΔV)
    ∂κ₁_∂T = ∂A_∂T(ΔV)
    ∂²κ₁_∂T² = ∂²A_∂T²(ΔV; dA = ∂κ₁_∂T)

    κ₂ = var(ΔV)
    ∂ΔV²_∂T = ∂A_∂T(ΔV²)
    ∂κ₂_∂T = ∂ΔV²_∂T - (2*κ₁*∂κ₁_∂T)
    ∂²κ₂_∂T² = ∂²A_∂T²(ΔV²; dA = ∂ΔV²_∂T) - 2*((∂κ₁_∂T^2) + (κ₁*∂²κ₁_∂T²))

    return Cumulants(κ₁, ∂κ₁_∂T, ∂²κ₁_∂T², κ₂, ∂κ₂_∂T, ∂²κ₂_∂T²)
end

function first_order(cd::CumulantData, T)
   
    F_correction = cd.κ₁
    S_correction = -cd.∂κ₁_∂T
    U_correction = cd.κ₁ - T*cd.∂κ₁_∂T
    Cv_correction = -T*cd.∂²κ₁_∂T²

    return F_correction, S_correction, U_correction, Cv_correction
end

function second_order(cd::CumulantData, kB, T, stochastic::Bool)

    pref = stochastic ? -1.0 : 1.0

    F_correction = pref * cd.κ₂ / (2*kB*T)
    S_correction =  pref * (cd.κ₂ - T*cd.∂κ₂_∂T) / (2*kB*T*T)
    U_correction =  pref * (cd.κ₂ - 0.5*T*cd.∂κ₂_∂T) / (kB*T)
    Cv_correction =  (pref / kB) * ((cd.∂κ₂_∂T /T) - (cd.κ₂/(T*T)) - (0.5 * cd.∂²κ₂_∂T²))

    return F_correction, S_correction, U_correction, Cv_correction
end

# name in lammps dump ==> name used in program
const header_dict = Dict(
    "DisplacementX" => "ux",
    "DisplacementY" => "uy",
    "DisplacementZ" => "uz",
)

function thermo_prop_checks(lammps_dump_path, order, dump_fields)
    if order ∉ [1,2]
        @error "Can only calculate first and second order cumulant corrections"
    end

    @info "Parsing LAMMPS dump file at $(lammps_dump_path)"
    ld = LammpsDump(lammps_dump_path)
    @info "Found $(n_samples(ld)) samples and $(n_atoms(ld))"

    required_fields = keys(dump_fields)
    actual_fields = fields(ld)
    @assert issubset(required_fields, actual_fields) "Dump file needs $(required_fields) fields, got $(actual_fields). Can re-name with dump_fields kwarg"

    return ld
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
- `limit::Limit` = Quantum(): Either `Quantum()` or `Classical()`.
- `order::Int` = 2: Order of the cumulant expansion (1 or 2).
- `dump_disp_names = ["ux", "uy", "uz"]` = header_dict: Name of dump columns corresponding to displacements. In x,y,z order!!

Returns:
-----------
- A tuple containing the free energy, entropy, internal energy, and heat capacity corrections.
"""
function estimate_thermo_properties(lammps_dump_path::String, stat_file_path::String, ifc2::Matrix{T},
                                     ω, kB, ħ, temperature::T; limit::Limit = Quantum(), order::Int = 2, 
                                     dump_disp_names::AbstractVector{String} = ["ux", "uy", "uz"], stochastic::Bool = false)

    D = length(dump_fields)

    ld = thermo_prop_checks(lammps_dump_path, order, dump_fields)
    bulk_properties = readdlm(stat_file_path, comments = true)

    E_total = bulk_properties[:, 3]
    V = bulk_properties[:, 4]
    T_sim = bulk_properties[:, 6]

    if mean(T_sim) - temperature > 1e-1
        @warn "Simulation temperature $(mean(T_sim)) different from target temperature $(temperature). Results may be inaccurate."
    end

    u = load_displacements(ld, dump_disp_names; D = D)
    V₂ = V_harmonic.(Ref(ifc2), eachcol(u))
    ΔV = V .- V₂

    F₀ = F_harmonic(ω, ħ, kB, temperature, limit)
    S₀ = S_harmonic(ω, ħ, kB, temperature, limit)
    U₀ = U_harmonic(ω, ħ, kB, temperature, limit)
    Cᵥ₀ = Cᵥ_harmonic(ω, kB, temperature, limit)

    cd = Cumulants(V, ΔV, kB, temperature)
    ΔF₁, ΔS₁, ΔU₁, ΔCᵥ₁ = first_order(cd, temperature) 
    ΔF₂, ΔS₂, ΔU₂, ΔCᵥ₂ = 0.0, 0.0, 0.0, 0.0

    if order >= 2
        ΔF₂, ΔS₂, ΔU₂, ΔCᵥ₂ = second_order(cd, kB, temperature, stochastic)
    end

    df = DataFrame(
        F = [F₀, ΔF₁, ΔF₂],
        S = [S₀, ΔS₁, ΔS₂],
        U = [U₀, ΔU₁, ΔU₂],
        Cv = [Cᵥ₀, ΔCᵥ₁, ΔCᵥ₂]
    )

    # Estimate true internal energy and heat capacity
    U_MD = mean(E_total)
    Cᵥ_MD = var(E_total) / (kB * temperature^2)

    return df, U_MD, Cᵥ_MD
end