"""
    first_order_harmonic_corrections(Ts, uc, ifc2s, mesh, L; n_threads) -> NamedTuple

Compute harmonic thermodynamic properties with temperature derivative corrections.

Takes multiple IFC2 objects fitted at different temperatures and computes:
1. Free energy F at each temperature (standard harmonic formula)
2. Temperature derivatives of frequencies dω/dT via Savitzky-Golay filter
3. Corrected S and U accounting for the temperature dependence of frequencies

The correction to heat capacity Cv requires second derivatives d²ω/dT² and
is not included here.

# Arguments
- `Ts::AbstractVector{<:Real}`: Temperatures in Kelvin (must be same length as ifc2s)
- `ucs::AbstractVector{<:CrystalStructure}`: Unit cells for each temperature
- `ifc2s::AbstractVector{<:IFC2}`: IFC2 objects, one for each temperature
- `mesh`: q-point mesh specification (e.g., (8,8,8))
- `L::Type{<:Limit}`: Either `Quantum` or `Classical`
- `n_threads::Integer`: Number of threads for parallel computation

# Returns
Named tuple with fields:
- `S`: Entropy with dω/dT correction (per atom)
- `U`: Internal energy with dω/dT correction (per atom)

# Notes
- All IFC2s must be built on the same unit cell
- Properties are normalized per atom in the primitive cell
- Temperatures should be evenly spaced for accurate Savitzky-Golay derivatives
- Uses quadratic polynomial fit (order=2) with window size = min(n_T, 9)
- Heat capacity is not returned as it requires d²ω/dT² for proper correction
- Assumes symmetry of all unit cells is the same for each temperature
"""
function first_order_harmonic_corrections(
        Ts::AbstractVector{<:Real},
        ucs::AbstractVector{<:CrystalStructure},
        ifc2s::AbstractVector{<:IFC2},
        mesh,
        ::Type{L};
        window_size = 5,
        order = 3,
        n_threads::Integer = Threads.nthreads()
    ) where {L <: Limit}

    n_T = length(Ts)
    @assert n_T == length(ifc2s) "Number of temperatures ($(n_T)) must match number of IFC2s ($(length(ifc2s)))"
    @assert n_T == length(ucs) "Number of temperatures ($(n_T)) must match number of unit cells ($(length(ucs)))"
    @assert n_T >= window_size "Need at least `window_size` ($(window_size)) temperatures to compute derivatives. Pass as kwarg"
    
    # Validate all IFCs match their corresponding unit cells
    for (i, (uc, ifc2)) in enumerate(zip(ucs, ifc2s))
        @assert length(uc) == ifc2.na "IFC2[$i] built on $(ifc2.na) atoms, but unitcell[$i] has $(length(uc)) atoms"
    end

    T_spacing = temperatures[2:end] .- temperatures[1:end-1]

    if !allequal(T_spacing)
        error("Temperature spacing is not uniform. Cannot use Savitzky-Golay filter to improve heat capacity estimate.")
    end

    # Build IBZ from first unit cell - use fractional q-points for all temperatures
    ibz = SimpleMesh(ucs[1], mesh)
    N = ibz.n_atoms_prim
    n_q = length(ibz.k_ibz)
    nb = 3 * length(ucs[1])

    # Storage for frequencies: freqs_all[i_q, i_T, i_mode]
    freqs_all = zeros(n_q, n_T, nb)
    
    # Compute frequencies at all q-points for all temperatures
    Threads.@threads for i_q in 1:n_q
        q = ibz.k_ibz[i_q]
        for i_T in 1:n_T
            D_q = LatticeDynamicsToolkit.dynmat_q(ifc2s[i_T], ucs[i_T], q)
            freqs_sq, _ = get_modes(D_q, Val{LatticeDynamicsToolkit.is_gamma(q)}())
            freqs_all[i_q, i_T, :] .= LatticeDynamicsToolkit.negsqrt.(freqs_sq)
        end
    end

    # Compute dω/dT for each mode at each q-point via Savitzky-Golay filter
    # dωdT_all[i_q, i_T, i_mode] = dω/dT
    dωdT_all = zeros(n_q, n_T, nb)

    for i_q in 1:n_q
        for i_mode in 1:nb
            # Extract ω(T) for this mode
            ω_vec = freqs_all[i_q, :, i_mode]
            # Savitzky-Golay derivative (w.r.t. index), then scale by 1/ΔT
            sg_deriv = savitzky_golay(ω_vec, window_size, order; deriv=1, rate = 1 / T_spacing[1])
            dωdT_all[i_q, :, i_mode] .= sg_deriv.y ./ ΔT
        end
    end

    # Compute properties at each temperature
    S_corr = zeros(n_T)
    U_corr = zeros(n_T)

    for i_T in 1:n_T
        T = Float64(Ts[i_T])
        kBT = kB_Hartree * T

        corr_sum = 0.0 
        
        for i_q in 1:n_q
            w = ibz.weights[i_q]
            
            for i_mode in 1:nb
                ω = freqs_all[i_q, i_T, i_mode]
                dωdT = dωdT_all[i_q, i_T, i_mode]
                
                if ω > lo_freqtol
                    # Correction term: (n + 1/2)(dω/dT) = (∂F/∂ω)(dω/dT)
                    corr_sum += w * dF_dω_single(ω, kBT, L) * dωdT
                end
            end
        end
        
        S_corr[i_T] = -corr_sum / N
        U_corr[i_T] = -T * corr_sum / N
    end

    return S_corr, U_corr
end

# Derivative of F with respect to ω: ∂F/∂ω
dF_dω_single(ω, kBT, ::Type{Quantum}) = (1 / (exp(ω/kBT) - 1)) + 0.5
dF_dω_single(ω, kBT, ::Type{Classical}) = kBT / ω



function estimate_cv_constant_correction(
    temperatures,
    outpath::Function;
    window_size = 5,
    order = 3
)

    T_spacing = temperatures[2:end] .- temperatures[1:end-1]

    if !allequal(T_spacing)
        @warn "Temperature spacing is not uniform. Cannot use Savitzky-Golay filter to improve heat capacity estimate."
        return nothing
    end

    # Load internal energy data across temperatures
    U_offsets = zeros(length(temperatures))
    for (i,T) in enumerate(temperatures)
        U_path = joinpath(outpath(T), "U_mean.h5")
        h5open(U_path, "r") do f
            U_offsets[i] = read(f, "U_offset")
        end
    end

    # Apply Savitzky-Golay filter to get smoothed derivative
    cv_offset_sg = savitzky_golay(U_offsets, window_size, order;
                                   deriv = 1, rate = 1 / T_spacing[1])

    return cv_offset_sg.y ./ kB_eV # Convert to [kB/atom]
end

function improve_constant_corrections(
        temperatures,
        outpath::Function,
        ifcs::AbstractVector{<:IFC2},
        ucs::AbstractVector{<:CrystalStructure},
        mesh,
        ::Type{L};
        kwargs...
    ) where {L <: Limit}

    S_corr, U_corr = first_order_harmonic_corrections(temperatures, ucs, ifcs, mesh, L; kwargs...)
    cv_offset_sg = estimate_cv_constant_correction(temperatures, outpath; kwargs...)

    return S_corr, U_corr, cv_offset_sg

end