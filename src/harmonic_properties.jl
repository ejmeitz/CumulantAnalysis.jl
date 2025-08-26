
function harmonic_properties(
    estim::ThermoEstimator,
    ω::AbstractVector,
    _kB, _ħ;
    normalization_factor = 1.0
)
    F₀ = F_harmonic(ω, _ħ, _kB, ustrip(estim.temperature), limit(estim)) ./ normalization_factor
    S₀ = S_harmonic(ω, _ħ, _kB, ustrip(estim.temperature), limit(estim)) ./ normalization_factor
    U₀ = U_harmonic(ω, _ħ, _kB, ustrip(estim.temperature), limit(estim)) ./ normalization_factor
    Cᵥ₀ = Cᵥ_harmonic(ω, _ħ, _kB, ustrip(estim.temperature), limit(estim)) ./ normalization_factor

    return F₀, S₀, U₀, Cᵥ₀
end

# Returns in eV/atom, eV/K/atom, and eV/K/atom
# function read_tdep_thermo_props(ifc_dir)
#     f = open(joinpath(ifc_dir, "outfile.free_energy"), "r")
#     data = readlines(f)
#     close(f)
#     _, F₀, S₀, Cᵥ₀ = parse.(Float64, strip(split(data)))
#     return F₀, S₀, Cᵥ₀
# end

function harmonic_properties(estim::ThermoEstimator, ifc_dir::String)

    pdr = PhononDispersionRelations(; dumpgrid = true, temperature = Float64(ustrip(estim.temperature)))
    isdir(ifc_dir) || error(ArgumentError("ifc_dir passed to harmonic_properties is not a valid directory: $(ifc_dir)"))
    
    cd(ifc_dir) do
        execute(pdr, ifc_dir, Threads.nthreads(), false)
    end

    p = joinpath(ifc_dir, "outfile.grid_dispersions.hdf5")
    isfile(p) || error("Could not find grid dispersion file to calculate harmonic properties")
    freqs_rad_s = 2pi .* h5read(p, "frequencies") ./ 1e12

    N_branch, N_full_q_point = size(freqs_rad_s)

    # Converts it to per-atom
    # N_branch / 3 == N_atoms_per_unitcell
    N = N_full_q_point / (N_branch / 3)

    return harmonic_properties(estim, reduce(vcat, freqs_rad_s), ustrip(kB), ustrip(ħ); normalization_factor = N)
end

#* FIX CONTRIBUTION FROM Zero-Point MOTION ON RIGID TRANSLATION MODES??
#* SEE HOW TDEP IMPLEMENTS THINGS LIKE FREE ENERGY

function V_harmonic(ifc2::AbstractMatrix, u::AbstractVector)
    return 0.5 * ((transpose(u) * ifc2) * u)
end

function sum_over_freqs(freqs, f::Function)
    res = 0.0
    for freq in freqs
        if freq > FREQ_TOL
            res += f(freq)
        end
    end
    return res
end

function U_harmonic(ω, ħ, kB, T, ::Type{Quantum})
    f = (freq) -> (ħ*freq) * ((1 / (exp(ħ*freq/(kB*T)) - 1)) + 0.5)
    return sum_over_freqs(ω, f)
end

function U_harmonic(ω, ħ, kB, T, ::Type{Classical})
    n_nonzero = count(freq -> freq > FREQ_TOL, ω)
    return n_nonzero*kB*T
end

function F_harmonic(ω, ħ, kB, T, ::Type{Quantum})
    kBT = kB * T
    f = (freq) -> (0.5*ħ*freq) + kBT * log(1 - exp(-ħ*freq/kBT))
    return sum_over_freqs(ω, f)
end

function F_harmonic(ω, ħ, kB, T, ::Type{Classical})
    kBT = kB * T
    f = (freq) -> log(ħ*freq/kBT)
    return kBT * sum_over_freqs(ω, f)
end

function S_harmonic(ω, ħ, kB, T, ::Type{Quantum})
    kBT = kB * T
    f = (freq) -> ((ħ*freq/kBT) / (exp(ħ*freq/kBT) - 1)) - log(1 - exp(-ħ*freq/kBT))
    return kB * sum_over_freqs(ω, f)
end

function S_harmonic(ω, ħ, kB, T, ::Type{Classical})
    f = (freq) -> (1 - log(ħ*freq/(kB * T)))
    return kB * sum_over_freqs(ω, f)
end

function Cᵥ_harmonic(ω, ħ, kB, T, ::Type{Quantum})
    tkBT =  2 * kB * T
    f = (freq) -> ((ħ*freq/tkBT)^2) * (csch(ħ*freq/tkBT)^2)
    return kB * sum_over_freqs(ω, f)
end

function Cᵥ_harmonic(ω, ħ, kB, T, ::Type{Classical})
    n_nonzero = count(freq -> freq > FREQ_TOL, ω)
    return n_nonzero*kB
end