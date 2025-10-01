using CumulantAnalysis
using Statistics
using DelimitedFiles

function parse_energies(path)

    data = readdlm(path, comments = true)

    # Parse n_atoms and T from header
    f = open(path, "r")
    readline(f) # skip
    T = parse(Float64, split(strip(readline(f)))[end])
    n_atoms = parse(Int, split(strip(readline(f)))[end])
    close(f)

    # stored as meV/atom, convert to eV
    conv = n_atoms / 1000 #n_atoms / 1000
    
    @views Vₚ = data[:, 3] .* conv
    @views V₂ = data[:, 4] .* conv
    @views V₃ = data[:, 5] .* conv
    @views V₄ = data[:, 6] .* conv

    return T, n_atoms, Vₚ, V₂ , V₃, V₄
end

function estimate_V0(V0_rv, V2, V3, V4, n_atoms, T)

    f = 3*n_atoms - 3

    μ₂ = 0.5*f*kB*T #/ n_atoms # ⟨V2⟩
    Z = V2 .- μ₂ #mean(V2)
    
    σ₂² = 0.5*f*(kB*T)^2 / n_atoms^2
    C1 = Z.^2 .- σ₂²
    C2 = Z .^ 3
    C3 = Z .^ 4 .- 3*σ₂²*Z.^2

    # use zero mean control variates for now
    return cv_estimate(V0_rv, Z; K=25, verbose=true)

end

# Custom implementation that takes advantage of the
# structure of the covariance to improve accuracy
function estimate_∂V0(V0_rv, V2, V3, V4, n_atoms, T)

    f = 3*n_atoms - 3

    μ₂ = 0.5*f*kB*T #/ n_atoms # ⟨V2⟩
    Z = V2 .- μ₂ #mean(V2)

    Y = V0_rv .* Z    
    
    σ₂² = 0.5*f*(kB*T)^2 / n_atoms^2
    C1 = Z.^2 .- σ₂²
    C2 = Z .^ 3
    C3 = Z .^ 4 .- 3*σ₂²*Z.^2

    # use zero mean control variates for now
    return cv_estimate(Y, C1, C2, C3, Z, V3 .* Z; K=25, verbose=true)

end

function V0_harm(V, V2, V3, V4)
    return V .- V2
end

function V0_quartic(V, V2, V3, V4)
    return V .- V2 .- V4 # technically has V3, but <V3> = 0
end

# This calculates V, V2, V3, V4. Can use this data for all possible V0 definitions
basepath = "Z:/emeitz/Projects/Cumulants/LJ_size_effects_QuarticV0/T80/6UC"
tep_energy_path = joinpath(basepath, "outfile.energies")
real_energy_path = joinpath(basepath, "outfile.true_potential_energy")
V0_func = V0_quartic

T, n_atoms, _, V2 , V3, V4 = parse_energies(tep_energy_path)
kB = 8.617333262145e-5 # eV/K

@assert T == 80.0

∂A_∂T(A, V, T_) = cov(A, V) / (kB * T_ * T_)


V = vec(readdlm(real_energy_path))
# V ./= n_atoms # eV / atom

V0_rv = V0_func(V, V2, V3, V4)
raw_V0 = mean(V0_rv)
raw_∂V0_∂T = ∂A_∂T(V0_rv, V2, T)
raw_F_offset = raw_V0 / n_atoms
raw_S_offset = -raw_∂V0_∂T / (n_atoms * kB)
raw_U_offset = (raw_V0 - T*raw_∂V0_∂T) / n_atoms

∂V0_opt_c, ∂V0_opt_c_se, red = CumulantAnalysis.derivative_with_opt_c(V0_rv, V2, T, n_atoms; K = 2)
println("var red: $(red)")

∂V0_cv = estimate_∂V0(V0_rv, V2, V3, V4, n_atoms, T)
V0_cv = estimate_V0(V0_rv, V2, V3, V4, n_atoms, T)


k = 1 / (n_atoms * kB)
println("Raw ∂V0/∂T estimate: $(raw_∂V0_∂T / k) [eV/atom/K]")
println("Optimal c ∂V0/∂T estimate: $(∂V0_opt_c / k) +/- $(∂V0_opt_c_se / k) [eV/atom//K]")
println("Control variates ∂V0/∂T estimate: $(∂V0_cv.mean_cv / k) [eV/atom/K]")

# F_corr = V₀
# S_corr = -∂V₀
# U_corr = V₀ - T*∂V₀
# Cv_corr = -T * ∂²V₀