export estimate, save

function calculate_cumulants(V, V₂, V₃, V₄, T, n_atoms, ce::CumulantEstimator{O}) where O

    ΔF = zeros(O+1); ΔS = zeros(O+1)
    ΔU = zeros(O+1); ΔCᵥ = zeros(O+1)

    ΔF[1], ΔS[1], ΔU[1], ΔCᵥ[1] = constant_corrections(ce, V, V₂, V₃, V₄, T, n_atoms)

    c1 = CumulantData(V, V₂, V₃, V₄, T, n_atoms, Val{1}(), ce)
    ΔF[2], ΔS[2], ΔU[2], ΔCᵥ[2] = first_order_corrections(c1, T) 

    if O >= 2
        c2 = CumulantData(V, V₂, V₃, V₄, T, n_atoms, c1, Val{2}(), ce)
        ΔF[3], ΔS[3], ΔU[3], ΔCᵥ[3] = second_order_corrections(c2, T, true)
    end

    # if O >= 3
    #     c3 = CumulantData(V, V₂, V₃, V₄, T, n_atoms, c1, Val{3}(), ce)
    #     ΔF[4], ΔS[4], ΔU[4], ΔCᵥ[4] = third_order_corrections(c3, T)
    # end

    return ΔF, ΔS, ΔU, ΔCᵥ

end

function bootstrap_corrections(V, V₂, V₃, V₄, T, outpath,
                                ce::CumulantEstimator{O, L}, Nat::Int) where {O, L <: Limit}

    # these are returned per-atom
    @info "Calculating Harmonic Properties"
    F₀, S₀, U₀, Cᵥ₀ = harmonic_properties(T, L, outpath)
    @info "Calculated Harmonic Properties"

    # Get point estimate of corrections
    ΔF, ΔS, ΔU, ΔCᵥ = calculate_cumulants(V, V₂, V₃, V₄, T, Nat, ce)
    F_total_point = sum(ΔF) + (F₀*Nat)
    S_total_point = sum(ΔS) + (S₀*Nat)
    U_total_point = sum(ΔU) + (U₀*Nat)
    Cᵥ_total_point = sum(ΔCᵥ) + (Cᵥ₀*Nat)

    # Estimate standard error by bootstrapping
    is = zeros(Int, length(V))
    ΔFs = zeros(O+1, ce.n_boot); ΔSs = zeros(O+1, ce.n_boot)
    ΔUs = zeros(O+1, ce.n_boot); ΔCᵥs = zeros(O+1, ce.n_boot)

    p = Progress(ce.n_boot, "Bootstrapping Corrections")
    for i in 1:ce.n_boot
        sample!(1:length(V), is; replace = true)
        ΔFs[:,i], ΔSs[:,i], ΔUs[:,i], ΔCᵥs[:,i] = calculate_cumulants(V[is], V₂[is], V₃[is], V₄[is], T, Nat, ce)
        next!(p)
    end
    finish!(p)

    F_totals = sum(ΔFs, dims = 1) .+ (F₀*Nat)
    S_totals = sum(ΔSs, dims = 1) .+ (S₀*Nat)
    U_totals = sum(ΔUs, dims = 1) .+ (U₀*Nat)
    Cᵥ_totals = sum(ΔCᵥs, dims = 1) .+ (Cᵥ₀*Nat)
    
    F_SEs = std(ΔFs, dims = 2); S_SEs = std(ΔSs, dims = 2)
    U_SEs = std(ΔUs, dims = 2); Cᵥ_SEs = std(ΔCᵥs, dims = 2)

    kBNat = ustrip(CumulantAnalysis.kB) * Nat

    F = BootstrapCumualantEstimate(
        F₀, SVector(ΔF...) ./ Nat, SVector(F_SEs...) ./ Nat,
        F_total_point / Nat, std(F_totals) / Nat, "F", "[eV/atom]"
    )

    S = BootstrapCumualantEstimate(
        S₀ / kB, SVector(ΔS...) ./ kBNat, SVector(S_SEs...) ./ kBNat,
        S_total_point / kBNat, std(S_totals) / kBNat, "S", "[kB / atom]"
    )

    U = BootstrapCumualantEstimate(
        U₀, SVector(ΔU...) ./ Nat, SVector(U_SEs...) ./ Nat,
        U_total_point / Nat, std(U_totals) / Nat, "U", "[eV/atom]"
    )

    Cᵥ = BootstrapCumualantEstimate(
        Cᵥ₀ / kB, SVector(ΔCᵥ...) ./ kBNat, SVector(Cᵥ_SEs...) ./ kBNat,
        Cᵥ_total_point / kBNat, std(Cᵥ_totals) / kBNat, "Cv", "[kB / atom]"
    )

    return F, S, U, Cᵥ

end

# Potentially useful for gauging convergence of different approaches
# Bootstrap estimates error on kappa and its derivatives for all orders
function do_size_study(ce::CumulantEstimator{O}, outpath, V, V₂, V₃, V₄, T, n_atoms) where O

    β = 1 / (kB*T)
    min_samples = (length(V) < 500) ? 10 : 100

    lg_pts = range(log10(min_samples), log10(length(V)), length = 12)
    Ns = round.(Int, 10 .^ lg_pts)

    κs = zeros(length(Ns), O+1, ce.n_boot)
    ∂κs = zeros(length(Ns), O+1, ce.n_boot)
    ∂²κs = zeros(length(Ns), O+1, ce.n_boot)

    p = Progress(length(Ns) * ce.n_boot, "Sampling Study")

    κ_point = zeros(length(Ns))
    ∂κ_point = zeros(length(Ns))
    ∂²κ_point = zeros(length(Ns))
    
    for (i,N) in enumerate(Ns)
        # pre-allocate things
        idxs = zeros(Int, N)

        # Samples are IID so our "fake" smaller dataset we'll
        # just take as the first N samples 
        V_sub  = @views V[1:N]
        V₂_sub = @views V₂[1:N]
        V₃_sub = @views V₃[1:N]
        V₄_sub = @views V₄[1:N]

        #!DO POINT ESTIMATE

        for j in 1:ce.n_boot

            sample!(1:N, idxs; replace = true)

            V_samples = V_sub[idxs]
            V₂_samples = V₂_sub[idxs]
            V₃_samples = V₃_sub[idxs]
            V₄_samples = V₄_sub[idxs]

            c1 = CumulantData(V_samples, V₂_samples, V₃_samples, V₄_samples, T, n_atoms, Val{1}(), ce)
            c2 = CumulantData(V_samples, V₂_samples, V₃_samples, V₄_samples, T, n_atoms, c1, Val{2}(), ce)
            # c3 = CumulantData(V_samples, V₂_samples, V₃_samples, V₄_samples, T, n_atoms, c1, Val{3}(), ce)

            cds = (c1, c2)#, c3)

            for co in 0:O         
                if co == 0
                    if ce isa EffectiveHamiltonianEstimator
                        κs[i, co + 1, j] = NaN
                        ∂κs[i, co + 1, j] = NaN
                        ∂²κs[i, co + 1, j] = NaN
                    else
                        #! UPDATE THIS TO USE CV SO IT MATCHES WHAT CODE ACTUALLY DOES
                        X = V₀_rv(ce, V_samples, V₂_samples, V₃_samples, V₄_samples)
                        κs[i, co + 1, j] = get_V₀(ce, V, V₂, V₃, V₄) #* DO NOT SET TO mean(X), BREAKS MixedEstimator
                        ∂κs[i, co + 1, j] = ∂A_∂T(X, V₂_samples, T)
                        ∂²κs[i, co + 1, j] = ∂²A_∂T²(X, V₂_samples, T) 
                    end
                else
                    κs[i, co + 1, j] = cds[co].κ
                    ∂κs[i, co + 1, j] = cds[co].∂κ_∂T
                    ∂²κs[i, co + 1, j] = cds[co].∂²κ_∂T²
                end
            end
            next!(p)
        end
    end
    finish!(p)

    #! TODO NON-DIMENSONALIZE
    #! TODO REPLACE WITH POINT ESTIMATE
    κ_estimates = mean(κs; dims = 3) ./ n_atoms
    ∂κ_estimates = mean(∂κs; dims = 3) ./ n_atoms
    ∂²κ_estimates = mean(∂²κs; dims = 3) ./ n_atoms

    κ_SEs = std(κs; dims = 3) ./ n_atoms
    ∂κ_SEs = std(∂κs; dims = 3) ./ n_atoms
    ∂²κ_SEs = std(∂²κs; dims = 3) ./ n_atoms

    data_fmt_str = (N) -> Printf.Format("%7d"*join(fill("%15.8f", N), " "))
    d_fmt = data_fmt_str(6)
    str_fmt_str = (N) -> Printf.Format("%7s"*join(fill("%15s", N-1), " "))

    for co in 0:O
        header = ["N" "k" "k_SE" "dk_dT" "dk_dT_SE" "d2k_dT2" "d2k_dT2_SE"]

        open(joinpath(outpath, "outfile.nsamples_study_order$(co)"), "w") do f
            println(f, "# Standard Error estimated from $(ce.n_boot) bootstraps of size N from origianl dataset with $(length(V)) samples")
            println(f, "# Temperature $(T), N_atoms $(n_atoms)")
            println(f, Printf.format(str_fmt_str(length(header)), header...))
            for i in eachindex(Ns)
                println(f, Printf.format(d_fmt, Ns[i], κ_estimates[i, co+1], κ_SEs[i, co+1],
                                                       ∂κ_estimates[i, co+1], ∂κ_SEs[i, co+1], 
                                                       ∂²κ_estimates[i, co+1], ∂²κ_SEs[i, co+1]))
            end
        end
    end

end

function parse_energies(path, is_hdf5::Val{false})

    data = readdlm(path, comments = true)

    # Parse n_atoms and T from header
    f = open(path, "r")
    readline(f) # skip
    T = parse(Float64, split(strip(readline(f)))[end])
    n_atoms = parse(Int, split(strip(readline(f)))[end])
    close(f)

    # energies in file are meV / atom, Converts to eV
    conv = n_atoms / 1000
    
    @views Vₚ = data[:, 3] .* conv
    @views V₂ = data[:, 4] .* conv
    @views V₃ = data[:, 5] .* conv
    @views V₄ = data[:, 6] .* conv

    return T, n_atoms, Vₚ, V₂ , V₃, V₄
end

function parse_energies(path, is_hdf5::Val{true})

    f = h5open(path, "r")

    Vₚ = read(f, "polar_potential_energy")
    V₂ = read(f, "secondorder_potential_energy")
    V₃ = read(f, "thirdorder_potential_energy")
    V₄ = read(f, "fourthorder_potential_energy")

    T = Float64(attrs(f)["temperature_thermostat"][1])
    n_atoms = Int64(attrs(f)["number_of_atoms"][1])

    close(f)
    
    return T, n_atoms, Vₚ, V₂, V₃, V₄
end

function calculate_true_energies(calc, nconf, ssposcar_path, outpath)
    isnothing(calc) && throw(ArgumentError("This CumulantEstimator requires a calculator for energies, but `calc` was nothing."))
    (AtomsCalculators.energy_unit(calc) != u"eV") && throw(ArgumentError("Expected calculator with eV energy units. Got $(AtomsCalculators.energy_unit(calc)) "))

    sys_ss = TDEPSystem(ssposcar_path)
    V = zeros(nconf)

    p = Progress(nconf, desc = "Calculating True Energies")
    L = typeof(1.0u"Å")

    configs_path = joinpath(outpath, "outfile.canonical_configs.hdf5")
    configs = h5read(configs_path, "positions")

    for i in 1:nconf   
        @views c = configs[:,:,i]
        #! CAN I REMOVE THIS ALLOCATION?
        new_sys = TDEPSystem(sys_ss, vec(collect(reinterpret(SVector{3, L}, c))))

        V[i] = ustrip(AtomsCalculators.potential_energy(new_sys, calc))

        next!(p)
    end
    finish!(p)
    
    open(joinpath(outpath, "outfile.true_potential_energy"), "w") do f
        writedlm(f, V)
    end

    return V
end

function estimate(
        ce::CumulantEstimator{O,L},
        T::Real, # kelvin
        outpath::String;
        ucposcar_path::String = joinpath(outpath, "infile.ucposcar"),
        ssposcar_path::String = joinpath(outpath, "infile.ssposcar"),
        nthreads::Int = Threads.nthreads(),
        rm_configs::Bool = true
    ) where {O, L <: Limit}

    @assert O <= 2 "Up to second order cumulant corrections are supported. Asked for $(O)."

    isfile(ucposcar_path) || throw(ArgumentError("ucposcar path is not a file: $(ucposcar_path)"))
    isfile(ssposcar_path) || throw(ArgumentError("ssposcar path is not a file: $(ssposcar_path)"))

    new_uc_path = joinpath(outpath, "infile.ucposcar")
    new_ss_path = joinpath(outpath, "infile.ssposcar")
    isfile(new_uc_path) || cp(ucposcar_path, new_uc_path; force = true)
    isfile(new_ss_path) || cp(ssposcar_path, new_ss_path; force = true)

    move_ifcs(ce, outpath)

    # Check for mpirun
    res = run(`which mpirun`; wait = false)
    found_mpirun = success(res)
    found_mpirun || @warn "Could not find mpirun on path, defaulting to 1 thread."

    flags = ""
    flags *= (ce isa HarmonicEstimator) ? "" : "--thirdorder --fourthorder "
    flags *= (L === Quantum) ? "--quantum " : ""
    flags *= needs_true_V(ce) ? "--dumpconfigs " : ""

    if found_mpirun
        cmd_str = `mpirun -np $(nthreads) effective_hamiltonian --nconf $(ce.nconf) --temperature $(T) $(split(flags))`
    else
        cmd_str = `effective_hamiltonian --nconf $(ce.nconf) --temperature $(T) $(split(flags))`
    end
    @info cmd_str
    cd(outpath) do 
        run(cmd_str)
    end

    # If configs are dumped to calculate V, use HDF5
    # this gurantees order of positions and energies matches.
    if needs_true_V(ce)
        tep_energies_path = joinpath(outpath, "outfile.canonical_configs.hdf5")
        is_hdf5 = Val{true}()
    else
        tep_energies_path = joinpath(outpath, "outfile.energies")
        is_hdf5 = Val{false}()
    end

    @info "Parsing Energies"
    T_file, n_atoms, Vₚ, V₂, V₃, V₄ = parse_energies(tep_energies_path, is_hdf5)

    if T_file != T
        @warn "You said tempearture was $T, but parsed temperature as $(T_file) from outfile.energies"
    end

    if !(sum(Vₚ) ≈ 0.0)
        error("Got non-zero polar term. I don't know what to do with that.")
    end

    if needs_true_V(ce)
        V = calculate_true_energies(ce.force_calculator, ce.nconf, ssposcar_path, outpath)
    else
        V = zeros(eltype(V₂), size(V₂))
    end

    if rm_configs
        rm(joinpath(outpath, "outfile.canonical_configs.hdf5"); force = true)
    end
    res = bootstrap_corrections(V, V₂, V₃, V₄, T, outpath, ce, n_atoms)

    save.(res, Ref(outpath), Ref(ce.n_boot))

    # Compute some statistics to assess convergence with N
    do_size_study(ce, outpath, V, V₂, V₃, V₄, T, n_atoms)

    return res

end

function save(bce::BootstrapCumualantEstimate{L}, outdir::String, n_boot) where L
    prop_name = bce.property
    unit_str = bce.unit_str

    outpath_mean = (ext) -> joinpath(outdir, prop_name * "_mean.$(ext)")
    mean_data = OrderedDict(prop_name*"0" => bce.harmonic)
    SE_data = OrderedDict(prop_name*"0_SE" => 0.0)

    for order in 0:(L-1)
        if order == 0
            mean_data[prop_name * "_offset"] = bce.corrections[order+1]
            SE_data[prop_name * "_offset_SE"] = bce.correction_SEs[order+1]
        else
            mean_data[prop_name * "$(order)"] = bce.corrections[order+1]
            SE_data[prop_name * "$(order)_SE"] = bce.correction_SEs[order+1]
        end
    end

    mean_data[prop_name*"_total"] = bce.total
    SE_data[prop_name*"_total_SE"] = bce.total_SE

    float_fmt_str = (N) -> Printf.Format(join(fill("%15.7f", N), " "))
    str_fmt_str = (N) -> Printf.Format(join(fill("%15s", N), " "))

    # Human Readable Version
    header = collect(keys(mean_data))
    mean_values = collect(values(mean_data))
    SE_values = collect(values(SE_data))
    N = length(header)
    open(outpath_mean("txt"), "w") do f
        println(f, "# $(prop_name) Units: $(unit_str)")
        println(f, "# Row 1: Values, Row 2: Standard Error estimated from $(n_boot) bootstraps")
        println(f, Printf.format(str_fmt_str(N), header...))
        println(f, Printf.format(float_fmt_str(N), mean_values...))
        println(f, Printf.format(float_fmt_str(N), SE_values...))
    end
    # Save to HDF5
    h5open(outpath_mean("h5"), "w") do file
        for (k,v) in mean_data
            write(file, k, v)
        end
        for (k,v) in SE_data
            write(file, k, v)
        end
    end
end