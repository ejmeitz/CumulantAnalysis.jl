export estimate, save

function calculate_cumulants(V, V₂, V₃, V₄, T, ce::CumulantEstimator{O}) where O

    ΔF = zeros(O+1); ΔS = zeros(O+1)
    ΔU = zeros(O+1); ΔCᵥ = zeros(O+1)

    ΔF[1], ΔS[1], ΔU[1], ΔCᵥ[1] = constant_corrections(ce,  V, V₂, V₃, V₄, T)

    c1 = CumulantData(V, V₂, V₃, V₄, T, Val{1}(), ce)
    ΔF[2], ΔS[2], ΔU[2], ΔCᵥ[2] = first_order_corrections(c1, T) 

    if O >= 2
        c2 = CumulantData(V, V₂, V₃, V₄, T, c1, Val{2}(), ce)
        ΔF[3], ΔS[3], ΔU[3], ΔCᵥ[3] = second_order_corrections(c2, T, true)
    end

    if O >= 3
        c3 = CumulantData(V, V₂, V₃, V₄, T, c1, Val{3}(), ce)
        ΔF[4], ΔS[4], ΔU[4], ΔCᵥ[4] = third_order_corrections(c3, T)
    end

    return ΔF, ΔS, ΔU, ΔCᵥ

end

function bootstrap_corrections(V, V₂, V₃, V₄, T, outpath,
                                ce::CumulantEstimator{O, L}, Nat::Int) where {O, L <: Limit}

    # these are returned per-atom
    @info "Calculating Harmonic Properties"
    F₀, S₀, U₀, Cᵥ₀ = harmonic_properties(T, L, outpath)
    @info "Calculated Harmonic Properties"

    is = zeros(Int, ce.boot_size)

    ΔFs = zeros(O+1, ce.n_boot); ΔSs = zeros(O+1, ce.n_boot)
    ΔUs = zeros(O+1, ce.n_boot); ΔCᵥs = zeros(O+1, ce.n_boot)

    p = Progress(ce.n_boot, "Bootstrapping Corrections")
    for i in 1:ce.n_boot
        sample!(1:length(V), is; replace = true)
        ΔFs[:,i], ΔSs[:,i], ΔUs[:,i], ΔCᵥs[:,i] = calculate_cumulants(V[is], V₂[is], V₃[is], V₄[is], T, ce)
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
        writedlm(f, [V])
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
        tep_energies_path::Union{String, Nothing} = nothing,
    ) where {O, L <: Limit}

    isfile(ucposcar_path) || throw(ArgumentError("ucposcar path is not a file: $(ucposcar_path)"))
    isfile(ssposcar_path) || throw(ArgumentError("ssposcar path is not a file: $(ssposcar_path)"))

    new_uc_path = joinpath(outpath, "infile.ucposcar")
    new_ss_path = joinpath(outpath, "infile.ssposcar")
    isfile(new_uc_path) || cp(ucposcar_path, new_uc_path; force = true)
    isfile(new_ss_path) || cp(ssposcar_path, new_ss_path; force = true)

    move_ifcs(ce, outpath)

    if isnothing(tep_energies_path)

        # Check for mpirun
        res = run(`which mpirun`; wait = false)
        found_mpirun = success(res)
        found_mpirun || @warn "Could not find mpirun on path, defaulting to 1 thread."

        flags = ""
        flags *= (L === Quantum) ? "--quantum " : ""
        flags *= needs_true_V(ce) ? "--dumpconfigs" : ""

        if found_mpirun
            cmd_str = `mpirun -np $(nthreads) effective_hamiltonian --thirdorder --fourthorder --nconf $(ce.nconf) --temperature $(T) $(flags)`
        else
            cmd_str = `effective_hamiltonian --thirdorder --fourthorder --nconf $(ce.nconf) --temperature $(T) $(flags)`
        end
        @info cmd_str
        cd(outpath) do 
            run(cmd_str)
        end

        # Need to use the HDF5 output, this gurantees
        # order of positions and energies matches.
        tep_energies_path = joinpath(outpath, "outfile.canonical_configs.hdf5")
        is_hdf5 = Val{true}()
    else
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

    #!TODO VARIANCE ANALYSIS ON RANDOM VARIABLE

    res = bootstrap_corrections(V, V₂, V₃, V₄, T, outpath, ce, n_atoms)

    save.(res, Ref(outpath), Ref(ce.n_boot), Ref(ce.boot_size))

    return res

end

function save(bce::BootstrapCumualantEstimate{L}, outdir::String, n_boot, boot_size) where L
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
        println(f, "# Row 1: Values, Row 2: Standard Error estimated from $(n_boot) bootstraps of size $(boot_size)")
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