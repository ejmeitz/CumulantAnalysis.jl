export estimate, save


function estimate(
        ce::SamplingCumulantEstimator{O},
        T::Real, # Kelvin
        outpath::String,
        lammps_potential_cmds::Union{String, Vector{String}};
        ucposcar_path::String = joinpath(outpath, "infile.ucposcar"),
        ssposcar_path::String = joinpath(outpath, "infile.ssposcar"),
        n_threads::Int = Threads.nthreads(),
        size_study::Bool = false,
        quantum::Bool = false,
        q_mesh::AbstractVector{<:Integer} = [30,30,30],
    ) where {O}

    @assert O <= 2 "Up to second order cumulant corrections are supported. Asked for $(O)."

    T = Float64(T)
    
    if !needs_true_V(ce)
        error("The provided SamplingCumulantEstimator, $(typeof(ce)), does not use the true potential energies V. You probably didnt meant to call this method.")
    end

    isfile(ucposcar_path) || throw(ArgumentError("ucposcar path is not a file: $(ucposcar_path)"))
    isfile(ssposcar_path) || throw(ArgumentError("ssposcar path is not a file: $(ssposcar_path)"))

    new_uc_path = joinpath(outpath, "infile.ucposcar")
    new_ss_path = joinpath(outpath, "infile.ssposcar")
    isfile(new_uc_path) || cp(ucposcar_path, new_uc_path; force = true)
    isfile(new_ss_path) || cp(ssposcar_path, new_ss_path; force = true)

    uc = CrystalStructure(ucposcar_path)
    sc = CrystalStructure(ssposcar_path)
    n_atoms = length(sc)

    L = ifelse(quantum, Quantum, Classical)
    settings = ConfigSettings(ce.nconf, T, L)

    move_ifcs(ce, outpath)
    ifcs = load_ifcs(ce, ucposcar_path, outpath)
    ifc_kwargs = LatticeDynamicsToolkit.build_kwargs(ifcs...)
    
    make_calc = (sc) -> LAMMPSCalculator(sc, lammps_potential_cmds)

    @info "Calculating Harmonic Properties"
    if is_amorphous(ce)
        F₀, S₀, U₀, Cᵥ₀ = Hartree_to_eV .* harmonic_properties(T, ifc_kwargs.ifc2, sc, L;
                                                 n_threads = n_threads)
    else
        F₀, S₀, U₀, Cᵥ₀ = Hartree_to_eV .* harmonic_properties(T, uc, ifc_kwargs.ifc2, q_mesh, L;
                                                 n_threads = n_threads)
    end


    tep_energies, V = make_energy_dataset(
        settings,
        uc,
        sc,
        make_calc;
        ifc_kwargs...,
        n_threads = n_threads
    )

    V₂ = getindex.(tep_energies, 1)
    V₃ = getindex.(tep_energies, 2)
    V₄ = getindex.(tep_energies, 3)
    Vₚ = zeros(length(V))

    header = ["V" "Vp" "V2" "V3" "V4"]
    str_fmt_str = (N) -> Printf.Format(join(fill("%15s", N), " "))
    open(joinpath(outpath, "outfile.all_energies"), "w") do f
        println(f, "# Energies in eV, N atoms = $(n_atoms), tempearture [K] = $(T)")
        println(f, Printf.format(str_fmt_str(length(header)), header...))
        writedlm(f, [V Vₚ V₂ V₃ V₄])
    end

    res = bootstrap_corrections(
            V, V₂, V₃, V₄, T,
            outpath,
            F₀, S₀, U₀, Cᵥ₀,
            ce,
            n_atoms,
            L
        )

    save.(res, Ref(outpath), Ref(ce.n_boot))

    # Compute some statistics to assess convergence with N
    size_study && do_size_study(ce, outpath, V, V₂, V₃, V₄, T, n_atoms)

    return res
end

function estimate(
        ce::AnalyticalEstimator,
        T::Real, # Kelvin
        outpath::String,
        lammps_potential_cmds::Union{String, Vector{String}};
        ucposcar_path::String = joinpath(outpath, "infile.ucposcar"),
        ssposcar_path::String = joinpath(outpath, "infile.ssposcar"),
        n_threads::Int = Threads.nthreads(),
        size_study::Bool = false,
        quantum::Bool = false,
        harmonic_q_mesh::AbstractVector{<:Integer} = [30,30,30],
        free_energy_q_mesh::AbstractVector{<:Integer} = [25,25,25]
    )

    T = Float64(T)
    
    if !needs_true_V(ce)
        error("The provided SamplingCumulantEstimator, $(typeof(ce)), does not use the true potential energies V. You probably didnt meant to call this method.")
    end

    isfile(ucposcar_path) || throw(ArgumentError("ucposcar path is not a file: $(ucposcar_path)"))
    isfile(ssposcar_path) || throw(ArgumentError("ssposcar path is not a file: $(ssposcar_path)"))

    new_uc_path = joinpath(outpath, "infile.ucposcar")
    new_ss_path = joinpath(outpath, "infile.ssposcar")
    isfile(new_uc_path) || cp(ucposcar_path, new_uc_path; force = true)
    isfile(new_ss_path) || cp(ssposcar_path, new_ss_path; force = true)

    uc = CrystalStructure(ucposcar_path)
    sc = CrystalStructure(ssposcar_path)
    n_atoms = length(sc)

    L = ifelse(quantum, Quantum, Classical)
    settings = ConfigSettings(ce.nconf, T, L)

    move_ifcs(ce, outpath)
    ifcs = load_ifcs(ce, ucposcar_path, outpath)
    ifc_kwargs = LatticeDynamicsToolkit.build_kwargs(ifcs...)
    
    make_calc = (sc) -> LAMMPSCalculator(sc, lammps_potential_cmds)

    @info "Calculating Harmonic Properties"
    F₀, S₀, U₀, Cᵥ₀ = Hartree_to_eV .* harmonic_properties(T, uc, ifc_kwargs.ifc2, harmonic_q_mesh, L;
                                                 n_threads = n_threads)

    tep_energies, V = make_energy_dataset(
        settings,
        uc,
        sc,
        make_calc;
        ifc_kwargs...,
        n_threads = n_threads
    )

    V₂ = getindex.(tep_energies, 1)
    V₃ = getindex.(tep_energies, 2)
    V₄ = getindex.(tep_energies, 3)
    Vₚ = zeros(length(V))

    header = ["V" "Vp" "V2" "V3" "V4"]
    str_fmt_str = (N) -> Printf.Format(join(fill("%15s", N), " "))
    open(joinpath(outpath, "outfile.all_energies"), "w") do f
        println(f, "# Energies in eV, N atoms = $(n_atoms), tempearture [K] = $(T)")
        println(f, Printf.format(str_fmt_str(length(header)), header...))
        writedlm(f, [V Vₚ V₂ V₃ V₄])
    end

    # Get analytical corrections
    analytical_corrections = free_energy_corrections(
        T,
        uc,
        ifc_kwargs.ifc2,
        ifc_kwargs.ifc3,
        ifc_kwargs.ifc4;
        mesh = free_energy_q_mesh,
        quantum = quantum,
        n_threads = n_threads
    )

    res = bootstrap_corrections(
            V, V₂, V₃, V₄, T,
            outpath,
            F₀, S₀, U₀, Cᵥ₀,
            analytical_corrections,
            ce,
            n_atoms,
            L
        )

    save.(res, Ref(outpath), Ref(ce.n_boot))

    # Compute some statistics to assess convergence with N
    size_study && do_size_study(ce, outpath, V, V₂, V₃, V₄, T, n_atoms)

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



# function parse_energies(path, is_hdf5::Val{false})

#     data = readdlm(path, comments = true)

#     # Parse n_atoms and T from header
#     f = open(path, "r")
#     readline(f) # skip
#     T = parse(Float64, split(strip(readline(f)))[end])
#     n_atoms = parse(Int, split(strip(readline(f)))[end])
#     close(f)

#     # energies in file are meV / atom, Converts to eV
#     conv = n_atoms / 1000
    
#     @views Vₚ = data[:, 3] .* conv
#     @views V₂ = data[:, 4] .* conv
#     @views V₃ = data[:, 5] .* conv
#     @views V₄ = data[:, 6] .* conv

#     return T, n_atoms, Vₚ, V₂ , V₃, V₄
# end

# function parse_energies(path, is_hdf5::Val{true})

#     f = h5open(path, "r")

#     Vₚ = read(f, "polar_potential_energy")
#     V₂ = read(f, "secondorder_potential_energy")
#     V₃ = read(f, "thirdorder_potential_energy")
#     V₄ = read(f, "fourthorder_potential_energy")

#     T = Float64(attrs(f)["temperature_thermostat"][1])
#     n_atoms = Int64(attrs(f)["number_of_atoms"][1])

#     close(f)
    
#     return T, n_atoms, Vₚ, V₂, V₃, V₄
# end


# function calculate_true_energies(calc, nconf, ssposcar_path, outpath)
#     isnothing(calc) && throw(ArgumentError("This SamplingCumulantEstimator requires a calculator for energies, but `calc` was nothing."))
#     (AtomsCalculators.energy_unit(calc) != u"eV") && throw(ArgumentError("Expected calculator with eV energy units. Got $(AtomsCalculators.energy_unit(calc)) "))

#     sys_ss = TDEPSystem(ssposcar_path)
#     V = zeros(nconf)

#     p = Progress(nconf, desc = "Calculating True Energies")
#     L = typeof(1.0u"Å")

#     configs_path = joinpath(outpath, "outfile.canonical_configs.hdf5")
#     configs = h5read(configs_path, "positions")

#     for i in 1:nconf   
#         @views c = configs[:,:,i]
#         #! CAN I REMOVE THIS ALLOCATION?
#         new_sys = TDEPSystem(sys_ss, vec(collect(reinterpret(SVector{3, L}, c))))

#         V[i] = ustrip(AtomsCalculators.potential_energy(new_sys, calc))

#         next!(p)
#     end
#     finish!(p)

#     return V
# end

# function estimate(
#         ce::SamplingCumulantEstimator{O,L},
#         T::Real, # kelvin
#         outpath::String;
#         ucposcar_path::String = joinpath(outpath, "infile.ucposcar"),
#         ssposcar_path::String = joinpath(outpath, "infile.ssposcar"),
#         nthreads::Int = Threads.nthreads(),
#         rm_configs::Bool = true,
#         use_control_variates::Bool = true
#     ) where {O, L <: Limit}

#     @assert O <= 2 "Up to second order cumulant corrections are supported. Asked for $(O)."

#     isfile(ucposcar_path) || throw(ArgumentError("ucposcar path is not a file: $(ucposcar_path)"))
#     isfile(ssposcar_path) || throw(ArgumentError("ssposcar path is not a file: $(ssposcar_path)"))

#     new_uc_path = joinpath(outpath, "infile.ucposcar")
#     new_ss_path = joinpath(outpath, "infile.ssposcar")
#     isfile(new_uc_path) || cp(ucposcar_path, new_uc_path; force = true)
#     isfile(new_ss_path) || cp(ssposcar_path, new_ss_path; force = true)

#     move_ifcs(ce, outpath)

#     # Check for mpirun
#     res = run(`which mpirun`; wait = false)
#     found_mpirun = success(res)
#     found_mpirun || @warn "Could not find mpirun on path, defaulting to 1 thread."

#     flags = ""
#     flags *= (ce isa HarmonicEstimator) ? "" : "--thirdorder --fourthorder "
#     flags *= (L === Quantum) ? "--quantum " : ""
#     flags *= needs_true_V(ce) ? "--dumpconfigs " : ""

#     if found_mpirun
#         cmd_str = `mpirun -np $(nthreads) effective_hamiltonian --nconf $(ce.nconf) --temperature $(T) $(split(flags))`
#     else
#         cmd_str = `effective_hamiltonian --nconf $(ce.nconf) --temperature $(T) $(split(flags))`
#     end
#     @info cmd_str
#     cd(outpath) do 
#         run(cmd_str)
#     end

#     # If configs are dumped to calculate V, use HDF5
#     # this gurantees order of positions and energies matches.
#     if needs_true_V(ce)
#         tep_energies_path = joinpath(outpath, "outfile.canonical_configs.hdf5")
#         is_hdf5 = Val{true}()
#     else
#         tep_energies_path = joinpath(outpath, "outfile.energies")
#         is_hdf5 = Val{false}()
#     end

#     @info "Parsing Energies"
#     T_file, n_atoms, Vₚ, V₂, V₃, V₄ = parse_energies(tep_energies_path, is_hdf5)

#     if T_file != T
#         @warn "You said tempearture was $T, but parsed temperature as $(T_file) from outfile.energies"
#     end

#     if !(sum(Vₚ) ≈ 0.0)
#         error("Got non-zero polar term. I don't know what to do with that.")
#     end

#     if needs_true_V(ce)
#         V = calculate_true_energies(ce.force_calculator, ce.nconf, ssposcar_path, outpath)   
#     else
#         V = NaN .* zeros(eltype(V₂), size(V₂))
#     end

#     #!TODO ADD TEMP/n_ATOMS DATA AND UNITS TO FILE
#     # This energy file will have the correct ordering it is not always correct 
#     # to load V from outfile.true_potential_energy and V2/V3/V4 from outfile.energies
#     header = ["V" "Vp" "V2" "V3" "V4"]
#     str_fmt_str = (N) -> Printf.Format(join(fill("%15s", N), " "))
#     open(joinpath(outpath, "outfile.all_energies"), "w") do f
#         println(f, Printf.format(str_fmt_str(length(header)), header...))
#         writedlm(f, [V Vₚ V₂ V₃ V₄])
#     end

#     if rm_configs
#         rm(joinpath(outpath, "outfile.canonical_configs.hdf5"); force = true)
#     end

#     res = bootstrap_corrections(V, V₂, V₃, V₄, T, outpath, ce, n_atoms, use_control_variates)
#     save.(res, Ref(outpath), Ref(ce.n_boot))

#     # Compute some statistics to assess convergence with N
#     do_size_study(ce, outpath, V, V₂, V₃, V₄, T, n_atoms, use_control_variates)

#     return res

# end


