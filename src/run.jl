export crystal_thermodynamic_properties

function crystal_thermodynamic_properties(
    temperatures::AbstractVector{<:Real},
    outpath::Function,
    ucposcar_path::Function,
    ssposcar_path::Function,
    ifc2_path::Function,
    ifc3_path::Function,
    ifc4_path::Function,
    pot_cmds::Union{String, Vector{String}};
    quantum::Bool = false,
    nconf::Int = 100_000,
    nboot::Int = 2500,
    size_study::Bool = false,
    harmonic_q_mesh::AbstractVector{<:Integer} = [30,30,30],
    free_energy_q_mesh::AbstractVector{<:Integer} = [25,25,25],
    n_threads::Integer = Threads.nthreads(),
    use_hot::Bool = false,
    kwargs...
)
    all_ifcs = Vector{IFC2}(undef, length(temperatures))
    all_ucs = Vector{CrystalStructure}(undef, length(temperatures))

    LIMIT = quantum ? Quantum : Classical

    for (i,T) in enumerate(temperatures)

        mkpath(outpath(T))

        estim = AnalyticalEstimator(
                ifc2_path(T), ifc3_path(T), ifc4_path(T), nconf, nboot
            )

        res, all_ifcs[i] = estimate(
            estim,
            Float64(T),
            outpath(T),
            pot_cmds;
            ucposcar_path = ucposcar_path(T),
            ssposcar_path = ssposcar_path(T),
            n_threads = n_threads,
            size_study = size_study,
            quantum = quantum,
            harmonic_q_mesh = harmonic_q_mesh,
            free_energy_q_mesh = free_energy_q_mesh,
            use_hot = use_hot
        )

        all_ucs[i] = CrystalStructure(ucposcar_path(T))

    end

    if length(temperatures) < 5
        @warn "Fewer than 5 temperatures provided. Temperature dependent corrections may be inaccurate."
    end

    @info "Calculating Temperature Dependent Corrections"

    # Post-Processing Across Temperatures
    # 1) Differentiate freuqencies to get dω/dT
    # 2) Differentiate internal energy constant correction
    #      to get better heat capacity predictions
    S_corr, U_corr, cv_offset_sg = improve_constant_corrections(
        temperatures,
        outpath,
        all_ifcs,
        all_ucs,
        harmonic_q_mesh,
        LIMIT;
        kwargs...
    )

    # Save corrections to files
    for (c, prop) in zip((U_corr, S_corr, cv_offset_sg), ("U", "S", "Cv"))
        p = joinpath(dirname(outpath(1.0)), prop*"_raw_corrections.txt")
        open(p, "w") do f
            writedlm(f, [temperatures c])
        end
    end

    # Make new S/U/Cv files with correction terms
    for (c, prop) in zip((U_corr, S_corr, cv_offset_sg), ("U", "S", "Cv"))
        for (i,T) in enumerate(temperatures)
            p = joinpath(outpath(T), prop*"_mean.h5")
            h5open(p, "r") do f
                harm = read(f, prop*"0")
                corrs = SVector{3}(
                    read(f, prop*"_offset"),
                    read(f, prop*"1"),
                    read(f, prop*"2"),
                )
                corr_SEs = SVector{3}(
                    read(f, prop*"_offset_SE"),
                    read(f, prop*"1_SE"),
                    read(f, prop*"2_SE"),
                )
                total = read(f, prop*"_total")
                total_SE = read(f, prop*"_total_SE")
                #! this definitely does not propagate error correctly
                if prop == "Cv"
                    new_total = total - corrs[1] + c[i]
                    bce = BootstrapCumualantEstimate(
                        harm, SVector(c[i], corrs[2], corrs[3]),
                        corr_SEs, new_total,
                        total_SE, prop * "_corrected", "[kB / atom]"
                    )
                else # U / S
                    units = prop == "U" ? "[eV/atom]" : "[kB / atom]"
                    bce = BootstrapCumualantEstimate(
                        harm + c[i], corrs, corr_SEs, total + c[i],
                        total_SE, prop * "_corrected", units
                    )
                end
                save(bce, outpath(T), nboot)
            end
        end
    end

end