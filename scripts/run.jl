using JLD2
using Plots
using CumulantAnalysis
using Statistics
using Measurements
using CSV

basepath = "/mnt/merged/emeitz/test/LJCumulantData6UC"
# basepath = "/mnt/merged/emeitz/test/SWCumulantLang5UC"
#tep_path = (T) -> "/mnt/mntsdb/emeitz/ForceConstants/SW_ALM/SW_$(T)K_residual.jld2"
tep_path = (T) -> "/mnt/mntsdb/emeitz/ForceConstants/LJ_TDEP_80K_6UC.jld2"
# tep_path = (T) -> "/mnt/mntsdb/emeitz/ForceConstants/SW_TDEP_5UC_1300K.jld2"
ld_path = (T,s) -> joinpath(basepath, "T$(T)", "seed$(s)", "dump.positions_unrolled")
stat_path = (T,s) -> joinpath(basepath, "T$(T)", "seed$(s)", "dump.stat")
eq_path = (T,s) -> joinpath(basepath, "T$(T)", "seed$(s)", "equilibrium.atom")
outdir = (T) -> joinpath(basepath, "T$(T)")

kB = 8.617333e-5 # eV / K
hbar = 6.582119569e-16 # eV-s
n_seeds = 25
order = 3

# Lennard-Jones Argon
ifc_conv = 23.060541945 # converts kcal/mol/Ang^2 to eV/Ang^2
f_conv = sqrt(418.4) * 1e12 # converts freqs from real units --> rad /s 
m = 39.95
temperatures = [80, 10]


# Stillinger-Weber Silicon
# ifc_conv = 1.0
# f_conv = 9.82269474855602e13 # converts freqs from metal units --> rad / s
# m = 28.085
# temperatures = [1300]


for T in temperatures

    F_corrections = CumulantCorrections{order}[]
    S_corrections = CumulantCorrections{order}[]
    U_corrections = CumulantCorrections{order}[]
    Cv_corrections = CumulantCorrections{order}[]

    F_running = Measurement[]
    S_running = Measurement[]
    U_running = Measurement[]
    Cv_running = Measurement[]

    U_actual = 0.0
    Cv_actual = 0.0

    freqs_sq, dynmat = load(tep_path(T), "freqs_sq", "dynmat");
    freqs = f_conv .* sqrt.(freqs_sq);
    F2 = dynmat .* m  ./ ifc_conv

    for s in range(0,n_seeds-1)
        F_c, S_c, U_c, Cv_c = 
            estimate_thermo_properties(
                ld_path(T,s),
                eq_path(T,s),
                stat_path(T,s),
                F2,
                outdir(T),
                freqs, kB, hbar, T;
                order = order
            )

        U_actual += U_c.ground_truth
        Cv_actual += Cv_c.ground_truth

        push!(F_corrections, F_c)
        F_current = measurement(CumulantAnalysis.mean(F_corrections), CumulantAnalysis.se(F_corrections))
        push!(F_running, F_current)

        push!(S_corrections, S_c)
        S_current = measurement(CumulantAnalysis.mean(S_corrections), CumulantAnalysis.se(S_corrections))
        push!(S_running, S_current)

        push!(U_corrections, U_c)
        U_current = measurement(CumulantAnalysis.mean(U_corrections), CumulantAnalysis.se(U_corrections))
        push!(U_running, U_current)

        push!(Cv_corrections, Cv_c)
        Cv_current = measurement(CumulantAnalysis.mean(Cv_corrections), CumulantAnalysis.se(Cv_corrections))
        push!(Cv_running, Cv_current)
    end

    U_actual /= n_seeds
    Cv_actual /= n_seeds
  
    CumulantAnalysis.save_errors(F_corrections, outdir(T))
    CumulantAnalysis.save_errors(S_corrections, outdir(T))
    CumulantAnalysis.save_errors(U_corrections, outdir(T))
    CumulantAnalysis.save_errors(Cv_corrections, outdir(T))    

    props = Dict(
        "F"  => F_running,
        "S"  => S_running,
        "U"  => U_running,
        "Cv" => Cv_running,
    )

    # Plot running average and SE for each observable
    for (name, arr) in props
        seeds = 1:length(arr)
        y     = getproperty.(arr, :val)
        yerr  = getproperty.(arr, :err)
    
        p = scatter(
          seeds, y;
          yerr      = yerr,
          xlabel    = "N Seeds",
          ylabel    = name,
          title     = "Running average of $name",
          legend    = :topright,
        )
    
        # add horizontal â€œtrueâ€ line for U and Cv
        if name == "U"
            hline!(p, [U_actual];
                label      = "True U",
                linestyle  = :dash,
                linewidth  = 2,
            )
        elseif name == "Cv"
            hline!(p, [Cv_actual];
                label      = "True Cv",
                linestyle  = :dash,
                linewidth  = 2,
            )
        end
    
        savefig(p, joinpath(outdir(T), "$(name)_running.png"))
    end
end
