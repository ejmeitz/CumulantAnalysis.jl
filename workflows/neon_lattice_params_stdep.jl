# Runs sTDEP for all the lattice parameter and temperature combinations

quantum = true
base_outpath = "/mnt/merged/emeitz/CumulantAnalysisTest/Neon_ThermalExpansion_IFCs_sTDEP"
getoutpath = (T) -> joinpath(base_outpath,  "T$(T)")

Ts = [4,6,8,10,12,14,16,18,20,22,24]
pot_cmds = ["pair_style lj/cut 6.955", "pair_coeff * * 0.0032135 2.782", "pair_modify shift yes"]


n_lattice_params = 8
n_iter = 9
maximum_frequency = 2.5

structures_basepath = (T,a) = "/mnt/merged/emeitz/CumulantAnalysisTest/thermal_expansion_series/T$(T)/a$(a)"


for i in 1:n_lattice_params

    outpath_a = joinpath(base_outpath, "a$(i)")

    for T in Ts
        sb = structures_basepath(T, i)
        outpath_T = joinpath(outpath_a, "T$(T)")
        mkpath(outpath_T)
        make_stdep_ifcs(
            joinpath(sb, "infile.ucposcar"),
            joinpath(sb, "infile.ucposcar"),
            outpath_T,
            pot_cmds,
            n_iter,
            maximum_frequency,
            quantum,
        )
    end
end