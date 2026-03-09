using CumulantAnalysis

# Could technically just remap the 4 UC IFCs, but this is easy enough

T = 80
Ns = [5,6] # already have 4 unitcell data


quantum = false
n_iter = 10
maximum_frequency = 2.5 # THz

base_outpath = (N) -> "/mnt/merged/emeitz/CumulantAnalysisTest/LJ_sTDEP_SIZE_EFFECTS/N$(N)_T$(T)"

ucposcar_path = (T) -> "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/infile.ucposcar"
ssposcar_path = (T) -> "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/infile.ssposcar" #4UC
pot_cmds = ["pair_style lj/cut 8.5", "pair_coeff * * 0.010423 3.4", "pair_modify shift yes"]

structure_basepath = (N) -> "/home/emeitz/scripts/TDEP/LJ/size_effects/80K/$(N)UC"

for T in Ts
    make_stdep_ifcs(
        joinpath(structure_basepath, "infile.ucposcar"),
        joinpath(structure_basepath, "infile.ssposcar"),
        base_outpath(T),
        pot_cmds,
        n_iter,
        r_cut,
        T,
        maximum_frequency,
        quantum
    )
end
