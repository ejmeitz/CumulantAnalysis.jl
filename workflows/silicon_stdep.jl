# Runs sTDEP for all the lattice parameter and temperature combinations
using CumulantAnalysis

Ts = [1400,1500,1600]
sw_pot = "/home/emeitz/software/lammps/potentials/Si.sw"
pot_cmds = ["pair_style sw", "pair_coeff * * \"$(sw_pot)\" Si"]
r_cut = 4.0

quantum = false
n_lattice_params = 10
n_iter = 10
maximum_frequency = 18 # THz

base_outpath = (T) -> "/mnt/merged/emeitz/CumulantAnalysisTest/Silicon_HighT_sTDEP_IFCs/T$(T)"
structure_basepath = "/mnt/merged/emeitz/SW_sTDEP"

for T in Ts
    
    make_stdep_ifcs(
        joinpath(structure_basepath, "infile.ucposcar"),
        joinpath(structure_basepath, "infile.ssposcar"),
        base_outpath,
        pot_cmds,
        n_iter,
        r_cut,
        T,
        maximum_frequency,
        quantum
    )
end
