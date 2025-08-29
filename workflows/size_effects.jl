using LAMMPS
using TDEP
using CumulantAnalysis
using SimpleCrystals
using SimpleCrystals


expansion_order = 3
# just pick a big number, havent looked at effect of nsteps yet
# probably needs to be bigger for LJ than SW
stdep_samples = 50_000 
is_quantum = false

# Lennard-Jones Argon
# make_outpath = (T,s) -> "/mnt/merged/emeitz/CumulantAnalysisTest/sTDEP_LJ_size_effects/T$(ustrip(T))/$(s)UC"
# temperatures = [10, 80]
# sizes = [3,4,5,6,7,8]
# ucposcar_path = "/home/emeitz/scripts/TDEP/LJ/infile.ucposcar_oneatom"
# make_ifc_path = (T) -> "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/IFCs/T$(T)_0/infile.forceconstant"
# pot_cmds = ["pair_style lj/cut 8.5", "pair_coeff * * 0.010423 3.4", "pair_modify shift yes"]
# make_crystal = (s) -> FCC(5.2468u"angstrom", :Ar, SVector(s,s,s))

# Stillinger-Weber Silicon
make_outpath = (T,s) -> "/mnt/merged/emeitz/CumulantAnalysisTest/sTDEP_SW_size_effects/T$(ustrip(T))/$(s)UC"
temperatures = [100, 1300]
sizes = [2,3,4,5,6,7]
ucposcar_path = "/home/emeitz/scripts/TDEP/SW/infile.ucposcar2"
make_ifc_path = (T) -> "/mnt/merged/emeitz/SW_IFC_NODES/IFCs/T$(T)_0/infile.forceconstant"

sw_pot = "/home/emeitz/software/lammps/potentials/Si.sw"
pot_cmds = ["pair_style sw", "pair_coeff * * \"$(sw_pot)\" Si"]
make_crystal = (s) -> Diamond(5.43u"angstrom", :Si, SVector(s,s,s))

LAMMPS.MPI.Init()


for T in temperatures
    for s in sizes

        outpath = make_outpath(T,s)

        @info "Temperature: $(T), Supercell: $(s)x$(s)x$(s)"
        mkpath(outpath)

        crys = make_crystal(s)

        ssposcar_path = joinpath(outpath, "infile.ssposcar")
        to_ssposcar(crys, ssposcar_path)

        s = TDEPSystem(ssposcar_path)
        calc = LAMMPSCalculator(s, pot_cmds)

        se = sTDEPEstimator(expansion_order, stdep_samples, T, is_quantum)
        
        F_c, S_c, U_c, Cv_c = estimate(
            se,
            calc,
            outpath;
            ucposcar_path = ucposcar_path,
            ssposcar_path = ssposcar_path,
            ifc_path = make_ifc_path(T),
            n_boot = 200,
        )

        save_errors(F_c,  outpath)
        save_errors(S_c,  outpath)
        save_errors(U_c,  outpath)
        save_errors(Cv_c, outpath)
    end
end