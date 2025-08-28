using LAMMPS
using TDEP
using CumulantAnalysis
using SimpleCrystals
using SimpleCrystals

basedir = (T,s) -> "/mnt/merged/emeitz/CumulantAnalysisTest/sTDEP_SW_TEST/T$(ustrip(T))/$(s)UC"

expansion_order = 3
# just pick a big number, havent looked at effect of nsteps yet
# probably needs to be bigger for LJ than SW
stdep_samples = 50_000 
is_quantum = false

# Lennard-Jones Argon
#temperatures = [80, 10]
#ucposcar_path = "/home/emeitz/scripts/TDEP/LJ/infile.ucposcar_oneatom"


# Stillinger-Weber Silicon
temperatures = [100, 1300]
sizes = [2,3,4,5,6,7,8]
ucposcar_path = "/home/emeitz/scripts/TDEP/SW/infile.ucposcar2"
make_ifc_path = (T) -> "/mnt/merged/emeitz/SW_IFC_NODES/IFCs/T$(T)_0/infile.forceconstant"

sw_pot = "/home/emeitz/software/lammps/potentials/Si.sw"
sw_cmds = ["pair_style sw", "pair_coeff * * \"$(sw_pot)\" Si"]

LAMMPS.MPI.Init()


for T in temperatures
    for s in sizes

        @info "Temperature: $(T), Supercell: $(s)x$(s)x$(s)"
        mkpath(basedir(T,s))

        crys = Diamond(5.43u"angstrom", :Si, SVector(s,s,s))
        # crys = FCC(5.2468u"angstrom", :Ar, SVector(s,s,s))

        ssposcar_path = joinpath(basedir(T,s), "infile.ssposcar")
        to_ssposcar(crys, ssposcar_path)
        s = TDEPSystem(ssposcar_path)
        calc = LAMMPSCalculator(s, sw_cmds)

        se = sTDEPEstimator(expansion_order, stdep_samples, T, is_quantum)
        
        F_c, S_c, U_c, Cv_c = estimate(
            se,
            calc,
            basedir(T,s);
            ucposcar_path = ucposcar_path,
            ssposcar_path = ssposcar_path,
            ifc_path = make_ifc_path(T),
            n_boot = 200,
        )

        save_errors(F_c,  basedir(T,s))
        save_errors(S_c,  basedir(T,s))
        save_errors(U_c,  basedir(T,s))
        save_errors(Cv_c, basedir(T,s))
    end
end