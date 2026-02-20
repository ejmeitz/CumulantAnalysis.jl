using CumulantAnalysis


nconf = 100_000 
nboot = 5000

quantum = true
base_outpath = "/mnt/merged/emeitz/CumulantAnalysisTest/Neon_ANALYTICAL_PIMD"
getoutpath = (T) -> joinpath(base_outpath,  "T$(T)")

Ts = [4,6,8,10,12,14,16,18,20,22,24]
pot_cmds = ["pair_style lj/cut 6.955", "pair_coeff * * 0.0032135 2.782", "pair_modify shift yes"]

for a in trial_lattice_constants

    ## sTDEP IFCs
    ifc_basepath = (T) -> "/mnt/merged/emeitz/Neon_sTDEP/results/T$(T)"
    ucposcar_path = (T) -> joinpath(ifc_basepath(T), "infile.ucposcar")
    ssposcar_path = (T) -> joinpath(ifc_basepath(T), "infile.ssposcar")
    ifc2_path = (T) -> joinpath(ifc_basepath(T), "infile.forceconstant")
    ifc3_path = (T) -> joinpath(ifc_basepath(T), "infile.forceconstant_thirdorder")
    ifc4_path = (T) -> joinpath(ifc_basepath(T), "infile.forceconstant_fourthorder")

    crystal_thermodynamic_properties(
        Ts,
        getoutpath,
        ucposcar_path,
        ssposcar_path,
        ifc2_path,
        ifc3_path,
        ifc4_path,
        pot_cmds;
        quantum = quantum,
        nconf = nconf,
        nboot = nboot,
        size_study = false
    )
end

