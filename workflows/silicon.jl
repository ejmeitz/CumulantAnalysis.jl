using CumulantAnalysis

nconf = 100_000 
nboot = 5_000

quantum = false
base_outpath = "/mnt/merged/emeitz/CumulantAnalysisTest/SW_ANALYTICAL_sTDEP"
getoutpath = (T) -> joinpath(base_outpath,  "T$(T)")

Ts = [100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600]
sw_pot = "/home/emeitz/software/lammps/potentials/Si.sw"
pot_cmds = ["pair_style sw", "pair_coeff * * \"$(sw_pot)\" Si"]

# sTDEP IFCs
stdep_basepath = (T) -> "/mnt/merged/emeitz/SW_sTDEP/RESULTS/T$(T)_0"
ucposcar_path = (T) -> joinpath(stdep_basepath(T), "infile.ucposcar")
ssposcar_path = (T) -> joinpath(stdep_basepath(T), "infile.ssposcar")
ifc2_path = (T) -> joinpath(stdep_basepath(T), "infile.forceconstant")
ifc3_path = (T) -> joinpath(stdep_basepath(T), "infile.forceconstant_thirdorder")
ifc4_path = (T) -> joinpath(stdep_basepath(T), "infile.forceconstant_fourthorder")


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
    size_study = true
)


# MD-TDEP IFCs
# ucposcar_path = raw"/mnt/merged/emeitz/SW_IFC_NODES/infile.ucposcar"
# ssposcar_path = raw"/mnt/merged/emeitz/SW_IFC_NODES/infile.ssposcar" #3UC
# ifc2_path = (T) -> "/mnt/merged/emeitz/SW_IFC_NODES/IFCs/T$(T)_0/infile.forceconstant"
# ifc3_path = (T) -> "/mnt/merged/emeitz/SW_IFC_NODES/IFCs/T$(T)_0/infile.forceconstant_thirdorder"
# ifc4_path = (T) -> "/mnt/merged/emeitz/SW_IFC_NODES/IFCs/T$(T)_0/infile.forceconstant_fourthorder"
