using CumulantAnalysis

nconf = 100_000 
nboot = 5_000

quantum = false
base_outpath = "/mnt/merged/emeitz/CumulantAnalysisTest/LJ_ANALYTICAL_sTDEP"
getoutpath = (T) -> joinpath(base_outpath,  "T$(T)")

# LENNARD JONES
Ts = collect(5:5:80)
ucposcar_path = (T) -> "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/infile.ucposcar"
ssposcar_path = (T) -> "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/infile.ssposcar" #4UC
pot_cmds = ["pair_style lj/cut 8.5", "pair_coeff * * 0.010423 3.4", "pair_modify shift yes"]

# MD-TDEP IFCs
# mdtdep_basepath = (T) -> "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/IFCs/T$(T)_0"
# ifc2_path = (T) -> "$(mdtdep_basepath(T))/infile.forceconstant"
# ifc3_path = (T) -> "$(mdtdep_basepath(T))/infile.forceconstant_thirdorder"
# ifc4_path = (T) -> "$(mdtdep_basepath(T))/infile.forceconstant_fourthorder"

# sTDEP IFCs
stdep_basepath = (T) -> "/mnt/merged/emeitz/LJ_sTDEP/RESULTS/T$(T)_0"
ifc2_path = (T) -> "$(stdep_basepath(T))/infile.forceconstant"
ifc3_path = (T) -> "$(stdep_basepath(T))/infile.forceconstant_thirdorder"
ifc4_path = (T) -> "$(stdep_basepath(T))/infile.forceconstant_fourthorder"

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

