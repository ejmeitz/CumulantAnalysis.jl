using CumulantAnalysis

nconf = 100_000 
nboot = 5_000

quantum = false
base_outpath = "/mnt/merged/emeitz/CumulantAnalysisTest/SW_ANALYTICAL"
getoutpath = (T) -> joinpath(base_outpath,  "T$(T)")

pot = "SW"
Ts = [100,300,500,700,900,1100,1300]
ucposcar_path = raw"/mnt/merged/emeitz/SW_IFC_NODES/infile.ucposcar"
ssposcar_path = raw"/mnt/merged/emeitz/SW_IFC_NODES/infile.ssposcar" #3UC
ifc2_path = (T) -> "/mnt/merged/emeitz/SW_IFC_NODES/IFCs/T$(T)_0/infile.forceconstant"
ifc3_path = (T) -> "/mnt/merged/emeitz/SW_IFC_NODES/IFCs/T$(T)_0/infile.forceconstant_thirdorder"
ifc4_path = (T) -> "/mnt/merged/emeitz/SW_IFC_NODES/IFCs/T$(T)_0/infile.forceconstant_fourthorder"

sw_pot = "/home/emeitz/software/lammps/potentials/Si.sw"
pot_cmds = ["pair_style sw", "pair_coeff * * \"$(sw_pot)\" Si"]

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
