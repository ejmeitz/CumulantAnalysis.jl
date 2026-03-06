using CumulantAnalysis


nconf = 100 # doesnt matter just looking at F1/F2 terms
nboot = 10 # doesnt matter just looking at F1/F2 terms
k_meshes = [5, 10, 15, 20, 25, 30, 35]

base_outpath = "/mnt/merged/emeitz/CumulantAnalysisTest/KMESH_CONVERGENCE"

#### NEON ####

quantum = true
Ts = [4, 24]
pot_cmds = ["pair_style lj/cut 6.955", "pair_coeff * * 0.0032135 2.782", "pair_modify shift yes"]

## sTDEP IFCs
ifc_basepath = (T) -> "/mnt/merged/emeitz/Neon_sTDEP/results/T$(T)"
ucposcar_path = (T) -> joinpath(ifc_basepath(T), "infile.ucposcar")
ssposcar_path = (T) -> joinpath(ifc_basepath(T), "infile.ssposcar")
ifc2_path = (T) -> joinpath(ifc_basepath(T), "infile.forceconstant")
ifc3_path = (T) -> joinpath(ifc_basepath(T), "infile.forceconstant_thirdorder")
ifc4_path = (T) -> joinpath(ifc_basepath(T), "infile.forceconstant_fourthorder")

for k in k_meshes
    neon_outpath = (T) -> joinpath(base_outpath, "Neon", "k$(k)", "T$(T)")
    mkpath(neon_outpath(T))

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
        size_study = false,
        free_energy_q_mesh = [k,k,k]
    )
end

#### ARGON ####

quantum = false
Ts = [5, 80]

ucposcar_path = (T) -> "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/infile.ucposcar"
ssposcar_path = (T) -> "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/infile.ssposcar" #4UC
pot_cmds = ["pair_style lj/cut 8.5", "pair_coeff * * 0.010423 3.4", "pair_modify shift yes"]

# sTDEP IFCs
stdep_basepath = (T) -> "/mnt/merged/emeitz/LJ_sTDEP/RESULTS/T$(T)_0"
ifc2_path = (T) -> "$(stdep_basepath(T))/infile.forceconstant"
ifc3_path = (T) -> "$(stdep_basepath(T))/infile.forceconstant_thirdorder"
ifc4_path = (T) -> "$(stdep_basepath(T))/infile.forceconstant_fourthorder"

for k in k_meshes
    argon_outpath = (T) -> joinpath(base_outpath, "Argon", "k$(k)", "T$(T)")
    mkpath(argon_outpath(T))

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
        size_study = false,
        free_energy_q_mesh = [k,k,k]
    )
end

#### SILICON ####

quantum = false
Ts = [100, 1600]

sw_pot = "/home/emeitz/software/lammps/potentials/Si.sw"
pot_cmds = ["pair_style sw", "pair_coeff * * \"$(sw_pot)\" Si"]

# sTDEP IFCs
stdep_basepath = (T) -> "/mnt/merged/emeitz/SW_sTDEP/RESULTS/T$(T)_0"
ucposcar_path = (T) -> joinpath(stdep_basepath(T), "infile.ucposcar")
ssposcar_path = (T) -> joinpath(stdep_basepath(T), "infile.ssposcar")
ifc2_path = (T) -> joinpath(stdep_basepath(T), "infile.forceconstant")
ifc3_path = (T) -> joinpath(stdep_basepath(T), "infile.forceconstant_thirdorder")
ifc4_path = (T) -> joinpath(stdep_basepath(T), "infile.forceconstant_fourthorder")


for k in k_meshes
    silicon_outpath = (T) -> joinpath(base_outpath, "Silicon", "k$(k)", "T$(T)")
    mkpath(silicon_outpath(T))

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
        size_study = false,
        free_energy_q_mesh = [k,k,k]
    )
end