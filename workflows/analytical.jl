using CumulantAnalysis


nconf = 100_000 
nboot = 5000

base_outpath = "/mnt/merged/emeitz/CumulantAnalysisTest/SW_ANALYTICAL"
getoutpath = (T) -> joinpath(base_outpath,  "T$(T)")

# LENNARD JONES
# pot = "LJ"
# Ts = [10, 20, 30, 40, 50, 60, 70, 80]
# ucposcar_path = "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/infile.ucposcar"
# ssposcar_path = "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/infile.ssposcar" #4UC
# ifc2_path = (T) -> "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/IFCs/T$(T)_0/infile.forceconstant"
# ifc3_path = (T) -> "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/IFCs/T$(T)_0/infile.forceconstant_thirdorder"
# ifc4_path = (T) -> "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/IFCs/T$(T)_0/infile.forceconstant_fourthorder"
# pot_cmds = ["pair_style lj/cut 8.5", "pair_coeff * * 0.010423 3.4", "pair_modify shift yes"]

pot = "SW"
Ts = [100,300,500,700,900,1100,1300]
ucposcar_path = raw"/mnt/merged/emeitz/SW_IFC_NODES/infile.ucposcar"
ssposcar_path = raw"/mnt/merged/emeitz/SW_IFC_NODES/infile.ssposcar" #3UC
ifc2_path = (T) -> "/mnt/merged/emeitz/SW_IFC_NODES/IFCs/T$(T)_0/infile.forceconstant"
ifc3_path = (T) -> "/mnt/merged/emeitz/SW_IFC_NODES/IFCs/T$(T)_0/infile.forceconstant_thirdorder"
ifc4_path = (T) -> "/mnt/merged/emeitz/SW_IFC_NODES/IFCs/T$(T)_0/infile.forceconstant_fourthorder"

sw_pot = "/home/emeitz/software/lammps/potentials/SW/Si.sw"
pot_cmds = ["pair_style sw", "pair_coeff * * \"$(sw_pot)\" Si"]

for T in Ts
    @info "T = $(T)"

    o = getoutpath(T)
    mkpath(o)

    estim = CumulantAnalysis.AnalyticalEstimator(
            ifc2_path(T), ifc3_path(T), ifc4_path(T), nconf, nboot
        )

    estimate(
            estim,
            T, 
            o,
            pot_cmds;
            ucposcar_path = ucposcar_path,
            ssposcar_path = ssposcar_path,
            size_study = true
        )
end