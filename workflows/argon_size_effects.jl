using CumulantAnalysis

# Could technically just remap the 4 UC IFCs, but this is easy enough

T = 80
# Ns = [5,6] # already have 4 unitcell data


quantum = false
# n_iter = 10
# maximum_frequency = 2.5 # THz
r_cut = 8.5

pot_cmds = ["pair_style lj/cut $(r_cut)", "pair_coeff * * 0.010423 3.4", "pair_modify shift yes"]

structure_basepath = (N) -> "/home/emeitz/scripts/TDEP/LJ/size_effects/80K/$(N)UC"
ucposcar_path = (N) -> joinpath(structure_basepath(N), "infile.ucposcar")
ssposcar_path = (N) -> joinpath(structure_basepath(N), "infile.ssposcar")

# ifc_outpath = (N) -> "/mnt/merged/emeitz/CumulantAnalysisTest/LJ_sTDEP_SIZE_EFFECTS/N$(N)_T$(T)"

# for N in Ns
#     make_stdep_ifcs(
#         ucposcar_path(N),
#         ssposcar_path(N),
#         ifc_outpath(N),
#         pot_cmds,
#         n_iter,
#         r_cut,
#         T,
#         maximum_frequency,
#         quantum
#     )
# end

### GET THERMO PROPERTIES
nconf = 100_000 
nboot = 5_000
Ns = [4,5,6]

# sTDEP IFCs
stdep_basepath = (N) -> "/mnt/merged/emeitz/CumulantAnalysisTest/LJ_sTDEP_SIZE_EFFECTS/N$(N)_T80/iter009"
ifc2_path = (N) -> "$(stdep_basepath(N))/infile.forceconstant"
ifc3_path = (N) -> "$(stdep_basepath(N))/infile.forceconstant_thirdorder"
ifc4_path = (N) -> "$(stdep_basepath(N))/infile.forceconstant_fourthorder"

thermo_outpath = (N) -> joinpath("/mnt/merged/emeitz/CumulantAnalysisTest/LJ_SIZE_EFFECTS_80K/", "N$(N)")
for N in Ns

    crystal_thermodynamic_properties(
        [T],
        (T) -> thermo_outpath(N),
        (T) -> ucposcar_path(N),
        (T) -> ssposcar_path(N),
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