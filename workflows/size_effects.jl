using LAMMPS
using TDEP
using CumulantAnalysis
using SimpleCrystals

# Explicitly loop over sizes, each simulation will also run a study
# to determine the effect of sample size on convergence

expansion_order = 2
limit = Classical()
samples = 100_000 
n_boot = 500
boot_size = 10_000

# Lennard-Jones Argon
make_outpath = (CT, T,s) -> "/mnt/merged/emeitz/CumulantAnalysisTest/LJ_size_effects_$(CT)/T$(ustrip(T))/$(s)UC"
# make_outpath = (CT, T, s) -> "/mnt/merged/emeitz/CumulantAnalysisTest/LJ_TEST_$(CT)/T$(ustrip(T))/$(s)UC"
temperatures = [10, 80]
sizes = [3,4,5,6]
U0 = Dict(10 => -0.0777399947, 80 => -0.0780397902)
ucposcar_path = "/home/emeitz/scripts/TDEP/LJ/infile.ucposcar_oneatom"
ifc2_path = (T) -> "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/IFCs/T$(T)_0/infile.forceconstant"
ifc3_path = (T) -> "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/IFCs/T$(T)_0/infile.forceconstant_thirdorder"
ifc4_path = (T) -> "/mnt/merged/emeitz/LJ_IFC_INTERPOLATION_NODES_FINE/IFCs/T$(T)_0/infile.forceconstant_fourthorder"

pot_cmds = ["pair_style lj/cut 8.5", "pair_coeff * * 0.010423 3.4", "pair_modify shift yes"]
make_crystal = (s) -> FCC(5.2468u"angstrom", :Ar, SVector(s,s,s))

# Stillinger-Weber Silicon
# make_outpath = (CT, T,s) -> "/mnt/merged/emeitz/CumulantAnalysisTest/SW_size_effects_$(CT)/T$(ustrip(T))/$(s)UC"
# temperatures = [100, 1300]
# sizes = [2,3,4,5,6]
# ucposcar_path = "/home/emeitz/scripts/TDEP/SW/infile.ucposcar2"
# ifc2_path = (T) -> "/mnt/merged/emeitz/SW_IFC_NODES/IFCs/T$(T)_0/infile.forceconstant"
# ifc3_path = (T) -> "/mnt/merged/emeitz/SW_IFC_NODES/IFCs/T$(T)_0/infile.forceconstant_thirdorder"
# ifc4_path = (T) -> "/mnt/merged/emeitz/SW_IFC_NODES/IFCs/T$(T)_0/infile.forceconstant_fourthorder"
# U0 = Dict(100 => -4.3365637, 1300 => -4.3247041)

# sw_pot = "/home/emeitz/software/lammps/potentials/Si.sw"
# pot_cmds = ["pair_style sw", "pair_coeff * * \"$(sw_pot)\" Si"]
# make_crystal = (s) -> Diamond(5.43u"angstrom", :Si, SVector(s,s,s))

names_dict = Dict(
    EffectiveHamiltonianEstimator => "MDV0",
    HarmonicEstimator => "HarmV0",
    FourthOrderEstimator => "QuarticV0"
)

function make_estimator(::Type{EffectiveHamiltonianEstimator}, calc, T, n_atoms)
    return EffectiveHamiltonianEstimator(
            expansion_order, 
            limit, 
            ifc2_path(T), 
            ifc3_path(T), 
            ifc4_path(T), 
            U0[T] * n_atoms,
            samples,
            n_boot,
            boot_size
        )
end

function make_estimator(::Type{HarmonicEstimator}, calc, T, n_atoms)
    return HarmonicEstimator(
            expansion_order,
            limit, 
            calc,
            ifc2_path(T),
            samples,
            n_boot,
            boot_size
        )
end

function make_estimator(::Type{FourthOrderEstimator}, calc, T, n_atoms)
    return FourthOrderEstimator(
            expansion_order,
            limit, 
            calc,
            ifc2_path(T),
            ifc3_path(T),
            ifc4_path(T),
            samples,
            n_boot,
            boot_size
        )
end

LAMMPS.MPI.Init()

for CE_TYPE in (EffectiveHamiltonianEstimator, HarmonicEstimator, FourthOrderEstimator)
    for T in temperatures
        for s in sizes

            crys = make_crystal(s)

            outpath = make_outpath(names_dict[CE_TYPE], T, s)

            @info "Temperature: $(T), Supercell: $(s)x$(s)x$(s)"
            mkpath(outpath)

            ssposcar_path = joinpath(outpath, "infile.ssposcar")
            to_ssposcar(crys, ssposcar_path)

            s = TDEPSystem(ssposcar_path)
            calc = LAMMPSCalculator(s, pot_cmds)
            
            e = make_estimator(CE_TYPE, calc, T, length(crys))
            
            estimate(
                e,
                T,
                outpath;
                ucposcar_path = ucposcar_path,
                ssposcar_path = ssposcar_path,
            )

        end
    end
end