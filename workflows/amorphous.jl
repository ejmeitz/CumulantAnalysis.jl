using CumulantAnalysis


nconf = 100_000 
nboot = 5000

base_outpath = "/mnt/merged/emeitz/CumulantAnalysisTest/aSi_Cumulant"
getoutpath = (T) -> joinpath(base_outpath,  "T$(T)")

pot = "SW"
Ts = [100]
# Ts = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
ucposcar_path = raw"/home/emeitz/scripts/aSi/N192/infile.ucposcar"
ssposcar_path = raw"/home/emeitz/scripts/aSi/N192/infile.ssposcar"
ifc2_path = (T) -> "/mnt/mntsdb/emeitz/ForceConstants/aSi/IFC2_T$(T)_N192.h5"

sw_pot = "/home/emeitz/software/lammps/potentials/Si.sw"
pot_cmds = ["pair_style sw", "pair_coeff * * \"$(sw_pot)\" Si"]

for T in Ts
    @info "T = $(T)"

    o = getoutpath(T)
    mkpath(o)

    estim = CumulantAnalysis.AmorphousEstimator(
            2, ifc2_path(T), nconf, nboot
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