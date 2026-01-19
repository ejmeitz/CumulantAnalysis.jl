using CumulantAnalysis


nconf = 100_000 
nboot = 5000

quantum = true
getoutpath = (T) -> joinpath(base_outpath,  "T$(T)")

# LENNARD JONES
Ts = [4,6,8,10,12,14,16,18,20,22,24]
ucposcar_path = (T) -> "/mnt/merged/emeitz/Neon_sTDEP/results/T$(T)/infile.ucposcar"
ssposcar_path = (T) -> "/mnt/merged/emeitz/Neon_sTDEP/results/T$(T)/infile.ssposcar"
ifc2_path = (T) -> "/mnt/merged/emeitz/Neon_sTDEP/IFCs/T$(T)/infile.forceconstant"
ifc3_path = (T) -> "/mnt/merged/emeitz/Neon_sTDEP/IFCs/T$(T)/infile.forceconstant_thirdorder"
ifc4_path = (T) -> "/mnt/merged/emeitz/Neon_sTDEP/IFCs/T$(T)/infile.forceconstant_fourthorder"
pot_cmds = ["pair_style lj/cut 6.955", "pair_coeff * * 0.0032135 2.782", "pair_modify shift yes"]


base_outpath = "/mnt/merged/emeitz/CumulantAnalysisTest/Neon_ANALYTICAL"


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
            ucposcar_path = ucposcar_path(T),
            ssposcar_path = ssposcar_path(T),
            size_study = true,
            quantum = quantum
        )
end