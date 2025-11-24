using CumulantAnalysis

# Run for a super long time, to try and compare to analytical
# solution implemented in TDEP
order = 2
nconf = 3_000_000 
nboot = 1000

base_outpath = raw"Z:\emeitz\Data\FreeEnergies\Sampled\LONG"
getoutpath = (pot, estim, T) -> joinpath(base_outpath, pot, estim, "T$(T)")

# LENNARD JONES
# pot = "LJ"
# Ts = [80]
# ucposcar_path = raw"C:\Users\ejmei\repos\TDEP_IFCs.jl\data\LJ\infile.ucposcar"
# ssposcar_path = raw"C:\Users\ejmei\repos\TDEP_IFCs.jl\data\LJ\infile.ssposcar" #4UC
# ifc2_path = (T) -> "Z:/emeitz/Data/ForceConstants/LJ_IFC_INTERPOLATION_NODES/IFCs/T$(T)/infile.forceconstant"
# ifc3_path = (T) -> "Z:/emeitz/Data/ForceConstants/LJ_IFC_INTERPOLATION_NODES/IFCs/T$(T)/infile.forceconstasnt_thirdorder"
# ifc4_path = (T) -> "Z:/emeitz/Data/ForceConstants/LJ_IFC_INTERPOLATION_NODES/IFCs/T$(T)/infile.forceconstasnt_fourthorder"
# pot_cmds = ["pair_style lj/cut 8.5", "pair_coeff * * 0.010423 3.4", "pair_modify shift yes"]


# STILLINGER WEBER
pot = "SW"
Ts = [1300]
ucposcar_path = raw"C:\Users\ejmei\repos\TDEP_IFCs.jl\data\SW\infile.ucposcar"
ssposcar_path = raw"C:\Users\ejmei\repos\TDEP_IFCs.jl\data\SW\infile.ssposcar" #3UC
ifc2_path = (T) -> "Z:/emeitz/Data/ForceConstants/SW_IFC_INTERPOLATION_NODES/IFCs/T$(T)/infile.forceconstant"
ifc3_path = (T) -> "Z:/emeitz/Data/ForceConstants/SW_IFC_INTERPOLATION_NODES/IFCs/T$(T)/infile.forceconstasnt_thirdorder"
ifc4_path = (T) -> "Z:/emeitz/Data/ForceConstants/SW_IFC_INTERPOLATION_NODES/IFCs/T$(T)/infile.forceconstasnt_fourthorder"

sw_pot = raw"C:\Users\ejmei\repos\TDEP_IFCs.jl\data\SW\Si.sw"
pot_cmds = ["pair_style sw", "pair_coeff * * \"$(sw_pot)\" Si"]

estims = ["V3+V4", "V-V2"]

for T in Ts
    for estim in estims
        @info "T = $(T), estim = $(estim)"

        o = getoutpath(pot, estim, T)
        mkpath(o)

        if estim == "V3+V4"
            e = FourthOrderEstimator(
                order, ifc2_path(T), ifc3_path(T), ifc4_path(T), nconf, nboot
            )
        elseif estim == "V-V2"
            e = HarmonicEstimator(
                order, ifc2_path(T), nconf, nboot
            )
        else 
            error("Unknown estimator $(estim)")
        end

        estimate(
                e,
                T, 
                o,
                pot_cmds;
                ucposcar_path = ucposcar_path,
                ssposcar_path = ssposcar_path,
                size_study = true
                # use_control_variates::Bool = false,
            )
    end
end