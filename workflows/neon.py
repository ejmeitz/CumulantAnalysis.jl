from os.path import join

from cumulant_analysis import crystal_thermodynamic_properties

nconf = 100_000
nboot = 5000

quantum = True
base_outpath = "/home/emeitz/Neon_ANALYTICAL_PIMD"
getoutpath = lambda T: join(base_outpath, f"T{T}")

Ts = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]
pot_cmds = [
    "pair_style lj/cut 6.955",
    "pair_coeff * * 0.0032135 2.782",
    "pair_modify shift yes",
]

## sTDEP IFCs
ifc_basepath = lambda T: f"/home/emeitz/Neon_sTDEP/results/T{T}"
ucposcar_path = lambda T: join(ifc_basepath(T), "infile.ucposcar")
ssposcar_path = lambda T: join(ifc_basepath(T), "infile.ssposcar")
ifc2_path = lambda T: join(ifc_basepath(T), "infile.forceconstant")
ifc3_path = lambda T: join(ifc_basepath(T), "infile.forceconstant_thirdorder")
ifc4_path = lambda T: join(ifc_basepath(T), "infile.forceconstant_fourthorder")

crystal_thermodynamic_properties(
    Ts,
    getoutpath,
    ucposcar_path,
    ssposcar_path,
    ifc2_path,
    ifc3_path,
    ifc4_path,
    pot_cmds,
    quantum=quantum,
    nconf=nconf,
    nboot=nboot,
    size_study=True,
)
