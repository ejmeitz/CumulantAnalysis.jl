# CumulantAnalysis.jl

### Installation

This package depends on LatticeDynamicsToolkit.jl which is an unregistered package. To install in your chosen environment run:

```julia
using Pkg
Pkg.add(; url = "https://github.com/ejmeitz/LatticeDynamicsToolkit.jl.git", rev = "v0.1.0")
```

Then install CumulantAnalysis.jl into the same environment as LatticeDynamicsToolkit.jl

```julia
Pkg.add(; url = "https://github.com/ejmeitz/CumulantAnalysis.jl.git", rev = "v0.1.0")
```

### Citation

Coming soon.

### Solid Neon Example

Before we can compute thermodynamic properties we need to get the 2nd, 3rd and 4th order force constants from sTDEP. The method implemented by CumulantAnalysis.jl expects self-consistent phonons (e.g. sTDEP or SSCHA). If you use a method like MD-TDEP or a finite-dispalcement method your results will be less accurate. An in-depth sTDEP tutorial can be found [here](https://github.com/tdep-developers/tdep-tutorials/tree/main/02_sampling), but I also provide a script to compute the IFCs automatically. A more compelx workflow (used in the paper) which loops over multiple volume, temperature pairs can be found [here]("/workflows/neon_lattice_params_stdep.jl").

The results I obtained from this code may be found in the "paper/neon_example" folder. Note that you will not get exactly the same numbers as sTDEP is stochastic, but the IFCs should converge to roughly the same values.


```julia
using CumulantAnalysis.jl

repo_root = "<path-to-repo-root>" #UPDATE
outpath = "<whatever-directory-you-want>" # UPDATE

T = 24 # Kelvin

# Potential Definition (LAMMPS Commands)
r_cut = 6.955
pot_cmds = ["pair_style lj/cut $(r_cut)", "pair_coeff * * 0.0032135 2.782", "pair_modify shift yes"]

# sTDEP Parameters
n_iter = 10 # Number of self-consistent iterations
maximum_frequency = 2.5 # Frequency used to create initial IFC guess for iteration 0
quantum = true

# Must be POSCAR format
ssposcar_path = joinpath(repo_root, "paper", "neon_example", "ssposcar_T24")
ucposcar_path = joinpath(repo_root, "paper", "neon_example", "ucposcar_T24")

make_stdep_ifcs(
    ucposcar_path,
    ssposcar_path,
    outpath,
    pot_cmds,
    n_iter,
    r_cut,
    T,
    maximum_frequency,
    quantum
)
```

The force-constants for this example may also be found in the "paper/neon_example" folder of the repo. A full workflow (used in the paper) which loops over multiple temperatures can be found [here]("/workflows/neon.jl"). This script will create an output file for each thermodynamic property broken down into the harmonic, 0th, 1st and 2nd order corrections. Only the 0th order correction has an associated error. If `size_study` is set to `true` an additional output will contain the 0-th order correction as a function of the number of samples. This can be useful to detect convergence. 

```julia
using CumulantAnalysis

repo_root = "<path-to-repo-root>" #UPDATE
outpath = "<whatever-directory-you-want>" # UPDATE
basepath = joinpath(repo_root, "paper", "neon_example")

T = 24 # Kelvin

# Potential Definition (LAMMPS Commands)
r_cut = 6.955
pot_cmds = ["pair_style lj/cut $(r_cut)", "pair_coeff * * 0.0032135 2.782", "pair_modify shift yes"]

# Algorithm Parameters:
# Number of configurations sampled to estimate the 0-th order term
nconf = 100_000
# Number of bootstraps done to estimate error of 0-th order term
nboot = 5000
# Optional study to asses convergence of 0-th order term. Increses compute time.
size_study = true
# K-Mesh used to integrate harmonic properties, 30x30x30 is the default
harmonic_q_mesh = [30, 30, 30]
# K-Mesh sued to integrate 1st and 2nd order corrections, 25x25x25 is default. Typically a small grid is enough.
free_energy_q_mesh = [15, 15, 15]
# Other kwargs: n_threads, automatically uses all available threads.

# Must be POSCAR format
ssposcar_path = joinpath(basepath "ssposcar_T24")
ucposcar_path = joinpath(basepath "ucposcar_T24")

# Expects IFCs in TDEP format
ifc2_path = joinpath(basepath, "infile.forceconstant")
ifc3_path = joinpath(basepath, "infile.forceconstant_thirdorder")
ifc4_path = joinpath(basepath, "infile.forceconstant_fourthorder")

crystal_thermodynamic_properties(
    [T],
    outpath,
    ucposcar_path,
    ssposcar_path,
    ifc2_path,
    ifc3_path,
    ifc4_path,
    pot_cmds;
    quantum = true,
    nconf = nconf,
    nboot = nboot,
    size_study = size_study,
    harmonic_q_mesh = harmonic_q_mesh
    free_energy_q_mesh = free_energy_q_mesh
)

```

