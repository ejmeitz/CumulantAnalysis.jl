# CumulantAnalysis.jl

### Installation

This package depends on LatticeDynamicsToolkit.jl which is an unregistered package and requires Linux (macOS might work, Windows will not). To install in your chosen environment run:

```julia
using Pkg
Pkg.add(; url = "https://github.com/ejmeitz/LatticeDynamicsToolkit.jl.git", rev = "v0.1.0")
```

Then install CumulantAnalysis.jl into the same environment as LatticeDynamicsToolkit.jl

```julia
Pkg.add(; url = "https://github.com/ejmeitz/CumulantAnalysis.jl.git", rev = "v0.1.0")
```

Note that this package automatically installs LAMMPS and if a GPU is detected it will install a GPU version of LAMMPS. If you have compilation errors related to this open a PR. The GPU is not used, but can still cause headaches at compile time if your Linux is too "old".

### Citation

Coming soon.

### Solid Neon Example

To begin clone the repo, it contains some input files.
```
git clone --depth 1 --branch v0.1.0 https://github.com/ejmeitz/CumulantAnalysis.jl.git
```

##### sTDEP IFCs
The first step is to get the 2nd, 3rd and 4th order force constants from sTDEP. Or skip to the [thermodynamic property](#thermodynamic-properties) section. The method implemented by CumulantAnalysis.jl expects self-consistent phonons (e.g. sTDEP or SSCHA). If you use a method like MD-TDEP or a finite-dispalcement method your results will be less accurate. An in-depth sTDEP tutorial can be found [here](https://github.com/tdep-developers/tdep-tutorials/tree/main/02_sampling), but I also provide a script to compute the IFCs automatically. A more compelx workflow (used in the paper) which loops over multiple volume, temperature pairs can be found [here](workflows/neon_lattice_params_stdep.jl).

The force constants obtained from this method can be found [here](data/stdep_results). Note that you will not get exactly the same numbers as sTDEP is stochastic, but the IFCs should converge to roughly the same values. The results can be found in the `iter009` folder, and take about a minute to obtain.


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
quantum = true # Whether to sample quantum or classical configs
# Other sTDEP kwargs: 
# - mix = true --> Determines whether configs from prior step are mixed with current step
# - nconf_init = 8 --> Number of configs used for iteration 0
# - max_configs = 512 --> Maximum number of configurations for any iteration (doubles every iter)

# Must be POSCAR format
basepath = joinpath(repo_root, "data", "stdep_inputs")
ssposcar_path = joinpath(basepath, "infile.ssposcar")
ucposcar_path = joinpath(basepath, "infile.ucposcar")

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

##### 3rd and 4th Order IFCs
This has only computed the second-order IFCs. Now you can use a separate TDEP installation or LatticeDynamicsToolkit.jl to get the 3rd and 4th order IFCs. This will take a couple of minutes depending how many threads you have.

```julia
import LatticeDynamicsToolkit.TDEPWrapper: execute, ExtractForceConstants

rc2 = 6.955
rc3 = 6.955
rc4 = 4.0

efc = ExtractForceConstants(secondorder_cutoff = rc2, thirdorder_cutoff = rc3, fourthorder_cutoff = rc4)
rundir = joinpath(outpath, "iter009") # Update if you changed n_iter above.
execute(efc, rundir, Threads.nthreads())
```

If you use the IFCs from this result, you might have to rename them from `outfile.*` to `infile.*` for the next part to work.

##### Thermodynamic Properties
The next step is to compute the thermodynamic properties. You can use the IFCs from above (make sure filenames match paths in this example) or use the IFCs [here](data/thermo_inputs). This script will create an output file for each thermodynamic property (F, U, S, Cv) broken down into the harmonic, 0th, 1st and 2nd order parts. Only the 0th order correction has an associated error. If `size_study` is set to `true` an additional output will contain the 0-th order correction as a function of the number of samples. This can be useful to detect convergence. A full workflow (used in the paper) which loops over multiple temperatures can be found [here](workflows/neon.jl).  

The free energy should be roughly -0.0179270 eV/atom. This value is stochastic, but the first few decimals should definitely match. On my computer this took ~3 minutes to run on 40 cores. The full set of expected results can be found [here](data/thermo_results). Note that these used a 25x25x25 mesh for the free energy.

```julia
using CumulantAnalysis

repo_root = "<path-to-repo-root>" #UPDATE
outpath = "<whatever-directory-you-want>" # UPDATE
basepath = joinpath(repo_root, "data", "thermo_inputs")

T = 24 # Kelvin

# Potential Definition (LAMMPS Commands)
r_cut = 6.955
pot_cmds = ["pair_style lj/cut $(r_cut)", "pair_coeff * * 0.0032135 2.782", "pair_modify shift yes"]

# Number of configurations sampled to estimate the 0-th order term
nconf = 100_000
# Number of bootstraps done to estimate error of 0-th order term
nboot = 5000
# Optional study to asses convergence of 0-th order term. Increses compute time.
size_study = true
# K-Mesh used to integrate harmonic properties, 30x30x30 is the default
harmonic_q_mesh = [30, 30, 30]
# K-Mesh sued to integrate 1st and 2nd order corrections, 25x25x25 is default. 
# Typically ~ 15x15x15 grid is sufficient.
free_energy_q_mesh = [15, 15, 15]
# Other kwargs: n_threads, automatically uses all available threads. Must set JULIA_NUM_THREADS env var prior to launching Julia.

# Must be POSCAR format
ssposcar_path = joinpath(basepath, "infile.ssposcar")
ucposcar_path = joinpath(basepath, "infile.ucposcar")

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
    harmonic_q_mesh = harmonic_q_mesh,
    free_energy_q_mesh = free_energy_q_mesh
)

```

