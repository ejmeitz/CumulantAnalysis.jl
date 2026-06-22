# CrystalCumulants.jl

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://ejmeitz.github.io/CrystalCumulants.jl)

A fast implementation of the free energy cumulant expansion for crystals. The code is written in Julia, but we also provide a Python wrapper. There are two entry points to the code:

- `make_stdep_ifcs` — computes sTDEP IFCs (2nd through 4th order) for a specific temperature
- `crystal_thermodynamic_properties` — using the sTDEP IFCs, computes quantum-anharmonic thermodynamic properties via the free energy cumulant expansion

> [!NOTE]
> The energy from polar interactions is not accounted for. Even if this interaction is present in the `infile.forceconstant` file it will be ignored. The cumulant expansion theory is easily modified to incorporate the polar contribution, but the corresponding code was not implemented or tested.

> [!TIP]
> 1. Be sure to set `JULIA_NUM_THREADS` or `PYTHON_JULIACALL_THREADS` in your environment to enable multi-threading of the code!
> 2. Always use the primitive cell. The number of atoms in the primitive cell dictates the computational cost and RAM usage (lower is better).
> 3. The `free_energy_q_mesh` and `nconf` dictate runtime for a given primitive cell. It is recommended to run a convergence study to assess what grid size and how many samples are needed to get converged results and minimize runtime.

Full documentation: [https://ejmeitz.github.io/CrystalCumulants.jl](https://ejmeitz.github.io/CrystalCumulants.jl)

## Installation

### Julia

CrystalCumulants.jl depends on [LatticeDynamicsToolkit.jl](https://github.com/ejmeitz/LatticeDynamicsToolkit.jl), which is unregistered. Linux is required; macOS may work; Windows is not supported.

**Requirements:** Julia 1.10+

```julia
using Pkg
Pkg.add(; url = "https://github.com/ejmeitz/LatticeDynamicsToolkit.jl.git", rev = "v0.1.2")
Pkg.add(; url = "https://github.com/ejmeitz/CrystalCumulants.jl.git", rev = "v0.1.1")
```

LAMMPS is installed automatically as a dependency. If a GPU is detected, a CUDA build of LAMMPS may be installed; it is not used at runtime but can cause compile issues on older systems. Open an issue if you run into LAMMPS build problems.

```julia
using CrystalCumulants
```

### Python

The [Python wrapper](python/) calls the Julia package through [juliacall](https://juliapy.github.io/PythonCall.jl/stable/juliacall/). Julia dependencies (including LAMMPS and CrystalCumulants.jl) are managed automatically via `juliapkg.json`. First use may take a few minutes while the Julia environment is set up.

**Requirements:** Python 3.10+, Linux (macOS may work; Windows not supported)

> [!WARNING]
> 1. Set `PYTHON_JULIACALL_HANDLE_SIGNALS=yes`, or Python cannot pass threads through to Julia. Ctrl-C may not stop the process; you may need to kill it manually.
> 2. Set `PYTHON_JULIACALL_THREADS=<n-threads>` to control Julia thread count (default is 1). You may also need `JULIA_NUM_THREADS`.

From the repository root:

```bash
pip install -e ./python
```

```python
from cumulant_analysis import make_stdep_ifcs, crystal_thermodynamic_properties
```

## Example

Worked examples for solid neon (sTDEP force constants and thermodynamic properties) in Julia and Python are in the documentation:

- [Example walkthrough](https://ejmeitz.github.io/CrystalCumulants.jl/stable/example/) — step-by-step workflow with bundled input files
- [API reference](https://ejmeitz.github.io/CrystalCumulants.jl/stable/api/) — argument lists and temperature-dependent paths

Clone the repository for bundled input data:

```bash
git clone --depth 1 --branch v0.1.1 https://github.com/ejmeitz/CrystalCumulants.jl.git
```

Additional scripts are in [`workflows/`](workflows/).

## Citation

Coming soon.
