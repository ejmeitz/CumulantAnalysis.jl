# [Installation](@id Installation)

CrystalCumulants.jl is written in Julia. A Python wrapper is also available via [juliacall](https://juliapy.github.io/PythonCall.jl/stable/juliacall/).

## Julia

CrystalCumulants.jl depends on [LatticeDynamicsToolkit.jl](https://github.com/ejmeitz/LatticeDynamicsToolkit.jl), which is unregistered. Linux is required; macOS may work; Windows is not supported.

### Requirements

- Julia 1.10+
- Linux (macOS may work; Windows is not supported)

!!! tip
    Set `JULIA_NUM_THREADS` in your environment **before** launching Julia. Thread count is fixed at startup and cannot be changed from within a running session.

    ```bash
    export JULIA_NUM_THREADS=40
    ```

### Install

```julia
using Pkg
Pkg.add(; url = "https://github.com/ejmeitz/LatticeDynamicsToolkit.jl.git", rev = "v0.1.2")
Pkg.add(; url = "https://github.com/ejmeitz/CrystalCumulants.jl.git", rev = "v0.1.1")
```

LAMMPS is installed automatically as a dependency. If a GPU is detected, a CUDA build of LAMMPS may be installed; it is not used at runtime but can cause compile issues on older systems. Open an issue if you run into LAMMPS build problems.

```julia
using CrystalCumulants
```

## Python

The [Python wrapper](https://github.com/ejmeitz/CrystalCumulants.jl/tree/main/python) calls the Julia package through juliacall. Julia dependencies (including LAMMPS and CrystalCumulants.jl) are managed automatically via `juliapkg.json`. First use may take a few minutes while the Julia environment is set up.

### Requirements

- Python 3.10+
- Linux (macOS may work; Windows is not supported)

### Environment variables

!!! warning
    1. Set `PYTHON_JULIACALL_HANDLE_SIGNALS=yes`, or Python cannot pass threads through to Julia. Ctrl-C may not stop the process; you may need to kill it manually.
    2. Set `PYTHON_JULIACALL_THREADS=<n-threads>` to control Julia thread count (default is 1). You may also need `JULIA_NUM_THREADS`.

### Install

From the repository root:

```bash
pip install -e ./python
```

```python
from cumulant_analysis import make_stdep_ifcs, crystal_thermodynamic_properties
```

See the [API Reference](@ref API) and [Example](@ref Example) pages to get started.
