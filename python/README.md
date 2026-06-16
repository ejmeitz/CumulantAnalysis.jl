# cumulant-analysis

Python wrapper for [CumulantAnalysis.jl](https://github.com/ejmeitz/CumulantAnalysis.jl).

Julia dependencies are managed automatically by [juliacall](https://juliapy.github.io/PythonCall.jl/stable/juliacall/) via `juliapkg.json` in this package. On first use, juliacall will set up a Julia environment, install `CumulantAnalysis` and its dependencies (including LAMMPS). This may take a few minutes. If you have NVIDIA GPUs this might install a CUDA artifact. That is not used, but I have not figured out a great way to prevent this from happening automatically.

> [!IMPORTANT]
> 1) Be sure to set `PYTHON_JULIACALL_HANDLE_SIGNALS=yes` in your environment or Python will not be able to pass through the threads to Julia. This will prevent Ctrl-C from working to kill the process. You will likely have to manually kill the process to end it.
> 2) To specify the number of threads set `PYTHON_JULIACALL_THREADS=<n-threads>` in your environment. The default will be 1. You might also need to set `JULIA_NUM_THREADS` in your environment.

## Requirements

- Python 3.10+
- Linux (macOS may work; Windows is not supported)

## Installation

```bash
pip install -e ./python
```

## Usage

The `neon.jl` example is ported in the `neon.py` file. A small example for a single tempearture is shown below.

```python
from cumulant_analysis import make_stdep_ifcs, crystal_thermodynamic_properties

pot_cmds = [
    "pair_style lj/cut 6.955",
    "pair_coeff * * 0.0032135 2.782",
    "pair_modify shift yes",
]

make_stdep_ifcs(
    "data/stdep_inputs/infile.ucposcar",
    "data/stdep_inputs/infile.ssposcar",
    "out/stdep",
    pot_cmds,
    n_iter=10,
    r_cut=6.955,
    T=24.0,
    maximum_frequency=2.5,
    quantum=True,
)

crystal_thermodynamic_properties(
    [24.0],
    "out/thermo",
    "data/thermo_inputs/infile.ucposcar",
    "data/thermo_inputs/infile.ssposcar",
    "data/thermo_inputs/infile.forceconstant",
    "data/thermo_inputs/infile.forceconstant_thirdorder",
    "data/thermo_inputs/infile.forceconstant_fourthorder",
    pot_cmds,
    quantum=True,
    nconf=100_000,
    size_study=True,
)

# Temperature-dependent paths
base = "/data/stdep/RESULTS"
crystal_thermodynamic_properties(
    [100, 200, 300],
    lambda T: f"/out/T{int(T)}",
    lambda T: f"{base}/T{int(T)}_0/infile.ucposcar",
    lambda T: f"{base}/T{int(T)}_0/infile.ssposcar",
    lambda T: f"{base}/T{int(T)}_0/infile.forceconstant",
    lambda T: f"{base}/T{int(T)}_0/infile.forceconstant_thirdorder",
    lambda T: f"{base}/T{int(T)}_0/infile.forceconstant_fourthorder",
    pot_cmds,
)
```

Results are written to disk (for example `F_mean.txt`, `U_mean.txt`, `S_mean.txt`, `Cv_mean.txt`) rather than returned to Python.
