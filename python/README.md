# cumulant-analysis

Python wrapper for [CumulantAnalysis.jl](https://github.com/ejmeitz/CumulantAnalysis.jl).

Julia dependencies are managed automatically by [juliacall](https://juliapy.github.io/PythonCall.jl/stable/juliacall/) via `juliapkg.json` in this package. On first use, juliacall will set up a Julia environment, install `CumulantAnalysis` and its dependencies (including LAMMPS), and may take a few minutes.

## Requirements

- Python 3.10+
- Linux (macOS may work; Windows is not supported — same constraint as the Julia package)

## Installation

From this repository:

```bash
pip install git+https://github.com/ejmeitz/CumulantAnalysis.jl.git#subdirectory=python
```

For local development:

```bash
pip install -e ./python
```

## Threading

`cumulant_analysis` uses Julia threads for parallel work. Configure the thread count **before starting Python** — it cannot be changed after Julia loads.

Either set the environment variable:

```bash
export PYTHON_JULIACALL_THREADS=8
python your_script.py
```

or pass the `-X` flag to the Python interpreter:

```bash
python -X juliacall-threads=8 your_script.py
```

If `-X juliacall-threads` is set, it takes precedence over `PYTHON_JULIACALL_THREADS`. Use `auto` to match the number of CPU cores.

On the first call that starts Julia, the package reports the configured setting and the active thread count to stderr. If `PYTHON_JULIACALL_THREADS` is unset (and no `-X juliacall-threads` flag was passed), a warning is emitted showing how many threads Julia is using.

You can also pass `n_threads` to `crystal_thermodynamic_properties`, but that only limits parallelism inside the calculation; the Julia process thread pool is still set at startup.

See the [juliacall configuration docs](https://juliapy.github.io/PythonCall.jl/stable/juliacall/#Configuration) for details.

## Usage

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
