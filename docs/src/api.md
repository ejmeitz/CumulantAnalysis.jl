# [API Reference](@id API)

CrystalCumulants.jl exposes two functions. The Julia and Python interfaces use the same argument names and semantics; only syntax differs (keyword vs `**kwargs`, `true` vs `True`, etc.).

## [`make_stdep_ifcs`](@id make_stdep_ifcs)

Generate self-consistent TDEP force constants. Second-order IFCs are iterated to convergence; 3rd- and 4th-order IFCs are fit from the final configuration set. Harmonic free energy per iteration is written to `harmonic_free_energy.txt` to aid in assesing convergence. DOS and dispersion are also calculated at each iteration.

**Arguments** (positional, in order):

| Argument | Description |
|----------|-------------|
| `ucposcar_path` | Path to unit-cell POSCAR (TDEP format) |
| `ssposcar_path` | Path to supercell POSCAR (TDEP format) |
| `outdir` | Output directory for sTDEP results |
| `pot_cmds` | LAMMPS potential commands (`Vector{String}` / `list[str]`) |
| `n_iter` | Number of self-consistent sTDEP iterations |
| `r_cut` | Pair potential cutoff for 2nd-order IFCs |
| `T` | Temperature (K) |
| `maximum_frequency` | Maximum frequency for the initial IFC guess |
| `quantum` | Sample quantum or classical configurations |
| `rc3` | Cutoff radius for 3rd-order IFC fitting |
| `rc4` | Cutoff radius for 4th-order IFC fitting |

**Optional keyword arguments** are forwarded to sTDEP (Julia: `kwargs...`; Python: `**kwargs`). Common options include `mix`, `nconf_init`, and `max_configs`.

Output files (2nd- through 4th-order IFCs and related sTDEP outputs) are written under `outdir`.

## [`crystal_thermodynamic_properties`](@id crystal_thermodynamic_properties)

Calculate thermodynamic properties across a temperature range using the free energy cumulant expansion.

**Arguments** (positional, in order):

| Argument | Description |
|----------|-------------|
| `temperatures` | Temperatures to evaluate (K), e.g. `[24]` or `[100, 200, 300]` |
| `outpath` | Output directory (or temperature-dependent path; see below) |
| `ucposcar_path` | Unit-cell POSCAR (TDEP format) |
| `ssposcar_path` | Supercell POSCAR (TDEP format) |
| `ifc2_path` | 2nd-order IFC file (`infile.forceconstant`) |
| `ifc3_path` | 3rd-order IFC file (`infile.forceconstant_thirdorder`) |
| `ifc4_path` | 4th-order IFC file (`infile.forceconstant_fourthorder`) |
| `pot_cmds` | LAMMPS potential commands for the true potential energy |

**Keyword arguments:**

| Keyword | Default | Description |
|---------|---------|-------------|
| `quantum` | `false` | Use quantum (`true`) or classical (`false`) statistics |
| `nconf` | `100_000` | Configurations sampled for the 0th-order correction |
| `nboot` | `2500` | Bootstrap samples for 0th-order standard error |
| `size_study` | `false` | If enabled, write 0th-order correction vs. sample count |
| `harmonic_q_mesh` | `[30, 30, 30]` | q-mesh for the harmonic contribution |
| `free_energy_q_mesh` | `[25, 25, 25]` | q-mesh for 1st- and 2nd-order corrections |
| `n_threads` | all Julia threads | Parallel thread count (be sure to set `JULIA_NUM_THREADS` / `PYTHON_JULIACALL_THREADS`) |

### Output files

For each property, the code writes `F_mean.txt`, `U_mean.txt`, `S_mean.txt`, and `Cv_mean.txt` with rows broken down into harmonic (`F0`, etc.), 0th-order offset, 1st-order, 2nd-order, and total. Only the 0th-order correction has an associated bootstrap standard error. HDF5 versions (`.h5`) are also written.

## Temperature-dependent paths

For `crystal_thermodynamic_properties`, the path arguments (`outpath`, `ucposcar_path`, `ssposcar_path`, `ifc2_path`, `ifc3_path`, `ifc4_path`) may be either:

- a fixed path string, when inputs are the same for every temperature, or
- a function of temperature `T` (in K), when structures or IFCs vary with `T`.

In Julia, pass an anonymous function `(T) -> path`. In Python, pass a `lambda`.

**Julia**

```julia
using CrystalCumulants

base = "/data/stdep/RESULTS"
temperatures = [100.0, 200.0, 300.0]

crystal_thermodynamic_properties(
    temperatures,
    (T) -> "/out/T$(Int(T))",
    (T) -> joinpath(base, "T$(Int(T))_0", "infile.ucposcar"),
    (T) -> joinpath(base, "T$(Int(T))_0", "infile.ssposcar"),
    (T) -> joinpath(base, "T$(Int(T))_0", "infile.forceconstant"),
    (T) -> joinpath(base, "T$(Int(T))_0", "infile.forceconstant_thirdorder"),
    (T) -> joinpath(base, "T$(Int(T))_0", "infile.forceconstant_fourthorder"),
    pot_cmds,
)
```

**Python**

```python
from os.path import join

from cumulant_analysis import crystal_thermodynamic_properties

base = "/data/stdep/RESULTS"
temperatures = [100.0, 200.0, 300.0]

crystal_thermodynamic_properties(
    temperatures,
    lambda T: f"/out/T{int(T)}",
    lambda T: join(base, f"T{int(T)}_0", "infile.ucposcar"),
    lambda T: join(base, f"T{int(T)}_0", "infile.ssposcar"),
    lambda T: join(base, f"T{int(T)}_0", "infile.forceconstant"),
    lambda T: join(base, f"T{int(T)}_0", "infile.forceconstant_thirdorder"),
    lambda T: join(base, f"T{int(T)}_0", "infile.forceconstant_fourthorder"),
    pot_cmds,
)
```

Each function receives the temperature `T` for the current evaluation and must return a path string.
