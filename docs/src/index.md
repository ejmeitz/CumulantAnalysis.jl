# [CrystalCumulants.jl](@id Home)

A fast implementation of the free energy cumulant expansion for crystals. The code is written in Julia, but we also provide a Python wrapper. There are two entry points to the code:

- `make_stdep_ifcs` — computes sTDEP IFCs (2nd through 4th order) for a specific temperature
- `crystal_thermodynamic_properties` — using the sTDEP IFCs, computes quantum-anharmonic thermodynamic properties via the free energy cumulant expansion

## Getting started

- [Installation](@ref Installation) — Julia and Python setup
- [API Reference](@ref API) — function arguments, keywords, and temperature-dependent paths
- [Theory](@ref Theory) — background on the cumulant expansion and how it maps to the implementation
- [Example](@ref Example) — solid Neon walkthrough (Julia and Python)


!!! note
    The energy from polar interactions is not accounted for. Even if this interaction is present in the `infile.forceconstant` file it will be ignored. The cumulant expansion theory is easily modified to incorporate the polar contribution, but the corresponding code was not implemented or tested.

!!! tip
    1. Be sure to set `JULIA_NUM_THREADS` or `PYTHON_JULIACALL_THREADS` in your environment to enable multi-threading of the code!
    2. Always use the primitive cell. The number of atoms in the primitive cell dictates the computational cost and RAM usage (lower is better).
    3. The `free_energy_q_mesh` and `nconf` dictate runtime for a given primitive cell. It is recommended to run a convergence study to assess what grid size and how many samples are needed to get converged results and minimize runtime.

## Citation

Coming soon.