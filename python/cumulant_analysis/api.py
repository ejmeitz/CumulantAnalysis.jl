from __future__ import annotations

from collections.abc import Callable, Sequence

from cumulant_analysis._julia import (
    call_kwargs,
    get_jl,
    to_julia_int_vector,
    to_julia_path,
    to_julia_pot_cmds,
    to_julia_vector,
)


def crystal_thermodynamic_properties(
    temperatures: Sequence[float],
    outpath: str | Callable[[float], str],
    ucposcar_path: str | Callable[[float], str],
    ssposcar_path: str | Callable[[float], str],
    ifc2_path: str | Callable[[float], str],
    ifc3_path: str | Callable[[float], str],
    ifc4_path: str | Callable[[float], str],
    pot_cmds: str | list[str],
    *,
    quantum: bool = False,
    nconf: int = 100_000,
    nboot: int = 2500,
    size_study: bool = False,
    harmonic_q_mesh: Sequence[int] = (30, 30, 30),
    free_energy_q_mesh: Sequence[int] = (25, 25, 25),
    n_threads: int | None = None,
) -> None:
    """Calculate thermodynamic properties across a temperature range.

    All path-like parameters may be passed as a string, or as a callable that
    takes a temperature (Kelvin) and returns a path string.

    Parameters
    ----------
    temperatures:
        Temperatures (Kelvin) at which to calculate thermodynamic properties.
    outpath:
        Output directory for results at each temperature.
    ucposcar_path:
        Path to the unit-cell POSCAR file (TDEP format).
    ssposcar_path:
        Path to the supercell POSCAR file (TDEP format).
    ifc2_path:
        Path to the second-order force constant file (TDEP format).
    ifc3_path:
        Path to the third-order force constant file (TDEP format).
    ifc4_path:
        Path to the fourth-order force constant file (TDEP format).
    pot_cmds:
        LAMMPS potential commands used to evaluate the true potential energy.
    quantum:
        Whether to use quantum statistics.
    nconf:
        Number of configurations sampled for the constant correction.
    nboot:
        Number of bootstrap samples for error estimation of the constant correction.
    size_study:
        Whether to compute the constant correction as a function of sample count.
    harmonic_q_mesh:
        q-mesh used for the harmonic contribution.
    free_energy_q_mesh:
        q-mesh used for cumulant corrections.
    n_threads:
        Number of Julia threads to use. Defaults to all available Julia threads.

    Notes
    -----
    Results are written to disk (for example ``F_mean.txt``, ``U_mean.txt``,
    ``S_mean.txt``, and ``Cv_mean.txt``) rather than returned to Python.
    """
    _, ca = get_jl()

    kwargs = {
        "quantum": quantum,
        "nconf": nconf,
        "nboot": nboot,
        "size_study": size_study,
        "harmonic_q_mesh": to_julia_int_vector(harmonic_q_mesh),
        "free_energy_q_mesh": to_julia_int_vector(free_energy_q_mesh),
    }
    if n_threads is not None:
        kwargs["n_threads"] = n_threads

    ca.crystal_thermodynamic_properties(
        to_julia_vector(temperatures),
        to_julia_path(outpath),
        to_julia_path(ucposcar_path),
        to_julia_path(ssposcar_path),
        to_julia_path(ifc2_path),
        to_julia_path(ifc3_path),
        to_julia_path(ifc4_path),
        to_julia_pot_cmds(pot_cmds),
        **kwargs,
    )


def make_stdep_ifcs(
    ucposcar_path: str,
    ssposcar_path: str,
    outdir: str,
    pot_cmds: list[str],
    n_iter: int,
    r_cut: float,
    T: float,
    maximum_frequency: float,
    quantum: bool,
    **kwargs,
) -> None:
    """Generate self-consistent TDEP force constants with sTDEP.

    Parameters
    ----------
    ucposcar_path:
        Path to the unit-cell POSCAR file.
    ssposcar_path:
        Path to the supercell POSCAR file.
    outdir:
        Output directory for sTDEP results.
    pot_cmds:
        LAMMPS potential commands.
    n_iter:
        Number of self-consistent sTDEP iterations.
    r_cut:
        Pair potential cutoff used by sTDEP.
    T:
        Temperature (Kelvin).
    maximum_frequency:
        Maximum frequency used for the initial IFC guess.
    quantum:
        Whether to sample quantum or classical configurations.
    **kwargs:
        Additional keyword arguments forwarded to the underlying sTDEP call.

    Notes
    -----
    Output files are written under ``outdir``.
    """
    _, ca = get_jl()

    ca.make_stdep_ifcs(
        ucposcar_path,
        ssposcar_path,
        outdir,
        to_julia_pot_cmds(pot_cmds),
        n_iter,
        r_cut,
        T,
        maximum_frequency,
        quantum,
        **call_kwargs(**kwargs),
    )
