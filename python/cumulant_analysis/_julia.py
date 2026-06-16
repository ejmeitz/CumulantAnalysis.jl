from __future__ import annotations

import os
import sys
import warnings
from collections.abc import Callable, Sequence
from typing import Any

_JL = None
_CA = None


def _configured_thread_setting() -> str | None:
    xopt = getattr(sys, "_xoptions", {}).get("juliacall-threads")
    if xopt is not None:
        return str(xopt)
    return os.environ.get("PYTHON_JULIACALL_THREADS")


def _report_julia_threads(jl) -> None:
    configured = _configured_thread_setting()
    nthreads = int(jl.Threads.nthreads())

    if configured is None:
        warnings.warn(
            f"cumulant_analysis: PYTHON_JULIACALL_THREADS is not set; "
            f"Julia is using {nthreads} thread(s). "
            "Set PYTHON_JULIACALL_THREADS or pass -X juliacall-threads=N when starting Python.",
            stacklevel=3,
        )
        return

    print(
        f"cumulant_analysis: PYTHON_JULIACALL_THREADS={configured!r}, "
        f"Julia is using {nthreads} thread(s)",
        file=sys.stderr,
    )


def get_jl():
    global _JL, _CA
    if _JL is not None:
        return _JL, _CA

    import juliacall

    jl = juliacall.newmodule("CumulantAnalysisPy")
    _report_julia_threads(jl)
    jl.seval("using CumulantAnalysis")
    jl.seval("""
    function _path_from_py(f)
        T -> begin
            using PythonCall
            pyconvert(String, f(T))
        end
    end
    """)

    _JL = jl
    _CA = jl.CumulantAnalysis
    return _JL, _CA


def to_julia_path(path: str | Callable[[float], str]):
    if isinstance(path, str):
        return path
    jl, _ = get_jl()
    return jl._path_from_py(path)


def to_julia_vector(values: Sequence[float | int]):
    jl, _ = get_jl()
    return jl.Vector(values)


def to_julia_pot_cmds(pot_cmds: str | Sequence[str]):
    if isinstance(pot_cmds, str):
        return pot_cmds
    jl, _ = get_jl()
    return jl.Vector([str(cmd) for cmd in pot_cmds])


def call_kwargs(**kwargs: Any) -> dict[str, Any]:
    return {key: value for key, value in kwargs.items() if value is not None}
