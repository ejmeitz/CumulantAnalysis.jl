from __future__ import annotations

import sys
from unittest.mock import MagicMock, patch

import pytest

import cumulant_analysis
from cumulant_analysis import _julia


def test_public_api_exports():
    assert hasattr(cumulant_analysis, "make_stdep_ifcs")
    assert hasattr(cumulant_analysis, "crystal_thermodynamic_properties")


def test_to_julia_path_string():
    assert _julia.to_julia_path("/tmp/out") == "/tmp/out"


def test_to_julia_path_callable():
    mock_jl = MagicMock()
    mock_jl._path_from_py.return_value = "julia-fn"
    with patch.object(_julia, "get_jl", return_value=(mock_jl, MagicMock())):
        result = _julia.to_julia_path(lambda T: f"/out/T{int(T)}")
    assert result == "julia-fn"
    mock_jl._path_from_py.assert_called_once()


def test_to_julia_pot_cmds_string():
    assert _julia.to_julia_pot_cmds("pair_style lj/cut 6.955") == "pair_style lj/cut 6.955"


def test_call_kwargs_drops_none():
    assert _julia.call_kwargs(a=1, b=None) == {"a": 1}


def test_lazy_julia_init():
    assert _julia._JL is None
    assert _julia._CA is None


def test_configured_thread_setting_from_env(monkeypatch):
    monkeypatch.delattr(sys, "_xoptions", raising=False)
    monkeypatch.setenv("PYTHON_JULIACALL_THREADS", "8")
    assert _julia._configured_thread_setting() == "8"


def test_configured_thread_setting_from_xoption(monkeypatch):
    monkeypatch.setattr(sys, "_xoptions", {"juliacall-threads": "4"})
    monkeypatch.setenv("PYTHON_JULIACALL_THREADS", "8")
    assert _julia._configured_thread_setting() == "4"


def test_report_julia_threads_warns_when_unset(monkeypatch):
    monkeypatch.delattr(sys, "_xoptions", raising=False)
    monkeypatch.delenv("PYTHON_JULIACALL_THREADS", raising=False)
    mock_jl = MagicMock()
    mock_jl.Threads.nthreads.return_value = 1
    with pytest.warns(UserWarning, match="PYTHON_JULIACALL_THREADS is not set"):
        _julia._report_julia_threads(mock_jl)


def test_report_julia_threads_prints_when_set(monkeypatch, capsys):
    monkeypatch.delattr(sys, "_xoptions", raising=False)
    monkeypatch.setenv("PYTHON_JULIACALL_THREADS", "8")
    mock_jl = MagicMock()
    mock_jl.Threads.nthreads.return_value = 8
    _julia._report_julia_threads(mock_jl)
    captured = capsys.readouterr()
    assert "PYTHON_JULIACALL_THREADS='8'" in captured.err
    assert "Julia is using 8 thread(s)" in captured.err
