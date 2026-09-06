"""Hardware probe, run configuration and the command-line prompt."""

import argparse
import io

import numpy as np
import pytest

from aa540fem import cli, hardware
from aa540fem.hardware import HardwareInfo, MPIInfo, RunConfig, configure, probe, recommend


def fake_hardware(cores=8, ram_gb=32.0, gpu=False, mpi=True):
    gpus = [hardware.GPUInfo(0, "Test GPU", 16 * 1024 ** 3)] if gpu else []
    return HardwareInfo(physical_cores=cores, logical_cores=2 * cores, usable_cores=cores,
                        ram_total=int(ram_gb * 1024 ** 3), ram_available=int(ram_gb * 1024 ** 3),
                        memory_limit=None, gpus=gpus,
                        mpi=MPIInfo(mpi, "mpirun" if mpi else None, 1),
                        backends={"numba": True, "petsc4py": False, "mumps": False,
                                  "mpi4py": mpi, "cuda": False}, blas="openblas")


def test_probe_never_raises_and_is_cached():
    hw = probe()
    assert hw.physical_cores >= 1 and hw.usable_cores >= 1 and hw.ram_total > 0
    assert probe() is hw and probe(refresh=True) is not hw
    assert "cores" in hw.describe()


def test_default_config_scales_with_the_fraction():
    hw = fake_hardware(cores=8, ram_gb=32)
    full = hardware.default_config(hw, 1.0)
    half = hardware.default_config(hw, 0.5)
    assert full.threads == 8 and half.threads == 4
    assert half.memory_budget == pytest.approx(0.5 * 0.8 * 32 * 1024 ** 3)
    assert full.memory_budget == pytest.approx(0.8 * 32 * 1024 ** 3)
    assert "50 %" in half.describe()


def test_configure_precedence(monkeypatch):
    hw = fake_hardware(cores=8)
    monkeypatch.setenv("AA540FEM_MACHINE_FRACTION", "50")
    monkeypatch.setenv("AA540FEM_ASSEMBLY", "numpy")
    cfg = configure(hw=hw)
    assert cfg.threads == 4 and cfg.assembly_backend == "numpy"
    cfg = configure(fraction=1.0, threads=3, assembly="numba", direct="superlu",
                    memory_budget=2.0, hw=hw)
    assert cfg.threads == 3 and cfg.assembly_backend == "numba"
    assert cfg.direct_backend == "superlu" and cfg.memory_budget == 2 * 1024 ** 3
    assert hardware.get_config() is cfg
    monkeypatch.delenv("AA540FEM_MACHINE_FRACTION")
    monkeypatch.delenv("AA540FEM_ASSEMBLY")
    configure(hw=hw)


def test_memory_check_warns_and_refuses():
    cfg = RunConfig(fraction=1.0, threads=1, memory_budget=10 * 1024 ** 2)
    import scipy.sparse as sp

    small = sp.eye(100, format="csr")
    assert cfg.check_memory(small) < cfg.memory_budget
    with pytest.raises(MemoryError):
        cfg.check_memory(10 ** 7, n=10 ** 6)
    assert cfg.check_memory(10 ** 7, n=10 ** 6, force=True) > cfg.memory_budget
    cfg.memory_budget = int(1.1 * cfg.check_memory(10 ** 7, n=10 ** 6, force=True))
    with pytest.warns(UserWarning):
        cfg.check_memory(10 ** 7, n=10 ** 6)


def parse(argv):
    parser = argparse.ArgumentParser()
    cli.add_run_arguments(parser)
    return parser.parse_args(argv)


def test_cli_flags_and_no_prompt(monkeypatch):
    monkeypatch.setattr(hardware, "probe", lambda refresh=False: fake_hardware(cores=8))
    monkeypatch.setattr(cli, "probe", hardware.probe)
    monkeypatch.delenv("AA540FEM_MACHINE_FRACTION", raising=False)
    cfg = cli.configure_from_args(parse(["--machine", "50", "--direct", "superlu"]),
                                  interactive=False)
    assert cfg.threads == 4 and cfg.fraction == 0.5 and cfg.direct_backend == "superlu"
    cfg = cli.configure_from_args(parse(["--threads", "2", "--no-prompt"]))
    assert cfg.threads == 2 and cfg.fraction == 1.0

    def no_input():
        raise AssertionError("stdin must not be read")
    cli.configure_from_args(parse(["--no-prompt"]), input_fn=no_input)


def test_cli_prompt(monkeypatch):
    monkeypatch.setattr(hardware, "probe", lambda refresh=False: fake_hardware(cores=8))
    monkeypatch.setattr(cli, "probe", hardware.probe)
    monkeypatch.delenv("AA540FEM_MACHINE_FRACTION", raising=False)
    out = io.StringIO()
    cfg = cli.configure_from_args(parse([]), interactive=True, stream=out, input_fn=lambda: "1")
    assert cfg.fraction == 0.5 and cfg.threads == 4
    assert "50 %" in out.getvalue() and "100 %" in out.getvalue()
    cfg = cli.configure_from_args(parse([]), interactive=True, stream=io.StringIO(),
                                  input_fn=lambda: "")
    assert cfg.fraction == 1.0 and cfg.threads == 8


def test_recommendations():
    text = recommend(fake_hardware(gpu=True), n_unknowns=5 * 10 ** 5)
    assert "mpirun" in text and "GPU" in text
    assert "SuperLU" in recommend(fake_hardware(mpi=False))


def test_apply_sets_thread_counts(monkeypatch):
    numba = pytest.importorskip("numba")
    RunConfig(fraction=1.0, threads=2, memory_budget=10 ** 9).apply()
    assert numba.get_num_threads() == min(2, numba.config.NUMBA_NUM_THREADS)
    assert np.zeros(2).size == 2
    RunConfig(fraction=1.0, threads=numba.config.NUMBA_NUM_THREADS, memory_budget=10 ** 9).apply()
