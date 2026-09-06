"""Hardware probe and run configuration.

:func:`probe` collects what the machine offers (cores, memory, NVIDIA GPUs,
MPI, installed backends); :class:`RunConfig` says how much of it a run may
use (a fraction of the cores, a memory budget, the assembly and solver
backends) and applies the thread counts to numba, the BLAS and MUMPS.
:func:`configure` builds the configuration from explicit arguments,
environment variables and the probe, never asking anything; the interactive
"50 % or 100 %" prompt of the examples lives in :mod:`aa540fem.cli`.

Environment variables: ``AA540FEM_MACHINE_FRACTION`` (0-1 or 50/100),
``AA540FEM_THREADS``, ``AA540FEM_ASSEMBLY`` (numpy/numba),
``AA540FEM_DIRECT`` (superlu/petsc), ``AA540FEM_LINEAR`` (scipy/petsc/
petsc-cuda), ``AA540FEM_MEMORY_BUDGET`` (GB), ``AA540FEM_NO_PROMPT``.
"""

from __future__ import annotations

import importlib.util
import os
import platform
import shutil
import subprocess
import warnings
from dataclasses import dataclass, field

# -- what the machine offers ------------------------------------------


@dataclass
class GPUInfo:
    index: int
    name: str
    memory_bytes: int


@dataclass
class MPIInfo:
    mpi4py: bool
    launcher: str | None
    world_size: int          # > 1 when already running under mpirun


@dataclass
class HardwareInfo:
    physical_cores: int
    logical_cores: int
    usable_cores: int        # affinity mask / cgroup quota
    ram_total: int
    ram_available: int
    memory_limit: int | None  # cgroup limit when set
    gpus: list[GPUInfo]
    mpi: MPIInfo
    backends: dict[str, bool]
    blas: str
    system: str = field(default_factory=platform.platform)

    def describe(self) -> str:
        gb = 1024 ** 3
        gpu = (", ".join(f"{g.name} ({g.memory_bytes / gb:.0f} GB)" for g in self.gpus)
               or "no NVIDIA GPU")
        mpi = (f"MPI: {self.mpi.launcher}" if self.mpi.launcher and self.mpi.mpi4py
               else "MPI: not available")
        if self.mpi.world_size > 1:
            mpi += f" ({self.mpi.world_size} ranks)"
        have = [k for k, v in self.backends.items() if v]
        return (f"{self.physical_cores} physical cores ({self.usable_cores} usable), "
                f"{self.ram_total / gb:.1f} GB RAM ({self.ram_available / gb:.1f} GB free), "
                f"{gpu}, {mpi}; installed: {', '.join(have) or 'numpy/scipy only'}")


def _read(path):
    try:
        with open(path) as f:
            return f.read()
    except OSError:
        return ""


def _physical_cores(logical):
    try:
        import psutil

        n = psutil.cpu_count(logical=False)
        if n:
            return int(n)
    except ImportError:
        pass
    if platform.system() == "Linux":
        pairs = set()
        phys = core = None
        for line in _read("/proc/cpuinfo").splitlines():
            if line.startswith("physical id"):
                phys = line.split(":")[1].strip()
            elif line.startswith("core id"):
                core = line.split(":")[1].strip()
                pairs.add((phys, core))
        if pairs:
            return len(pairs)
    elif platform.system() == "Darwin":
        try:
            out = subprocess.run(["sysctl", "-n", "hw.physicalcpu"], capture_output=True,
                                 text=True, timeout=2)
            return int(out.stdout.strip())
        except (OSError, ValueError, subprocess.SubprocessError):
            pass
    return logical


def _usable_cores(logical):
    n = logical
    try:
        n = min(n, len(os.sched_getaffinity(0)))
    except AttributeError:
        pass
    quota = _read("/sys/fs/cgroup/cpu.max").split()
    if len(quota) == 2 and quota[0] != "max":
        try:
            n = min(n, max(1, int(round(int(quota[0]) / int(quota[1])))))
        except ValueError:
            pass
    return max(1, n)


def _memory():
    total = available = None
    try:
        import psutil

        vm = psutil.virtual_memory()
        total, available = int(vm.total), int(vm.available)
    except ImportError:
        info = {}
        for line in _read("/proc/meminfo").splitlines():
            parts = line.split()
            if len(parts) >= 2 and parts[0].endswith(":"):
                info[parts[0][:-1]] = int(parts[1]) * 1024
        total = info.get("MemTotal")
        available = info.get("MemAvailable", total)
        if total is None:
            try:
                total = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
                available = total
            except (ValueError, OSError, AttributeError):
                total = available = 8 * 1024 ** 3
    limit = None
    text = _read("/sys/fs/cgroup/memory.max").strip()
    if text and text != "max":
        try:
            limit = int(text)
        except ValueError:
            limit = None
    return total, available, limit


def _gpus():
    gpus = []
    try:
        import pynvml

        pynvml.nvmlInit()
        for i in range(pynvml.nvmlDeviceGetCount()):
            h = pynvml.nvmlDeviceGetHandleByIndex(i)
            name = pynvml.nvmlDeviceGetName(h)
            name = name.decode() if isinstance(name, bytes) else str(name)
            gpus.append(GPUInfo(i, name, int(pynvml.nvmlDeviceGetMemoryInfo(h).total)))
        pynvml.nvmlShutdown()
        return gpus
    except Exception:
        pass
    if shutil.which("nvidia-smi"):
        try:
            out = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total",
                                  "--format=csv,noheader,nounits"], capture_output=True,
                                 text=True, timeout=2)
            for i, line in enumerate(out.stdout.strip().splitlines()):
                name, mem = [x.strip() for x in line.split(",")[:2]]
                gpus.append(GPUInfo(i, name, int(float(mem)) * 1024 ** 2))
        except Exception:
            pass
    return gpus


def _mpi():
    has = importlib.util.find_spec("mpi4py") is not None
    launcher = next((x for x in ("mpirun", "mpiexec", "srun") if shutil.which(x)), None)
    size = 1
    for var in ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "SLURM_NTASKS"):
        if os.environ.get(var):
            try:
                size = int(os.environ[var])
                break
            except ValueError:
                pass
    return MPIInfo(has, launcher, size)


def _backends():
    found = {name: importlib.util.find_spec(name) is not None
             for name in ("numba", "petsc4py", "mpi4py", "pymetis", "pyamg", "psutil",
                          "threadpoolctl", "meshio", "gmsh")}
    if found["petsc4py"]:
        from aa540fem.linalg.direct import mumps_available, petsc_available

        found["petsc4py"] = petsc_available()
        found["mumps"] = mumps_available()
        try:
            from petsc4py import PETSc

            found["cuda"] = bool(PETSc.Sys.hasExternalPackage("cuda"))
        except Exception:
            found["cuda"] = False
    return found


def _blas():
    try:
        from threadpoolctl import threadpool_info

        libs = {i["internal_api"] for i in threadpool_info()}
        return ", ".join(sorted(libs)) or "unknown"
    except Exception:
        return "unknown"


_probe_cache: HardwareInfo | None = None


def probe(refresh: bool = False) -> HardwareInfo:
    """Detect the hardware and the installed backends (cached; each detector
    is guarded, so this never raises)."""
    global _probe_cache
    if _probe_cache is not None and not refresh:
        return _probe_cache
    logical = os.cpu_count() or 1
    total, available, limit = _memory()
    _probe_cache = HardwareInfo(
        physical_cores=_physical_cores(logical), logical_cores=logical,
        usable_cores=_usable_cores(logical), ram_total=total, ram_available=available,
        memory_limit=limit, gpus=_gpus(), mpi=_mpi(), backends=_backends(), blas=_blas())
    return _probe_cache


# -- how much of it a run uses ----------------------------------------


@dataclass
class RunConfig:
    fraction: float                   # of the usable cores
    threads: int
    memory_budget: int                # bytes
    assembly_backend: str = "auto"    # numpy | numba
    direct_backend: str = "auto"      # superlu | petsc
    linear_backend: str = "scipy"     # scipy | petsc | petsc-cuda (Krylov path)
    ranks: int = 1                    # MPI ranks the run was launched with
    hardware: HardwareInfo | None = None

    def resolved_assembly(self) -> str:
        from aa540fem.backends import assembly_backends

        return assembly_backends()[0] if self.assembly_backend == "auto" else self.assembly_backend

    def resolved_direct(self) -> str:
        from aa540fem.linalg.direct import mumps_available

        if self.direct_backend == "auto":
            return "petsc" if mumps_available() else "superlu"
        return self.direct_backend

    def apply(self) -> RunConfig:
        """Set the thread counts of numba, the BLAS and OpenMP/MUMPS.

        The pools never run concurrently, so they all get ``threads``.
        """
        n = str(self.threads)
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "NUMBA_NUM_THREADS"):
            os.environ[var] = n
        try:
            import numba

            numba.set_num_threads(max(1, min(self.threads, numba.config.NUMBA_NUM_THREADS)))
        except (ImportError, ValueError):
            pass
        try:
            from threadpoolctl import threadpool_limits

            threadpool_limits(self.threads)
        except ImportError:
            pass
        return self

    def describe(self) -> str:
        gb = 1024 ** 3
        return (f"{self.fraction * 100:.0f} % of the machine: {self.threads} threads"
                + (f" x {self.ranks} ranks" if self.ranks > 1 else "")
                + f", {self.memory_budget / gb:.1f} GB budget, assembly {self.resolved_assembly()},"
                f" direct solver {self.resolved_direct()}, linear algebra {self.linear_backend}")

    def check_memory(self, matrix_or_nnz, n: int | None = None, force: bool = False):
        """Warn or refuse when a direct factorisation would exceed the budget."""
        from aa540fem.linalg.direct import estimate_direct_memory

        need = estimate_direct_memory(matrix_or_nnz, n)
        gb = 1024 ** 3
        if need > self.memory_budget and not force:
            raise MemoryError(f"the direct factorisation needs about {need / gb:.1f} GB, more "
                              f"than the {self.memory_budget / gb:.1f} GB budget: coarsen the "
                              f"mesh, raise the budget (--memory-budget / --machine 100), run "
                              f"under mpirun with the PETSc backend, or pass force=True")
        if need > 0.8 * self.memory_budget:
            warnings.warn(f"the direct factorisation needs about {need / gb:.1f} GB of the "
                          f"{self.memory_budget / gb:.1f} GB budget", stacklevel=2)
        return need


_config: RunConfig | None = None


def _fraction_from(text):
    if text is None:
        return None
    value = float(str(text).rstrip("%"))
    return value / 100.0 if value > 1.0 else value


def default_config(hw: HardwareInfo | None = None, fraction: float = 1.0,
                   threads: int | None = None) -> RunConfig:
    hw = hw or probe()
    fraction = min(1.0, max(0.0, fraction))
    ranks = max(1, hw.mpi.world_size)
    if threads is None:
        threads = max(1, int(round(fraction * hw.usable_cores / ranks)))
    available = hw.ram_available if hw.memory_limit is None else min(hw.ram_available,
                                                                     hw.memory_limit)
    budget = int(max(fraction, 0.1) * 0.8 * available)
    return RunConfig(fraction=fraction, threads=threads, memory_budget=budget, ranks=ranks,
                     hardware=hw)


def configure(fraction=None, threads=None, assembly=None, direct=None, linear=None,
              memory_budget=None, hw: HardwareInfo | None = None) -> RunConfig:
    """Build and apply the run configuration: arguments > environment > defaults."""
    global _config
    env = os.environ
    fraction = _fraction_from(fraction if fraction is not None
                              else env.get("AA540FEM_MACHINE_FRACTION"))
    if threads is None and env.get("AA540FEM_THREADS"):
        threads = int(env["AA540FEM_THREADS"])
    cfg = default_config(hw, 1.0 if fraction is None else fraction, threads)
    cfg.assembly_backend = assembly or env.get("AA540FEM_ASSEMBLY") or "auto"
    cfg.direct_backend = direct or env.get("AA540FEM_DIRECT") or "auto"
    cfg.linear_backend = linear or env.get("AA540FEM_LINEAR") or "scipy"
    budget = memory_budget if memory_budget is not None else env.get("AA540FEM_MEMORY_BUDGET")
    if budget is not None:
        cfg.memory_budget = int(float(budget) * 1024 ** 3)
    _config = cfg.apply()
    return _config


def get_config() -> RunConfig:
    """The current configuration (built from the environment and defaults on first use)."""
    if _config is None:
        configure()
    return _config


def recommend(hw: HardwareInfo | None = None, n_unknowns: int | None = None) -> str:
    """Advice on the best way to run on this machine (used by the CLI prompt)."""
    hw = hw or probe()
    lines = []
    if hw.backends.get("numba"):
        lines.append("threaded assembly (numba) is available")
    else:
        lines.append("install numba for threaded assembly: pip install aa540fem[numba]")
    if hw.backends.get("mumps"):
        lines.append("PETSc/MUMPS direct solver is available")
    elif hw.backends.get("petsc4py"):
        lines.append("PETSc is installed without MUMPS; SuperLU (one thread) does the "
                     "factorisations")
    else:
        lines.append("factorisations run on one thread (SuperLU); a PETSc build with MUMPS "
                     "uses all cores: pip install aa540fem[petsc]")
    if n_unknowns and n_unknowns > 2e5 and hw.mpi.mpi4py and hw.mpi.launcher:
        ranks = max(2, hw.usable_cores // 2)
        lines.append(f"{n_unknowns:.0f} unknowns: consider {hw.mpi.launcher} -n {ranks} "
                     f"python ... (distributed MUMPS)")
    if hw.gpus and hw.backends.get("cuda"):
        lines.append("PETSc has CUDA: the Krylov path can run on the GPU (--linear petsc-cuda)")
    elif hw.gpus:
        lines.append("an NVIDIA GPU is present but PETSc was built without CUDA; the GPU is "
                     "not used")
    return "; ".join(lines)
