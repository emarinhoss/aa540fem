"""Command-line helpers shared by the examples: run-configuration flags and
the "run at 50 % or 100 % of the machine" prompt.

    parser = argparse.ArgumentParser()
    add_run_arguments(parser)
    args = parser.parse_args()
    config = configure_from_args(args)      # prompts once on a terminal

Precedence: flags > environment (``AA540FEM_*``) > prompt > defaults (all
usable cores, 80 % of the free memory).  The prompt is skipped with
``--no-prompt``, ``AA540FEM_NO_PROMPT=1``, when a flag already fixes the
fraction, when stdin is not a terminal (batch jobs, CI, tests) or under
``mpirun``.
"""

from __future__ import annotations

import argparse
import os
import sys

from aa540fem.hardware import HardwareInfo, RunConfig, configure, default_config, probe, recommend

FRACTIONS = {"1": 0.5, "2": 1.0}


def add_run_arguments(parser: argparse.ArgumentParser) -> None:
    g = parser.add_argument_group("run configuration")
    g.add_argument("--machine", choices=["50", "100"], default=None,
                   help="use 50 %% or 100 %% of the cores and memory (asked on a terminal)")
    g.add_argument("--threads", type=int, default=None, help="threads for assembly and solvers")
    g.add_argument("--assembly", choices=["auto", "numpy", "numba"], default=None)
    g.add_argument("--direct", choices=["auto", "superlu", "petsc"], default=None,
                   help="sparse direct solver")
    g.add_argument("--linear", choices=["scipy", "petsc", "petsc-cuda"], default=None,
                   help="linear-algebra backend of the Krylov path")
    g.add_argument("--memory-budget", type=float, default=None, help="GB")
    g.add_argument("--no-prompt", action="store_true", help="never ask; use the defaults")


def prompt_fraction(hw: HardwareInfo, stream=None, input_fn=input) -> float:
    """Ask for 50 % or 100 % of the machine; Enter takes 100 %."""
    stream = stream or sys.stdout
    half = default_config(hw, 0.5)
    full = default_config(hw, 1.0)
    gb = 1024 ** 3
    stream.write(f"aa540fem: {hw.describe()}\n          {recommend(hw)}\n")
    stream.write(f"Run at [1] 50 % ({half.threads} threads, {half.memory_budget / gb:.1f} GB "
                 f"budget) or [2] 100 % ({full.threads} threads, "
                 f"{full.memory_budget / gb:.1f} GB)?  [2]: ")
    stream.flush()
    try:
        answer = input_fn().strip()
    except EOFError:
        answer = ""
    return FRACTIONS.get(answer or "2", 1.0)


def configure_from_args(args, interactive: bool | None = None, stream=None,
                        input_fn=input) -> RunConfig:
    """Run configuration from the parsed flags (see :func:`add_run_arguments`)."""
    hw = probe()
    fraction = None if args.machine is None else float(args.machine) / 100.0
    env_fraction = os.environ.get("AA540FEM_MACHINE_FRACTION")
    if interactive is None:
        interactive = (not args.no_prompt and fraction is None and env_fraction is None
                       and not os.environ.get("AA540FEM_NO_PROMPT")
                       and hw.mpi.world_size == 1 and sys.stdin is not None
                       and sys.stdin.isatty())
    if interactive and fraction is None:
        fraction = prompt_fraction(hw, stream, input_fn)
    return configure(fraction=fraction, threads=args.threads, assembly=args.assembly,
                     direct=args.direct, linear=args.linear, memory_budget=args.memory_budget,
                     hw=hw)


def main(argv=None):
    """``python -m aa540fem.cli``: print what the machine offers and the defaults."""
    parser = argparse.ArgumentParser(description="aa540fem hardware probe")
    add_run_arguments(parser)
    args = parser.parse_args(argv)
    cfg = configure_from_args(args, interactive=False)
    print(probe().describe())
    print(recommend(probe()))
    print(cfg.describe())


if __name__ == "__main__":
    main()
