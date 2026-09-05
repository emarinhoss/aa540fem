"""ParaView time-series (``.pvd`` collection) output."""

from __future__ import annotations

import pathlib


def write_series(prefix, times, write_step):
    """Write one file per stored time plus a ``.pvd`` collection.

    ``write_step(path, index)`` writes the ``index``-th stored state to
    ``path`` (a ``.vtu``); the files are ``<prefix>_0000.vtu`` ... and the
    collection ``<prefix>.pvd`` (open the latter in ParaView to animate).
    Returns the path of the ``.pvd`` file.
    """
    prefix = pathlib.Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    entries = []
    for i, t in enumerate(times):
        name = f"{prefix.name}_{i:04d}.vtu"
        write_step(prefix.parent / name, i)
        entries.append(f'    <DataSet timestep="{float(t)!r}" group="" part="0" file="{name}"/>')
    pvd = prefix.with_suffix(".pvd")
    pvd.write_text(
        '<?xml version="1.0"?>\n'
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n'
        "  <Collection>\n" + "\n".join(entries) + "\n  </Collection>\n</VTKFile>\n")
    return pvd
