"""Drive DIgSILENT PowerFactory for Porthos: model dump and RMS fault runs.

This package only drives PowerFactory through its Python API and writes what PowerFactory
computes (its CSV export) plus JSON records. It does no numerics: every comparison and
analysis is done in Julia (`src/io/powerfactory.jl`). It uses the standard library only and
runs with the Python version PowerFactory ships its module for (`pf/config.json`).
"""

__all__ = ["session", "inspect_model", "simulate"]
