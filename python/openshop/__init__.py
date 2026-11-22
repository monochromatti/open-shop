"""Nonlinear hydropower scheduling. Physical units are m³, m³/s, m, h and MW."""

from ._openshop import __version__, solve

__all__ = ["solve", "__version__"]
