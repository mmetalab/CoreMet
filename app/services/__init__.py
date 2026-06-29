"""
Services package for CoreMet.

Imports are lazy (PEP 562): the heavy optional dependencies — torch/dgl (via
cache_service and prediction_service) and rdkit (via molecular_service) — are NOT
loaded at package import time. They load only when the corresponding class is first
accessed (i.e., when a prediction is actually run). This keeps the database web app's
startup footprint small (~0.3 GB) so it stays well within the 2 GB instance limit;
loading them eagerly previously added ~0.6 GB and contributed to OOM-killed workers.
"""
import importlib

# public name -> submodule providing it
_LAZY = {
    "DataService": "data_service",
    "MolecularService": "molecular_service",
    "DockingService": "docking_service",
    "CacheService": "cache_service",
    "OptimizedDataLoader": "cache_service",
    "PerformanceMonitor": "performance_service",
    "MemoryOptimizer": "performance_service",
    "CacheOptimizer": "performance_service",
    "monitor_performance": "performance_service",
    "PredictionService": "prediction_service",
}

__all__ = list(_LAZY)


def __getattr__(name):
    target = _LAZY.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = importlib.import_module(f".{target}", __name__)
    return getattr(module, name)
