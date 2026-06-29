"""
Services package for CoreMet Application
"""

# Import services that don't depend on DGL first
from .data_service import DataService
from .molecular_service import MolecularService
from .docking_service import DockingService
from .cache_service import CacheService, OptimizedDataLoader
from .performance_service import PerformanceMonitor, MemoryOptimizer, CacheOptimizer, monitor_performance

# Import prediction service with DGL dependency handling
try:
    from .prediction_service import PredictionService
    PREDICTION_SERVICE_AVAILABLE = True
except ImportError as e:
    print(f"Warning: PredictionService not available due to DGL dependency: {e}")
    PredictionService = None
    PREDICTION_SERVICE_AVAILABLE = False

__all__ = [
    'DataService', 
    'MolecularService',
    'DockingService',
    'CacheService',
    'OptimizedDataLoader',
    'PerformanceMonitor',
    'MemoryOptimizer',
    'CacheOptimizer',
    'monitor_performance'
]

# Add PredictionService to __all__ only if available
if PREDICTION_SERVICE_AVAILABLE:
    __all__.append('PredictionService')
