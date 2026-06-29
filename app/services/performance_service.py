"""
Performance monitoring service for CoreMet Application
"""

import time
import psutil
import gc
from typing import Dict, Any, Optional
from functools import wraps
import logging
from datetime import datetime

from app.config import Config


class PerformanceMonitor:
    """Service for monitoring application performance"""
    
    def __init__(self):
        self.config = Config()
        self._metrics: Dict[str, Any] = {}
        self._start_time = time.time()
        
        # Setup logging
        self._setup_logging()
    
    def _setup_logging(self):
        """Setup performance logging"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler('performance.log'),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger('performance')
    
    def get_system_metrics(self) -> Dict[str, Any]:
        """Get current system performance metrics"""
        try:
            # CPU metrics
            cpu_percent = psutil.cpu_percent(interval=1)
            cpu_count = psutil.cpu_count()
            
            # Memory metrics
            memory = psutil.virtual_memory()
            memory_percent = memory.percent
            memory_available = memory.available / (1024**3)  # GB
            memory_used = memory.used / (1024**3)  # GB
            
            # Disk metrics
            disk = psutil.disk_usage('/')
            disk_percent = disk.percent
            disk_free = disk.free / (1024**3)  # GB
            
            return {
                'timestamp': datetime.now().isoformat(),
                'cpu': {
                    'percent': cpu_percent,
                    'count': cpu_count
                },
                'memory': {
                    'percent': memory_percent,
                    'available_gb': round(memory_available, 2),
                    'used_gb': round(memory_used, 2)
                },
                'disk': {
                    'percent': disk_percent,
                    'free_gb': round(disk_free, 2)
                },
                'uptime_seconds': round(time.time() - self._start_time, 2)
            }
        except Exception as e:
            self.logger.error(f"Error getting system metrics: {e}")
            return {}
    
    def get_python_metrics(self) -> Dict[str, Any]:
        """Get Python-specific performance metrics"""
        try:
            import sys
            
            # Python memory usage
            memory_info = psutil.Process().memory_info()
            memory_mb = memory_info.rss / (1024**2)
            
            # Garbage collection stats
            gc_stats = gc.get_stats()
            
            return {
                'python_version': sys.version,
                'memory_usage_mb': round(memory_mb, 2),
                'gc_stats': gc_stats,
                'gc_count': gc.get_count()
            }
        except Exception as e:
            self.logger.error(f"Error getting Python metrics: {e}")
            return {}
    
    def log_performance(self, operation: str, duration: float, **kwargs):
        """Log performance metrics for an operation"""
        metrics = {
            'operation': operation,
            'duration_seconds': round(duration, 3),
            'timestamp': datetime.now().isoformat(),
            **kwargs
        }
        
        self.logger.info(f"Performance: {operation} took {duration:.3f}s")
        self._metrics[operation] = metrics
    
    def get_performance_summary(self) -> Dict[str, Any]:
        """Get performance summary"""
        return {
            'system_metrics': self.get_system_metrics(),
            'python_metrics': self.get_python_metrics(),
            'operation_metrics': self._metrics
        }


def monitor_performance(operation_name: Optional[str] = None):
    """Decorator to monitor function performance"""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            start_time = time.time()
            operation = operation_name or f"{func.__module__}.{func.__name__}"
            
            try:
                result = func(*args, **kwargs)
                duration = time.time() - start_time
                
                # Log performance
                monitor = PerformanceMonitor()
                monitor.log_performance(operation, duration, success=True)
                
                return result
            except Exception as e:
                duration = time.time() - start_time
                
                # Log error performance
                monitor = PerformanceMonitor()
                monitor.log_performance(operation, duration, success=False, error=str(e))
                
                raise
        
        return wrapper
    return decorator


class MemoryOptimizer:
    """Service for memory optimization"""
    
    def __init__(self):
        self.monitor = PerformanceMonitor()
    
    def optimize_memory(self):
        """Perform memory optimization"""
        try:
            # Force garbage collection
            collected = gc.collect()
            
            # Log memory optimization
            self.monitor.log_performance(
                "memory_optimization", 
                0, 
                objects_collected=collected
            )
            
            return collected
        except Exception as e:
            self.monitor.logger.error(f"Error optimizing memory: {e}")
            return 0
    
    def get_memory_usage(self) -> Dict[str, Any]:
        """Get detailed memory usage information"""
        try:
            import tracemalloc
            
            # Start tracing if not already started
            if not tracemalloc.is_tracing():
                tracemalloc.start()
            
            # Get current memory usage
            current, peak = tracemalloc.get_traced_memory()
            
            return {
                'current_mb': round(current / (1024**2), 2),
                'peak_mb': round(peak / (1024**2), 2),
                'gc_count': gc.get_count(),
                'gc_stats': gc.get_stats()
            }
        except Exception as e:
            self.monitor.logger.error(f"Error getting memory usage: {e}")
            return {}
    
    def check_memory_threshold(self, threshold_mb: float = 1000) -> bool:
        """Check if memory usage exceeds threshold"""
        try:
            memory_info = psutil.Process().memory_info()
            memory_mb = memory_info.rss / (1024**2)
            
            if memory_mb > threshold_mb:
                self.monitor.logger.warning(
                    f"Memory usage ({memory_mb:.2f} MB) exceeds threshold ({threshold_mb} MB)"
                )
                return True
            
            return False
        except Exception as e:
            self.monitor.logger.error(f"Error checking memory threshold: {e}")
            return False


class CacheOptimizer:
    """Service for cache optimization"""
    
    def __init__(self):
        self.monitor = PerformanceMonitor()
    
    def optimize_cache(self, cache_service):
        """Optimize cache performance"""
        try:
            start_time = time.time()
            
            # Get cache statistics
            stats = cache_service.get_cache_stats()
            
            # Clear expired entries
            cache_service.clear_cache()
            
            duration = time.time() - start_time
            self.monitor.log_performance(
                "cache_optimization",
                duration,
                cache_size_before=stats.get('cache_size', 0),
                memory_usage_before=stats.get('memory_usage', 0)
            )
            
        except Exception as e:
            self.monitor.logger.error(f"Error optimizing cache: {e}")
    
    def get_cache_efficiency(self, cache_service) -> Dict[str, Any]:
        """Get cache efficiency metrics"""
        try:
            stats = cache_service.get_cache_stats()
            
            return {
                'cache_size': stats.get('cache_size', 0),
                'memory_usage': stats.get('memory_usage', 0),
                'cached_keys': stats.get('cached_keys', []),
                'efficiency_score': self._calculate_efficiency_score(stats)
            }
        except Exception as e:
            self.monitor.logger.error(f"Error getting cache efficiency: {e}")
            return {}
    
    def _calculate_efficiency_score(self, stats: Dict[str, Any]) -> float:
        """Calculate cache efficiency score"""
        try:
            cache_size = stats.get('cache_size', 0)
            memory_usage = stats.get('memory_usage', 0)
            
            if cache_size == 0:
                return 0.0
            
            # Simple efficiency score based on cache size vs memory usage
            efficiency = min(1.0, cache_size / max(1, memory_usage / 1000))
            return round(efficiency, 3)
        except Exception:
            self.monitor.logger.debug("Could not calculate efficiency score", exc_info=True)
            return 0.0
