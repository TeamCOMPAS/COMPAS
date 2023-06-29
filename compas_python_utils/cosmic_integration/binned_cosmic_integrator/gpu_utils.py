try:
    import cupy as cp
    gpu_available = True

except (ModuleNotFoundError, ImportError):
    import numpy as np
    gpu_available = False

if gpu_available:
    xp = cp

else:
    xp = np