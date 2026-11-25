import importlib
import platform
import sys

print(platform.platform(), sys.executable)
for name in ("numpy", "scipy", "mpmath", "sympy", "matplotlib"):
    try:
        mod = importlib.import_module(name)
        print(name, mod.__version__)
    except ImportError as exc:
        print(exc)
