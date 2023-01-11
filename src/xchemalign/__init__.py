from importlib.metadata import version

__version__ = version("xchemalign")
del version

__all__ = ["__version__"]
