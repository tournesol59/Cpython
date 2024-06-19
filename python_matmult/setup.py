from setuptools import setup, Extension

setup(
        name="matclass",
        version="0.0.1",
        description="matclass module",
        ext_modules=[Extension("matclass", sources=["matclassmodule.c"],
           include_dirs=["opt/OpenBLAS/include/lapacke.h"], 
           library_dirs=["/usr/lib/x86_64-linux-gnu","/opt/OpenBLAS/lib"],
           libraries=["lapack", "openblas"])]
)
