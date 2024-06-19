from setuptools import setup, Extension

setup(
        name="matmult",
        version="0.0.1",
        description="matmult module",
        ext_modules=[Extension("matmult", sources=["matmultmodule.c"],
           include_dirs=["opt/OpenBLAS/include/lapacke.h"], 
           library_dirs=["/usr/lib/x86_64-linux-gnu","/opt/OpenBLAS/lib"],
           libraries=["lapack", "openblas"])]
)
