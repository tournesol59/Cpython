from setuptools import setup, Extension

setup(
        name="matmult",
        version="0.0.1",
        description="matmult module",
        ext_modules=[Extension("matmult", sources=["matmultmodule.c"])]
)
