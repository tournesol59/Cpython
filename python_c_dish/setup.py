from setuptools import setup, Extension

setup(
        name="dish",
        version="0.0.1",
        description="dish module",
        ext_modules=[Extension("dish", sources=["dishmodule.c"])]
)
