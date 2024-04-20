from setuptools import setup, Extension

setup(
        name="lhlib",
        version="0.0.1",
        description="lhlib module",
        ext_modules=[Extension("lhlib", sources=["lhlibmodule.c"])]
)
