from setuptools import setup, Extension

setup(
        name="basicfopenquery",
        version="0.0.1",
        description="basicfopenquery",
        ext_modules=[Extension("basicfopenquery", sources=["basicfopenquerymodule.c"])]
)
