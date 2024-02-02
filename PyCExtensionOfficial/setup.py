from setuptools import setup, Extension

setup(
   name="spam",
   version="0.0.1",
   description="spam module",
   ext_modules= [Extension("spam", sources=["spam.c"])]
   
)
