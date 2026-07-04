from setuptools import setup

setup(
   name='basicsumstd',
   version='0.0.1',
   install_requires=[
      'requests',
      'importlib-metadata; python_version<"3.10"',
   ],

   packages=find_packages(
       where='src',
       include=['basicsumstd*'],
    ),
)
