#!/usr/bin/env python
from setuptools import setup, find_packages

version = '0.1'
author = 'Pascal Wolf'
author_email = 'pascal.wolf@posteo.de'

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='pyradiate',
    version=version,
    description='Gamma spectroscopy and (equivalent) dose calculations for radioctive sources (and samples) in Python',
    url='https://github.com/leloup314/pyradiate',
    license='MIT License',
    long_description='',
    author=author,
    maintainer=author,
    author_email=author_email,
    maintainer_email=author_email,
    packages=find_packages(),
    setup_requires=['setuptools'],
    install_requires=required,
    include_package_data=True,  # accept all data files and directories matched by MANIFEST.in or found in source control
    package_data={'': ['README.*', 'VERSION'], 'docs': ['*'], 'examples': ['*']},
    keywords=['spectroscopy', 'gamma', 'dose', 'sievert', 'radioactive', 'source', 'radiation', 'safety'],
    platforms='any'
)
