#!/usr/bin/env python3
from sys import version_info, stderr
from setuptools import setup, find_packages

REQUIRED_PYTHON: tuple = (3, 10)
CURRENT_PYTHON = version_info[:2]

if CURRENT_PYTHON < REQUIRED_PYTHON:
    stderr.write(
        f"invpg requires Python {REQUIRED_PYTHON[0]}.{REQUIRED_PYTHON[1]} or higher and your current version is {CURRENT_PYTHON}.")
    exit(1)

setup(
    name="invpg",
    version='0.0.1',
    description='A tool to annotate inversions from pangenome graph bubbles.',
    url='https://github.com/SandraLouise/INVPG_annot',
    author='Sandra Romain',
    author_email='sandra.romain@inria.fr',
    packages=['invpg'],
    include_package_data=True,
    zip_safe=False,
    license="LICENSE.md",
    long_description=open("README.md", encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    install_requires=open("requirements.txt", encoding='utf-8').readlines(),
    entry_points={'console_scripts': ['invpg=invpg.main:main']}
)
