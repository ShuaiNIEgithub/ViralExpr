#!/usr/bin/env python
from setuptools import setup, find_packages
from viralexpr.__version__ import version

with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='viralexpr',
    version=version,
    description='ViralExpress: universal viral gene expression quantification from host RNA-seq data',
    url='https://github.com/your-repo/viralexpr',
    author='Your Name',
    license='MIT',
    python_requires='>=3.8',
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'viralexpr = viralexpr.pipe:main',
        ],
    },
    install_requires=[
        'pandas>=1.0',
        'numpy>=1.18',
        'matplotlib>=3.2',
        'seaborn>=0.11',
        'scipy>=1.5',
    ],
)