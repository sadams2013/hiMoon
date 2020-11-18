from setuptools import setup, find_packages

base = None

setup(
    name="hiMoon",
    version="0.12.3",
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: POSIX",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry"
    ],
    description="Define names haplotypes from data stored in VCF",
    author="Solomon M. Adams, PharmD, PhD",
    author_email="sadams07@su.edu",
    packages=find_packages(),
    package_data={"hiMoon.tests": ["*"], "": ["template*.vcf"]},
    entry_points={"console_scripts": ["hiMoon = hiMoon.__main__:main"]},
    install_requires=[
        "pulp",
        "pandas",
        "pysam",
        "numpy"
    ]
)

