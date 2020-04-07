from setuptools import setup, find_packages

base = None

setup(
    name="hiMoon",
    version="0.02",
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Console",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry"
    ],
    description="Match NGS data with defined haplotypes",
    author="Solomon M. Adams, PharmD, PhD",
    author_email="sadams07@su.edu",
    packages=find_packages(),
    package_data={"hiMoon.tests": ["*"]},
    entry_points={"console_scripts": ["hiMoon = hiMoon.__main__:main"]},
    install_requires=[
        "pulp",
        "pandas",
        "pysam"
    ]
)

