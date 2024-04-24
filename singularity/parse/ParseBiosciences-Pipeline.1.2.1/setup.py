#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

from setuptools import setup

# From https://stackoverflow.com/questions/40031422/python-recursively-include-package-data-in-setup-py
from pathlib import Path

datadir = Path(__file__).parent / "splitpipe"
files = [str(p.relative_to(datadir)) for p in datadir.rglob("*")]

setup(
    name="splitpipe",
    version="1.2.1",
    # Generic url, description
    url="https://www.parsebiosciences.com/",
    description="Parse Bioscience's data analysis pipeline",
    # Package and top-level script
    packages=["splitpipe"],
    scripts=["split-pipe"],
    # Requirements
    python_requires=">= 3.8",
    install_requires=[
        "numpy >= 1.23",
        "pandas",
        "scipy",
        "matplotlib",
        "natsort",
        "h5py",
        "pysam >= 0.20.0",
        "scanpy >= 1.8",
        "anndata",
        "louvain",
        "leidenalg",
        "python-igraph",
        "jinja2",
        "psutil",
        # "openpyxl",
        # "openpyxl >= 3.1.0",  # for SampleLoadingTable; Maybe issues with pandas
        "python-calamine",      # Replacement for openpyxl
    ],
    # Data for package
    package_data={"splitpipe": files},
    # package_data={
    #    'splitpipe': ['config/*', 'barcodes/*', 'templates/*', 'scripts/*', 'TCR/*'],
    # },
    include_package_data=True,
)
