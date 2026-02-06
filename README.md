# gris-iceocean-process
Processing workflow for the Greenland ocean forcing datasets for ISMIP7, developed by the Greenland Ocean Forcing Focus Group.

# Overview
The notebooks and Python modules in this repository are used to create ocean forcing datasets for standalone ice sheet model experiments for the Ice Sheet Model Intercomparison for CMIP7 (ISMIP7). Two datasets are produced that are used to calculate ocean forcing within an ice sheet model:
1. Ocean thermal forcing
2. Runoff

Together, these two datasets can be used to calculate submarine melt rates at outlet glacier termini around the ice sheet. The ice sheet models then simulate iceberg calving processes and melt+calving together can be used to calculate the retreat or advance of outlet glacier termini through time.

# Organization of the repository
In order to run the code within this repository, several required Python modules must be installed. All dependencies are listed in `requirements.txt`. One option is to setup a virtual environment and use this file to install all dependencies like this:
```
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Notebooks and code are organized into two directories:
1. [oceanTF/](oceanTF): this directory contains code to process ocean thermal forcing datasets
2. [runoff/](runoff): this directory contains code to process runoff datasets

The oceanTF code must be run first because the code within the runoff workflow relies on the output of the oceanTF code.
