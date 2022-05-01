---
sidebar_position: 2
---

# Installation Instructions

SpaRTaNS is a Python 3 package with a custom c++ binding (compiled using [pybind11](https://pybind11.readthedocs.io/en/stable/index.html)), which relies heavily on MPI and HDF5, through numpy, mpi4py, and h5py.

## Building from source

Currently, SpaRTaNS is distributed from source and built using the [scikit-build](https://scikit-build.readthedocs.io/en/latest/#) package.

To build and install the most recent version of SpaRTaNS, you can clone the main branch from the source repository and install:

``` bash
git clone https://github.com/narang-lab/spartans
cd spartans
pip install .
```

### Troubleshooting

The scikit-build should handle all dependencies automatically.
If you find this is not the case, you can try installing all the dependencies manually before installing SparTaNS:

- "setuptools"
- "wheel"
- "pybind11"
- "scikit-build"
- "cmake"
- "ninja"
- "mpi4py==3.0.3"
- "numpy"
- "h5py"
- "PyYAML"
- "parse"

:::note

Notice the version-pinning on mpi4py, which is likely the installation issues culprit.

:::

