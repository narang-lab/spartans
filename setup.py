# -*- coding: utf-8 -*-

import sys, os

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

setup(
    name="spartans",
    version="1.0.0",
    description="SpaRTaNS",
    author="Georgios Varnavides",
    license="MIT",
    packages=['spartans'], # find_packages(where = 'src'),
    package_dir={"": "src"},
    cmake_install_dir="src/spartans",
    cmake_args=['-DCMAKE_C_COMPILER=gcc','-DCMAKE_CXX_COMPILER=g++'],
    include_package_data = True,
)
