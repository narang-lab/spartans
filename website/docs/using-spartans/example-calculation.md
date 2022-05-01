---
sidebar_position: 1
---

# Example Calculation

## Directory Structure

Navigate to the `SO3-cube-periodic-diffusive-absorbing` directory under `examples` in the source repo and list the directory's contents:

``` bash 
cd examples/SO3-cube-periodic-diffusive-absorbing
tree
```

You should see the following structure:

``` bash
.
├── clean.sh
├── in-files
│   └── spartans_test_SO3_cube-periodic_dataset-compressed.h5
├── log-files
├── notebooks
│   └── spartans-examples_SO3-cube-periodic-diffusive-absorbing.nb
├── out-files
├── parser.py
└── visualizations
    └── steady-state-horizontal-flux.png

5 directories, 5 files
```

:::info

The `notebooks` and `visualizations` directory are optional, but we find it's a natural way to organize the files used to generate the inputs and visualize the outputs. You can find an example walkthrough of how we generate various input files under [Tutorials](../tutorials/double-chamber-flow/intro.mdx).

:::

The default directory structure consists of the following three directories:
  - `in-files`: Directory SpaRTaNS looks in for the `*hdf5`/`*.h5` user-input files
  - `out-files`: Directory SpaRTaNS will export output files
  - `log-files`: Directory SpaRTaNS will write log files

You can customize these directories in `parser.py`, which we elaborate on below.

Finally, there's a utility `clean.sh` script, which you might find useful when running calculations:

``` bash
#!/bin/bash
rm out-files/*
rm log-files/*
```

## Understanding `parser.py`

SpaRTaNS is meant to be executed from the command line as follows:

``` bash
mpirun [-n num_processes] python parser.py
```

Running the line above should result in output similar to
``` yaml
2022-05-01 11:41:44,786 : 
______             ____ _____      _   _______
\  ___)           |  _ (_   _)    | \ | \  ___)
 \ \  ______ __  _| |_) )| | __  _|  \| |\ \
  > >(  __  )  \/ /  __/ | |/  \/ /     | > >
 / /__| || ( ()  <| |    | ( ()  <| |\  |/ /__
/_____)_||_|\__/\_\_|    |_|\__/\_\_| \_/_____)


2022-05-01 11:41:44,786 : Setting up Database dictionary.
2022-05-01 11:41:44,786 : Building Database object for files in ./in-files
2022-05-01 11:41:44,786 : Preparing Database dry run.
2022-05-01 11:41:44,786 : Parsing file ./in-files/spartans_test_SO3_cube-periodic_dataset-compressed.h5
2022-05-01 11:41:44,792 : SpaRTaNS will import the following datasets:
--------------------------------------------------------------------------------
connectivity_info:
  bounce_tensors:
    bounce_00: (48, 48)
    bounce_01: (48, 48)
    bounce_02: (48, 48)
  connectivity: (1, 1)
  connectivity(structure(mesh(000)--material(A))--structure(mesh(000)--material(A))):
    bounce_00: (1452, 4)
    bounce_01: (1452, 4)
    bounce_02: (1452, 4)
  connectivity.ordering:
  - '000'
material(A):
  diagonal: (48)
  frequencies: (48)
  mixing_matrix: (48, 48)
  velocities: (48, 3)
mesh(000):
  surface_normals: (1452, 3)
  tetrahedra_indices: (7986, 4)
  triangle_indices: (1452, 3)
  vertices: (1728, 3)
structure(mesh(000)--material(A)):
  body_injection: (48, 1728)
--------------------------------------------------------------------------------

ABORT now if this is not what you expected.
```

The output above is SpaRTaNS inferring the user-input dataset shapes by looking into the `database_directory`.
This, along with the `output_directory` and `log_directory` can be changed by chaning the following lines in `parser.py`:

``` python
    from spartans.logger_builder import setup_logger
    from spartans.parser_utils import print_spartans_logo

# highlight-start
    config_dict = {
            'directories': {
                'database_directory': './in-files',
                'output_directory': './out-files',
                'log_directory': './log-files'
                }
            }
# highlight-end

    # Set up log files
    log_directory = config_dict['directories']['log_directory']
    setup_logger('main_log', '{}/main_log.txt'.format(log_directory))
    main_log = logging.getLogger('main_log')

