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
```

SpaRTaNS then proceeds to set up the necessary `Structure` and `Interface` objects (see [API Design section](./api-design.md)).  

SpaRTaNS solves the BTE iteratively (see [Formalism section](../formalism/boltzmann-transport-theory.md)).
As such, inputs/outputs are stored internally for each structure and iteration number.

However, to avoid large output size, and since the steady-state distribution is given by the accumulated distribution functions at each iteration, SpaRTaNS only writes the accumulated distribution and surface fluxes in `{outDir}/accumulated_{structureId}.h5` by default.
If individual distribution outputs are required (e.g. for debugging), the output files 'verbosity' can be controlled by changing the following lines:

```python
    # Set up Structure objects
    main_log.info('Setting up Structure objects.')
    structure_keys = [struct_key for struct_key in database_dict.keys() if 'structure(' in struct_key]
# highlight-start
    structures = [Structure(comm,config_dict,database_dict,struct_key, write_debug_outputs=False, 
                  debug_frequency=1) for struct_key in structure_keys]
# highlight-end
```

Finally, SpaRTaNS uses a `Scheduler` object (see [API Design section](./api-design.md)) to handle the various collision, and bounce operators in a logical sequence.
The default `Scheduler` uses the following pseudo-code:

```python
for struct in structures:
    # First push from structures interfaces
    for interface_id, interface_object in struct.interfaces.items():
        interface_object.push_to(struct)

    # Then, run structure
    struct.run()

    # Finally, pull to structures interfaces
    for interface_id, interface_object in struct.interfaces.items():
        interface_object.pull_from()
```

This can be changed by modifying the following lines in `parser.py`:
```python
    # Set up Scheduler
    # Note: change this to fit your needs, and/or define new schedulers in scheduler.py
    main_log.info('Setting up Scheduler.')

# highlight-start
    from spartans.scheduler import Unweighted_simple as Scheduler
    scheduler = Scheduler(
            min_condition=50,
            structures=structures
            )
#highlight-end
```
## C++ Streaming Operator Output

In addition to the collision and bounce operator logs from the python side, SpaRTaNS also outputs the streaming operator logs from the c++ side, with output similar to:

``` bash
# highlight-next-line
nTets: 7986  nUniqueTets: 6  expectedSpeedup: 1331x
nFaces: 16698  nExtFaces: 1452
Explicit event mode.
nBins: 48  nEvents: 48  nVerts: 1728
nEvents by process: [ 12 12 12 12 ]

Processing events: 8% 17% 25% 33% 42% 50% 58% 67% 75% 83% 92% 100% done.
Applying Oinv: done.

# highlight-start
Sum rule check:
Largest sum rule error: 5.551115e-15, at ie: 8
# highlight-end
```

Two points are of interest in the output above:  
- SpaRTaNS recognizes the mesh is very regular and only performs streaming in the 6 `nUniqueTets` tetrahedra, rotating and translating where appropriate to achieve an expected speedup of 1331x. As such, when the geometry allows it, it's beneficial to use a regular mesh.
- SpaRTaNS reports the largest sum-rule error ((Integrated `Body_in` + Integrated `Surface_in`) - (Integrated `Body_out` + Integrated `Surface_out`)) across states. You'll want to check this is reasonably small across all runs.
