Table of Contents
=================

1.  [Libraries](#libraries)
2.  [Main Function](#main-function)
    1.  [Syntax Check](#syntax-check)
    2.  [Seeds](#seeds)
    3.  [Temperature](#temperature)
    4.  [File Input](#file-input)
    5.  [Other Operations](#other-operations)
3.  [Other Functions](#other-functions)
4.  [run\_simulation.sh Script](#run_simulation.sh-script)
    1.  [Introduction](#introduction)
    2.  [SETTINGS](#settings)
        1.  [GPU architecture](#gpu-architecture)
        2.  [PATH FOR GSL LIBRARY](#path-for-gsl-library)
        3.  [Temperature range](#temperature-range)
        4.  [Number of replicas](#number-of-replicas)
        5.  [Number of Monte Carlo iterations](#number-of-monte-carlo-iterations)
    3.  [RUN](#run)
    4.  [OUTPUT](#output)
    5.  [POST-PROCESSING](#post-processing)
    6.  [EXAMPLES](#examples)
    
    
SMrandomTetrads
===============

SMrandomTetrads is a program written by many people, and has a history of almost ten years. The goal of this program is to simulate the 4-phasor random laser model, which is known to be a disordered spin glass system. Here we explain how to use this tool.

Libraries
---------

The following libraries are included in the code:

*   `iostream`
*   `fstream`
*   `vector`
*   `queue`
*   `algorithm`
*   `cmath`
*   `iomanip`
*   `sys/time.h`
*   `stdlib.h`
*   `stdio.h`
*   `cuda.h`
*   `curand.h`
*   `curand_kernel.h`
*   `thrust/host_vector.h`
*   `thrust/device_vector.h`
*   `thrust/device_ptr.h`
*   `thrust/generate.h`
*   `thrust/reduce.h`
*   `thrust/functional.h`
*   `thrust/fill.h`
*   `thrust/copy.h`
*   `thrust/execution_policy.h`
*   `SMrandomTetradsRUNCHECK.h`
*   `SMrandomTetradsSettings.h`
*   `gtcStructures.h`
*   `SMrandomTetrads_structures.h`
*   `hashing.cpp`
*   `graph.cpp`
*   `generate4plets.cpp`
*   `gsl_randist.h`
*   `gsl_rng.h`
*   `tetrads.cpp`
*   `generateQuadsFully.cpp`
*   `parallelMCstep.cu`
*   `functionsSM_CPU.cpp`
*   `SMrandomTetrads_CPU_GPU_initializations_v2.h`
*   `gtcTools.h`
*   `time.h`
*   `resumeToolsGTC.h`

Main Function
-------------

The main function of the code takes in command line arguments, including the number of arguments, and performs various operations based on these arguments.

### Syntax Check

The function `checkSyntax()` is called to check that the correct number of arguments have been passed in.

### Seeds

The code generates and sets random seeds for the program, including a master seed and seeds for the graph and replicas.

### Temperature

The code also handles the temperature range and intervals, and stores the values in an array.

### File Input

A file named "input\_analysis.dat" is opened and data is written to it.

### Other Operations

Other operations are performed in the code, such as checking for certain conditions and exiting the program if they are not met.

Other Functions
---------------

In addition to the main function, the code includes various other functions for performing specific tasks, such as generating and handling seeds, checking syntax, and writing to files. These functions are called from within the main function to complete the program's overall functionality.

run\_simulation.sh Script
=========================

Introduction
------------

The `run_simulation.sh` script is a Bash script that is used to configure and run a Monte Carlo simulation of a physical system. It sets various parameters for the simulation, such as the GPU architecture, the temperature range, the number of replicas, and the number of Monte Carlo iterations.

SETTINGS
--------

### GPU architecture

The script starts by setting the GPU architecture that will be used for the simulation. This is done by setting the value of the `arch` variable to the appropriate value. The values that are currently supported are 30 (VISNU/DURGA), 35 (KRAKEN), and 70 (ECATON). It is important to note that before running the script, the value of this variable should be checked and set to the correct value according to the README.md file, SEC-I.

### PATH FOR GSL LIBRARY

The script also sets the path for the GNU Scientific Library (GSL), which is a numerical library that provides a wide range of mathematical functions. The path is set by the `gsl_path` variable. For example, on the developer's laptop, the path is set to `/usr/local/include/gsl`, while on a cluster, it is set to `/usr/include/gsl`.

### Simulation settings

The script then sets various simulation settings, such as the size of the system, the temperature range, the number of PT replicas, and the number of real replicas. These settings include:

*   `Size`: Number of modes of the system.
*   `t_min` and `t_max`: Temperature range of the system.
*   `number_of_PT_replicas`: Number of temperatures between `t_min` and `t_max`. They will be linearly spaced.
*   `number_of_real_replicas`: Number of independent replicas of the system.
*   `PT_flag`: Flag for the parallel exchange algorithm. If set to 1, only the equilibrium properties will be looked at, otherwise it should be set to 0.
*   `number_of_PT_step`: How often the exchange between two nearby replicas in energy is proposed, in Monte Carlo iterations.
*   `power_of_iterations`: The number of iterations that we want to run the Monte Carlo simulation for. The real number of iterations will be calculated as `(2^power_of_iterations)/number_of_PT_step`.

### Tuning settings

The script also includes settings that should not be modified as doing so is not recommended. These settings include `print_config`, which controls how often the configuration is printed, and `frequency_mode`, which indicates how the frequencies are generated.

Script Execution
----------------

The script then calculates the number of iterations to run the simulation for, and sets flags for the resume protocol (BETA) and the backup flag. It also generates the `SMrandomTetradsRUNCHECK.h` header file which is used in the simulation code.

Finally, the script runs the simulation by executing the appropriate code. It is important to note that editing this part of the script is not recommended, as it may cause the simulation to not run as intended.

It is also important to note that running the script with the argument "man" will give you further information about the script.

This script is useful for automating the process of running simulations and it allows the user to easily change the simulation parameters without having to manually change them in the code. This also makes it easy to run multiple simulations with different parameters without having to manually change the code each time. This can save a lot of time and effort, especially when running large simulations that take a long time to complete.

It is also important to note that the script uses the parallel exchange algorithm, which allows the system to explore a wide range of temperatures in a relatively short amount of time. This can be useful for studying the equilibrium properties of the system over a wide range of temperatures.

Overall, the `run_simulation.sh` script is a useful tool for automating the process of running Monte Carlo simulations, and it allows the user to easily change the simulation parameters without having to manually change them in the code. This can save a lot of time and effort, especially when running large simulations that take a long time to complete.


The arch parameter.
-------------------

The first important thing that needs to be set is the architecture of your NVIDIA GPU. The -arch option of the nvcc command must match both the architecture and your version of CUDA. Here are dropped some examples:

#### CUDA 3.2 - 8.x (FERMI ARCH)

*   SM\_20, GeForce 400, 500, 600, GT-630.

#### CUDA 5.x - 10.x (KEPLER ARCH)

*   SM\_30, GeForce 700, GT-730.
*   SM\_35, Tesla K40.
*   SM37, Tesla K80.

#### CUDA 6.x - 11.x (MAXWELL ARCH)

*   SM\_50, Tesla/Quadro M series.
*   SM\_52, Quadro M6000 , GeForce 900, GTX-970, GTX-980, GTX Titan X.
*   SM\_53, Tegra (Jetson) TX1 / Tegra X1, Drive CX, Drive PX, Jetson Nano.

#### CUDA 8 and later (PASCAL ARCH)

*   SM\_60, Quadro GP100, Tesla P100, DGX-1 (Generic Pascal)
*   SM\_61, GTX 1080, GTX 1070, GTX 1060, GTX 1050, GTX 1030 (GP108), GT 1010 (GP108) Titan Xp, Tesla P40, Tesla P4, Discrete GPU on the NVIDIA Drive PX2
*   SM\_62, Integrated GPU on the NVIDIA Drive PX2, Tegra (Jetson) TX2

#### CUDA 9 and later (VOLTA ARCH)

*   SM\_70, DGX-1 with Volta, Tesla V100, GTX 1180 (GV104), Titan V, Quadro GV100
*   SM\_72, Jetson AGX Xavier, Drive AGX Pegasus, Xavier NX

#### CUDA 10 and later (TURING ARCH)

*   SM\_75, GTX/RTX Turing – GTX 1660 Ti, RTX 2060, RTX 2070, RTX 2080, Titan RTX, Quadro RTX 4000, Quadro RTX 5000, Quadro RTX 6000, Quadro RTX 8000, Quadro T1000/T2000, Tesla T4

#### CUDA 11.1 and later (AMPERE ARCH)

*   SM\_80, NVIDIA A100–GA100, NVIDIA DGX-A100
*   SM\_86, Tesla GA10x cards, RTX Ampere – RTX 3080, GA102 – RTX 3090, RTX A2000, A3000, RTX A4000, A5000, A6000, NVIDIA A40, GA106 – RTX 3060, GA104 – RTX 3070, GA107 – RTX 3050, RTX A10, RTX A16, RTX A40, A2 Tensor Core GPU
*   SM\_87, (from CUDA 11.4 onwards, introduced with PTX ISA 7.4 / Driver r470 and newer) for Jetson AGX Orin and Drive AGX Orin only

#### Lovelace (CUDA 11.8 and later)

*   SM\_89, NVIDIA GeForce RTX 4090, RTX 4080, RTX 6000, Tesla L40

Temperatures
------------

As is known, the critical temperature of the system depends on its size. This means that for each size, a range of temperatures should be chosen in order to make the transition to the glassy state visible. The maximum temperature for all sizes should be 1.6, while we report this table, which should be taken only as a rough reference, as based on the model's finishes (edge conditions, the way frequencies are distributed), they can also change significantly.

```
| size | t\_min |
| --- | --- |
| 18 | 0.40 |
| 32 | 0.45 |
| 48 | 0.50 |
| 62 | 0.55 |
| 76 | 0.70 |
| 96 | 0.80 |
| 120 | 0.80 |
| 150 | 0.85 |
```



Usage Guide for Analysis Class
==============================

This code is used to perform an analysis on a set of configurations and frequencies. It uses the following modules:

```python
import numpy as np
import os
import sys
import loadingModule 
import concurrent.futures as multithreading
```

The main class is called `Analysis`, and it should be initialized with the following parameters:

*   `param`: a list of 8 parameters that define the analysis
*   `path`: the path to the directory where the configurations and frequencies are stored

The `Analysis` class has the following methods:

*   `print_parameters()`: prints the parameters passed during initialization
*   `print_path()`: prints the path passed during initialization
*   `get_parameters()`: returns the parameters passed during initialization
*   `get_path()`: returns the path passed during initialization
*   `LoadWholeSample()`: loads all the configurations and stores them in the `configurations` attribute
*   `LoadFrequencies()`: loads the frequencies and stores them in the `frequencies` attribute
*   `print_frequencies()`: prints the frequencies
*   `print_config()`: prints the configurations

To use the code, create an instance of the `Analysis` class, passing the appropriate parameters and path, then call the appropriate methods to perform the analysis.

For example:

```scss
param = [10, 5, 0.1, 10000, 100, 100, 0, 100]
path = '/path/to/data'
analyzer = Analysis(param, path)
analyzer.LoadWholeSample()
analyzer.LoadFrequencies()
analyzer.print_frequencies()
```

This will initialize the analyzer with the given parameters and path, load the configurations and frequencies, and print the frequencies.






   
