SMrandomTetrads
===============

SMrandomTetrads is a program written by many people, and has a history of almost ten years. The goal of this program is to simulate the 4-phasor random laser model, which is known to be a disordered spin glass system. Here we explain how to use this tool.

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

#) SEC-5 The directory tree





   
