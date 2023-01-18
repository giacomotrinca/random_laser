
SMrandomTetrads is a program written by many people, and has a history of 
almost ten years. The goal of this program is to simulate the 4-phasor random 
laser model, which is known to be a disordered spin glass system.
Here we explain how to use this tool.

#) SEC-1: The arch parameter.
The first important thing that needs to be set is the architecture of your NVIDIA GPU. 
The -arch option of the nvcc command must match both the architecture and your version of CUDA.
Here are dropped some examples:
#-------------------------------------------------------------------#
° CUDA 3.2 - 8.x (FERMI ARCH)
#SM_20,
    GeForce 400, 500, 600, GT-630.
#-------------------------------------------------------------------#
° CUDA 5.x - 10.x (KEPLER ARCH)
#SM_30, 
    GeForce 700, GT-730.
#SM_35, 
    Tesla K40.
#SM37, 
    Tesla K80.
#-------------------------------------------------------------------#
° CUDA 6.x - 11.x (MAXWELL ARCH)
#SM_50,
    Tesla/Quadro M series.
#SM_52,
    Quadro M6000 , GeForce 900, GTX-970, GTX-980, GTX Titan X.
#SM_53,
    Tegra (Jetson) TX1 / Tegra X1, Drive CX, Drive PX, Jetson Nano.
#-------------------------------------------------------------------#
° CUDA 8 and later (PASCAL ARCH)
#SM_60,
    Quadro GP100, Tesla P100, DGX-1 (Generic Pascal)
#SM_61,
    GTX 1080, GTX 1070, GTX 1060, GTX 1050, GTX 1030 (GP108), 
    GT 1010 (GP108) Titan Xp, Tesla P40, Tesla P4, Discrete GPU on the NVIDIA Drive PX2
#SM_62, 
    Integrated GPU on the NVIDIA Drive PX2, Tegra (Jetson) TX2
#-------------------------------------------------------------------#
° CUDA 9 and later (VOLTA ARCH)
#SM_70,
    DGX-1 with Volta, Tesla V100, GTX 1180 (GV104), Titan V, Quadro GV100
#SM_72,
    Jetson AGX Xavier, Drive AGX Pegasus, Xavier NX
#-------------------------------------------------------------------#
° CUDA 10 and later (TURING ARCH)
#SM_75,
    GTX/RTX Turing – GTX 1660 Ti, RTX 2060, RTX 2070, RTX 2080, Titan RTX, 
    Quadro RTX 4000, Quadro RTX 5000, Quadro RTX 6000, Quadro RTX 8000, 
    Quadro T1000/T2000, Tesla T4 
#-------------------------------------------------------------------#
° CUDA 11.1 and later (AMPERE ARCH)
#SM_80,
    NVIDIA A100–GA100, NVIDIA DGX-A100
#SM_86,
    Tesla GA10x cards, RTX Ampere – RTX 3080, GA102 – RTX 3090, RTX A2000, A3000, 
    RTX A4000, A5000, A6000, NVIDIA A40, GA106 – RTX 3060, GA104 – 
    RTX 3070, GA107 – RTX 3050, RTX A10, RTX A16, RTX A40, A2 Tensor Core GPU
#SM_87, 
    (from CUDA 11.4 onwards, introduced with PTX ISA 7.4 / Driver r470 and newer)
    for Jetson AGX Orin and Drive AGX Orin only

#-------------------------------------------------------------------#
Lovelace (CUDA 11.8 and later)
#SM_89,
    NVIDIA GeForce RTX 4090, RTX 4080, RTX 6000, Tesla L40
#-------------------------------------------------------------------#
Hopper (CUDA 12 and later)
#SM_90,
    NVIDIA H100 (GH100)
#-------------------------------------------------------------------#

For further informations visit 
https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/
#) SEC-2: Simulation settings.
#-------------------------------------------------------------------#
° Size
  This value represent the random laser's modes number. It must be an integer. The same value must be 
  set in the right section of simulation_code/SMrandomTetradsRUNCHECK.h
#-------------------------------------------------------------------#
° Temperature

  T_max might be 1.6 for all sizes, and here we drop some reference values for t_min
  # -------------- #
  # size  |  t_min #
  # -------------- #
  # 18    |   0.40 #
  # 32    |   0.45 #
  # 48    |   0.50 #
  # 62    |   0.55 #
  # 76    |   0.70 #
  # 96    |   0.80 #
  # 120   |   0.80 #
  # 150   |   0.85 #
  # -------------- #
  
  These reference values are acceptable for the model corresponding to the "comb" frequencies. 
  By choosing these temperatures, the specific heat graph will display the peak, which corresponds 
  to the critical temperature for that specific sample.
#-------------------------------------------------------------------#
° Resume a simulation (BETA)
  Sometimes, simulations can crash. If you don't want to start the simulation from scratch 
  for that specific sample, you can try restarting the Monte Carlo simulation from the last saved step. 
  To do this, simply set resume=1 in the "Simulation settings" section in the run_simulation.sh script.

#) SEC-3 Run simulations and checks
To launch the program, once all the necessary libraries and the nvcc compiler 
are installed, just run the run_simulation.sh script with:

$ ./run_simulation.sh index_of_GPU initial_sample number_of_samples

where:
° index_of_GPU = label of GPU device in which we want to perform the simulation
° initial_sample = index of the starting sample
° number_of_samples = number of independent MonteCarlo simulation we want to 
  perform for each GPU device. 


#) SEC-4 Run analysis (BETA)

#) SEC-5 The directory tree


   
