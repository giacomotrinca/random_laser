#!/bin/bash


#---------------------------------------------#
#             GPU architecture                #
#---------------------------------------------#
#         BEFORE RUN CHECK THIS VALUE 	      #
#	It must be the right one, according   #
#	to the README.md file, SEC-I          #
#---------------------------------------------#
# 30 (VISNU/DURGA) 35 (KRAKEN) 70(ECATON)
arch=86

#---------------------------------------------#
#             PATH FOR GSL LIBRARY            #
#---------------------------------------------# 
gsl_path=/usr/local/include/gsl #for my laptop
#gsl_path=/usr/include/gsl	#for cluster
#---------------------------------------------#
#              Simulation settings            #
#---------------------------------------------#
Size=48
t_min=0.40
t_max=1.1
number_of_PT_replicas=20
number_of_real_replicas=10
PT_flag=1
number_of_PT_step=64
power_of_iterations=20      #the real number of iteration will be (2^number_of_iteration)/number_of_PT_step
number_of_iterations=$((2**$power_of_iterations/$number_of_PT_step))
print_config=4		                            #configurations will be printed every 
						    #$print_config * $number_of_PT_step steps
						    
#iter_start_print=$number_of_iterations/2	    #typical values= 0, 
						    #$number_of_iterations/2, $number_of_iteration/4
			                            #if PT_flag = 1 then iter_start_print should be 
			                            #different from zero
if [ "$PT_flag" == "1" ]; then
	#echo "Warning: you have chosen PT_flag=1 but you are not printing all configurations..."
	#echo "Setting iter_start_print = 0"
	iter_start_print=0
elif [ "$PT_flag" == "0" ]; then
	iter_start_print=$number_of_iterations/2
else 
	echo "PT_flag has an invalid value. Exiting."
	exit 1
fi
		                            
print_all=0
frequency_bandwidth=1.
frequency_mode=3	    #if its 1 -> combliek frequencies
			    #if it's 2 -> list of frequencies loaded from file (to be implemented)
			   #if it's 3 -> random frequencies taken from 0 to $frequency_bandwidth
ntjumps=$((number_of_PT_replicas -1))
#flags for resume protocol (BETA)
backup_flag=0
print_loading=0


#printing parameters settings in a file for the analysis input
echo "#This file was automatically generated by run_simulation.sh script"

echo "$Size" > analysis_input.info
echo "$number_of_real_replicas" >> analysis_input.info
echo "$number_of_PT_step" >> analysis_input.info
echo "$number_of_iterations" >> analysis_input.info
echo "$iter_start_print" >> analysis_input.info
echo "$number_of_PT_replicas" >> analysis_input.info
echo "$PT_flag" >> analysis_input.info

#printing the SMrandomTetradsRUNCHECK.h header file
echo "//This file was automatically generated by run_simulation.sh script" > simulation_code/SMrandomTetradsRUNCHECK.h

echo "#define Size $Size" >> simulation_code/SMrandomTetradsRUNCHECK.h
echo "#define NREPLICAS	$number_of_real_replicas" >> simulation_code/SMrandomTetradsRUNCHECK.h
echo "#define NSTEP $number_of_PT_step" >> simulation_code/SMrandomTetradsRUNCHECK.h
echo "#define NITERATIONS $number_of_iterations" >> simulation_code/SMrandomTetradsRUNCHECK.h
echo "#define NITER_PRINT_CONF $print_config" >> simulation_code/SMrandomTetradsRUNCHECK.h
echo "#define NITER_MIN_PRINT $iter_start_print" >> simulation_code/SMrandomTetradsRUNCHECK.h
echo "#define PRINTALLCONFIG $print_all" >> simulation_code/SMrandomTetradsRUNCHECK.h
echo "#define NPT $number_of_PT_replicas" >> simulation_code/SMrandomTetradsRUNCHECK.h
echo "#define NTJUMPS $ntjumps" >> simulation_code/SMrandomTetradsRUNCHECK.h
echo "#define PT_EXCHANGE $PT_flag" >> simulation_code/SMrandomTetradsRUNCHECK.h
echo "#define DW $frequency_bandwidth" >> simulation_code/SMrandomTetradsRUNCHECK.h
echo "#define FREQ_ENABLE $frequency_mode" >> simulation_code/SMrandomTetradsRUNCHECK.h
echo "#define BACKUP $backup_flag" >> simulation_code/SMrandomTetradsRUNCHECK.h
echo "#define PRINTLOADING $print_loading" >> simulation_code/SMrandomTetradsRUNCHECK.h






#flag for resuming a simulation (BETA)
resume=0    #check the last sample by looking 
            # at simulation.info file








# DO NOT EDIT THIS PART OF THE SCRIPT
#---------------------------------------------#
#             some useful argument            #
#---------------------------------------------#
if [ "$1" == "-v" ] || [ "$1" == "--version" ]; then
    clear
    echo "SMrandomTetrads - v 0.1.0-ResumeBeta1"
    echo "SMrandomTetrads - Analysis - v 3.0.0"
    exit 0
fi

if [ "$1" == "man" ]; then
	clear
	more README.md
fi

if [ "$1" == "--help" ] || [ "$1" == "-h" ] || [ -z "$1" ]; then
	echo -e "\033[31mUsage\033[0m"
	echo -e "\033[31m./run_simulation.sh index_of_GPU initial_sample number_of_samples\033[0m"
	#echo "./run_simulation.sh index_of_GPU initial_sample number_of_samples"
	echo "List of argument accepted:"
	echo " "
	echo "--version or -v"
	echo "Print the current version of the simulation and analysis programs"
	echo "man"
	echo "Print the whole instructions (README.md file)"
	echo "--help or -h"
	echo "Print this help page"
	exit 0
fi

#---------------------------------------------#
#          checking the right syntax          #
#---------------------------------------------#

if [ -z "$1" ]; then
    echo "You must insert the index of GPU device, and it must be an integer"
    echo "Usage:"
    echo "./run_simulation.sh index_of_GPU initial_sample number_of_samples"
    exit 1
else
    if [ -z "$2" ]; then
    	echo "You must insert the index of starting sample, and it must be an integer"
    	echo "Usage:"
    	echo "./run_simulation.sh index_of_GPU initial_sample number_of_samples"
    	exit 1
    else
    	if [ -z "$3" ]; then
    		echo "You must insert the number of samples you want simulate, and it must be an integer"
    		echo "Usage:"
    		echo "./run_simulation.sh index_of_GPU initial_sample number_of_samples"
    		exit 1
    	fi
    fi
fi


clear
echo "##############################################################################"
echo "#  °°°°°°°      °°°°°°°°     °°°°°°°°   °°°°°°°°°°°°°°      °°°°°°°°°°°°°°°° #"
echo "# °°°°          °°°°°°°°°   °°°°°°°°°   °°°°°      °°°°     °°°°°°°°°°°°°°°° #"
echo "# °°°°          °°°°°°  °° °°  °°°°°°   °°°°°°°°°°°°°            °°°°°°      #"
echo "#    °°°°°      °°°°°°    °    °°°°°°   °°°°°°°°°°°              °°°°°°      #"
echo "#       °°°°°   °°°°°°         °°°°°°   °°°°°      °°°°°         °°°°°°      #"
echo "#     °°°°°°    °°°°°°         °°°°°°   °°°°°      °°°°°         °°°°°°      #"
echo "# °°°°°°°       °°°°°°         °°°°°°   °°°°°       °°°°         °°°°°°      #"
echo "##############################################################################"




#---------------------------------------------#
#   variables given when you run the script   #
#---------------------------------------------#

#first argument of the script call
device=$1 
#second argument of the script call
sample_init=$2
#third argument of the script call                      
number_of_samples=$3



#computing the number of last sample
sample_end=$[ $number_of_samples + $sample_init -1 ]



#---------------------------------------------#
# creating the directory  tree for simulation #
#---------------------------------------------#

mkdir -p N$Size
mkdir -p N$Size/DEV$device
if [ "$PT_flag" == "1" ]; then
	path=N$Size/DEV$device/$(date +%Y%m%d_%H%M%S)_PT
elif [ "$PT_flag" == "0" ]; then
	path=N$Size/DEV$device/$(date +%Y%m%d_%H%M%S)
else 
	echo "the value chosen for pt_flag = $PT_flag is invalid."
	exit 1
fi
   
mkdir -p $path
#---------------------------------------------#

echo "The simulation started at $(date)" > $path/simulation.info

for i in $(seq $sample_init $sample_end); do
    echo "Compiling source..."
    mkdir -p $path/sample$i
    
    init_time=$(date +%s)

    nvcc -I $gsl_path -lcurand -arch=compute_$arch -use_fast_math -w -m=64 -O3 ./simulation_code/SMrandomTetrads.cu -o SMrandomTetrads_sample$i

    mv SMrandomTetrads_sample$i $path/sample$i
    cd $path/sample$i
    
    random_seed=$RANDOM
    
    ./SMrandomTetrads_sample$i $random_seed $t_min $t_max $device $resume
    resume=0
    
    end_time=$(date +%s)
    time_elapsed=$((end_time - init_time))
    num_files=$(ls -1 | wc -l)
    cd ../../../..
    time_string=$(printf "%02dh %02dm %02ds" $(($time_elapsed/3600)) $(($time_elapsed%3600/60)) $(($time_elapsed%60)))
    echo "Sample $i (seed $random_seed) finished at $(date) ($time_string) and sample$i contains $num_files files" >> $path/simulation.info
done

echo "The simulation finished at $(date)" >> $path/simulation.info
