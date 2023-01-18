#!/bin/bash
#This script checks if all sample directories have the same number of files. 

N=10 # replace with desired number of directories

# initialize variable to store number of elements in first directory
elements_in_first_dir=$(ls -1 sample1 | wc -l)

# loop through remaining directories and check if they have the same number of elements
for i in $(seq 2 $N); do
  elements_in_dir=$(ls -1 sample$i | wc -l)
  if [ $elements_in_first_dir -ne $elements_in_dir ]; then
    echo "sample$i has a different number of elements"
    
  fi
done

echo "All directories have the same number of elements"
