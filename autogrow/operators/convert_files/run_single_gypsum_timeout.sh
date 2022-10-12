#!/bin/bash

# This bash script will execute a timeout run of gypsum.
# This is necesary to be able to create a set of gypsum logs for each ligand.
# These logs are used to determine which ligands timed out and which ones converted.

# Arguments:
# executegypsum_path: the 1st uservar and is the path to the run_gypsum.py file. 
#           This is the Gypsum executable
# json_path: the 2nd uservar and is the path to the json file to submit to gypsum.
#           This provides gypsum with information like number of conformers to make, 
#               the SMILES string to convert, the number of processors (1), output path.
#
# This will echo Timeout if it timed out
# Echo completed the gypsum run if it completed the ligand conversion to 3D
# echo "FAILED ERRORS WITH THE LIGAND" if it failed conversion for any other reason.
executegypsum_path=$1
json_path=$2
timeout_limit=$3
python_path=$4

timeout $timeout_limit $python_path $executegypsum_path -j $json_path

# Test if it timed out. If it did it will append TIMEOUT to the end of the log file
exit_status=$?

echo "exit_status "$exit_status
if [ $exit_status -eq 124 ]; then 
    echo "TIMEOUT"
elif [ $exit_status -eq 0 ]; then 
    echo "Completed the gypsum run"  
elif [ $exit_status -eq 1 ]; then 
    echo "FAILED ERRORS WITH THE LIGAND"  
fi
