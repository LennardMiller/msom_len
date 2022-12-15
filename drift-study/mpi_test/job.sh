##!/bin/bash            
##OAR -l /nodes=4/cores=64
##OAR --project pr-data-ocean

## Ensure Nix is loaded. The following line can be into your ~/.bashrc file.
source /applis/site/nix.sh

## Run the program
mpirun -np `cat $OAR_FILE_NODES|wc -l` --machinefile $OAR_NODE_FILE -mca plm_rsh_agent "oarsh" --prefix $HOME/.nix-profile $OAR_WORKING_DIRECTORY/qg.e
