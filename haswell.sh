#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --constraint=haswell
srun -n 4 -c 8 shifter --env OMP_NUM_THREAD=8 --image=joezuntz/txpipe python shearz1z2.py

## check https://docs.nersc.gov/systems/cori/ for available nodes and cores
