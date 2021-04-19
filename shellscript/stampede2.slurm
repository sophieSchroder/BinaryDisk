#!/bin/bash
#SBATCH -N 8
#SBATCH -n 512
#SBATCH -t 40:00:00
#SBATCH -p normal
#SBATCH -A TG-AST150042
#SBATCH -J g11hllen


module purge
module load intel
module load impi
module load phdf5

ibrun ./code/bin/athena -i athinput.binarydisk_stream
