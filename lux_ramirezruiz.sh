#!/bin/bash
#SBATCH --job-name=q03      # Job name
#SBATCH --partition=ramirez-ruiz             # queue for job submission
#SBATCH --account=ramirez-ruiz               # queue for job submission
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sschrode@ucsc.edu   # Where to send mail
#SBATCH --ntasks=160                  # Number of MPI ranks
#SBATCH --nodes=4                    # Number of nodes
#SBATCH --ntasks-per-node=40         # How many tasks on each node
#SBATCH --time=6-23:59:00              # Time limit hrs:min:sec
#SBATCH --output=test_%j.log     # Standard output and error log

pwd; hostname; date

echo "Running program on $SLURM_JOB_NUM_NODES nodes with $SLURM_NTASKS total tasks, with each node getting $SLURM_NTASKS_PER_NODE running on cores."

module load intel
module load intel/impi
module load hdf5/1.10.6-parallel



mpirun -n 160 --ppn 40 /home/sschrode/Athena/DISK/BinaryDisk/code/bin/athena -i athinput.binarydisk_stream

date
