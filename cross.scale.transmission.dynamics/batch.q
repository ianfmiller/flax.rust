#!/bin/bash
#
#SBATCH -o Rtest.out  # will be empty if there's no output
#SBATCH -e Rtest.err        # useful for troubleshooting, will be empty if there are no errors
#SBATCH -N 1
#SBATCH --cores-per-socket=10
#SBATCH --array=1-1
#SBATCH -t 0-12:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ifmiller@princeton.edu

srun R CMD BATCH simulate.epidemic.R