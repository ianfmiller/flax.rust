#!/bin/bash
#
#SBATCH -N 1
#SBATCH -J "flax.epi"
#SBATCH --cores-per-socket=10
#SBATCH -t 0-12:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ifmiller@princeton.edu

srun simulate.epidemic.R
