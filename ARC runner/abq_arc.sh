#!/bin/bash

#SBATCH --clusters=arc
#SBATCH --partition=devel
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --time=00:10:00
#SBATCH --job-name=AbaqusJob

module purge

module load Abaqus/2022
module load iimpi/2020a


. abaqus.sh


abaqus input=Job-1.inp job=AbaqusJob user=DBFcode.f cpus=${SLURM_NTASKS} interactive
