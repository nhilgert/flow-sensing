#!/bin/sh -l
# FILENAME: hilgertsim
#PBS -V
#PBS -q amugler
#PBS -l nodes=1:ppn=1,naccesspolicy=singleuser
#PBS -l walltime=00:05:00
#PBS -N hilgertsim

cd $PBS_O_WORKDIR

./partsim $seed $T
