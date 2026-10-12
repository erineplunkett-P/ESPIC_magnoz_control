#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=384G
#SBATCH --time=168:00:00
#SBATCH --mail-user=sra2ja@virginia.edu
#SBATCH --mail-type=all
#SBATCH --partition=standard
#SBATCH --constraint=afton
#SBATCH --account=CHEN_CUI
#SBATCH --output=out_ESPIC_debug_%j.out
#SBATCH --error=out_ESPIC_debug_%j.err

# Avoid nounset because Intel setvars.sh references unset vars.
# set -eo pipefail

WORK_HOME=/scratch/sra2ja/20251230_H_0857_VDF_large_continue
cd $WORK_HOME

module purge

module load gcc
module load openmpi
module load eigen
module load fftw

source /apps/software/standard/mpi/gcc/11.4.0/openmpi/4.1.4/imkl/2023.1.0/setvars.sh

ulimit -c unlimited
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd src/
make clean
make

cd ../
cd run/
sh run.sh 

# srun -n $SLURM_NTASKS  ./empic
# mpirun -n $SLURM_NTASKS ./empic
