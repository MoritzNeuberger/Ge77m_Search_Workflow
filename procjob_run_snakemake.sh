#!/bin/bash
#SBATCH --nodes=1
#SBATCH --qos=shared
#SBATCH --time=08:00:00
#SBATCH --constraint=cpu
#SBATCH --account=m2676
#SBATCH --export=HDF5_USE_FILE_LOCKING=FALSE
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=moritz.neuberger@tum.de
#SBATCH --chdir=/global/cfs/cdirs/m2676/users/neuberger/Ge77m_dc_search/v04/Ge77m_Search_Workflow
#SBATCH --output=/global/cfs/cdirs/m2676/users/neuberger/Ge77m_dc_search/v04/Ge77m_Search_Workflow/.log/perlmutter-%j.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30GB

echo "Job Start:"
date
echo "Node(s):  "$SLURM_JOB_NODELIST
echo "Job ID:  "$SLURM_JOB_ID
echo "Ncpu:  "$SLURM_CPUS_PER_TASK
scontrol show job $SLURM_JOB_ID | grep -i cpus

export TERM=xterm-256color

module load python
source activate ge77m_snakemake

cd /global/cfs/cdirs/m2676/users/neuberger/Ge77m_dc_search/v04/Ge77m_Search_Workflow

snakemake -c 10 --unlock
snakemake -c 10
