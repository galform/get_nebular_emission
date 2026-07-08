#!/bin/sh
# Template for SLURM job submission - GNE processing

#SBATCH -A 16cores
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=__GNE_JOB_NAME__
#SBATCH --error=__GNE_LOG_DIR__/%x_ivol%a_%A.err
#SBATCH --output=__GNE_LOG_DIR__/%x_ivol%a_%A.out
#SBATCH --partition=all
#SBATCH --exclude=epi
#SBATCH --time=05:00:00
#SBATCH --array=__GNE_VOLS__%30
#
export GNE_SUBVOL_INDEX=$SLURM_ARRAY_TASK_ID
srun python << 'EOF_PYTHON_SCRIPT'
__GNE_PARAM_CONTENT__
EOF_PYTHON_SCRIPT
