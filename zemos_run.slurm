#!/bin/bash
#-----------------------------------------------------------------
# Example SLURM job script adapted from MOGON. ZEMOS has 2
# partition,QC and MD. QC has 40 dual socket, 16 core/nodes,
# 84GBRam, while MD has 188 dual socket, 24 core/node,192GB Ram.
# It is important to consider -n and -c does not exceed the
# amount of core/nodes
#
#-----------------------------------------------------------------
#SBATCH -J POST       # Job name
#SBATCH -p MD                   # Partition name QC or MD
#SBATCH -N 1                   # Total number of nodes
#SBATCH -n 24                    # Total number of tasks
#SBATCH -c 1                    # Total cpus required per task
#SBATCH --propagate=STACK


# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1

# Load all necessary modules if needed (these are examples)
# Loading modules in the script ensures a consistent environment.
module load mpi/openmpi/gnu/4.11
module load  python/intelpython3/2022.1.2
export OMP_STACKSIZE=10000M

# Launch the executable

run_types=( "spce" "tip3p" "pure_spce" "pure_tip3p" )

for run_type in ${run_types[@]}; do
    parent_directory="/data/esmith/water_models/$run_type"
    for run in `ls $parent_directory`; do
        data_file="$parent_directory/$run/outdata/merged_posteq.data"
        xyz_file="$parent_directory/$run/traj/merged_prod.xyz"

        python3 zemos_tetra.py $xyz_file $data_file "./plots/$run"&
    done
done
