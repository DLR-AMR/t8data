#!/bin/bash

# This is a batch-file used on CARA. Cara has 2168 nodes, a 2x(AMD EPYC (32 Cores))


#SBATCH --time=00:10:00
#SBATCH --ntasks=128
#SBATCH --output=Mesh_%j
#SBATCH --error=Mesh_err_%j
NUM_PROCS="128"

# d = Dimension
# l = initial level
# r = number of refinements
# x = left wall of the mesh
# t = thickness of the refinement region
# D = distance of the refinement region to travel
# s = number of steps to travel
# n = number of reruns
# b = balance after each refinement step
# g = ghost cells

PART_ARGS="-d3 -l5 -r2 -x-0.5 -t0.2 -D2 -s5 -n2 -g"


JOBFILE="/scratch_fast/ws/0/knap_da-t8code_benchmarks/benchmark_build/benchmark"

MSH_FILE="/scratch_fast/ws/0/knap_da-t8code_benchmarks/t8data/paper_files/New_for_hybrid/tonne_100k"

JOB_CMD="$JOBFILE -f $MSH_FILE $PART_ARGS"
for PROCS in $NUM_PROCS ; do
    echo "-------------  Running: $JOB_CMD with $PROCS procs ------------"
    srun -n $PROCS $JOB_CMD &
done

