#!/bin/bash

# This is a batch-file used on CARA. Cara has 2168 nodes, a 2x(AMD EPYC (32 Cores))


#SBATCH --time=00:10:00
#SBATCH --ntasks=64
#SBATCH --output=New_for_hybrid_PYRA_%j
#SBTACH --error=New_for_hybrid_PYRA_err_%j
NUM_PROCS="1 2 4 8 16 32 64"

# e = Element typ
# l = initial level
# r = number of refinements
# n = number of reruns
# x = left wall
# X = right wall
# T = time
# C = CFL
ARGS="-e3 -l3 -r3 -n3 -x-0.1 -X0 -T1 -C1"

JOBFILE="/scratch/ws/4/knap_da-t8code_timings/benchmark/benchmark"

for PROCS in $NUM_PROCS ; do
	JOB_CMD="$JOBFILE $ARGS"
	echo "-------------  Running: $JOB_CMD with $PROCS procs ------------"
	srun -n $PROCS $JOB_CMD
done
