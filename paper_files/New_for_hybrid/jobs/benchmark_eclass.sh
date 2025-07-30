#!/bin/bash

# This is a batch-file used on CARA. Cara has 2168 nodes, a 2x(AMD EPYC (32 Cores))


#SBATCH --time=00:10:00
#SBATCH --ntasks=64
#SBATCH --output=New_for_hybrid_ELEM_%j
#SBATCH --error=New_for_hybrid_ELEM_err_%j
NUM_PROCS="8 16 32 64"
ELEMENT_NAMES=("TETRAHEDRON" "HEXAHEDRON" "PRISM" "PYRAMID")

# e = Element typ
# l = initial level
# n = number of reruns

PART_ARGS="-l3 -n3"


JOBFILE="/scratch/ws/4/knap_da-t8code_timings/benchmark/benchmark"

for i in {0..3}; do
	ELEMENT=${ELEMENT_NAMES[$i]}
	echo "TEST $ELEMENT"
	ARGS="$PART_ARGS -e$i"
	JOB_CMD="$JOBFILE $ARGS"
	for PROCS in $NUM_PROCS ; do
		echo "-------------  Running: $JOB_CMD with $PROCS procs ------------"
		srun -n $PROCS $JOB_CMD &
	done
done

