# README

## How to run the benchmarks
Install the benchmarks via CMake and provide Paths to t8code (and SC & P4est, which can be installed with t8code alltogether)
Run the benchmark_elems example to evaluate the element-wise performance. With -n you can define the number of repeated runs, with -e you define the type of element used and with -l the initial uniform refinement level of the forest. 

Run the benchmark example to evaluate the performance on a large hybrid mesh. Provide a mesh file via -f, and describe the dimension of the mesh via -d. The initial uniform refinement level is given via -l and -r defines the number of refinement levels. -g toggles on the usage of ghost-cells, -b enables the computation of a 2:1 balancing. -x is the minimum coordinate of the mesh, -t describes the thickness of the wall and -D the distance of the wall to travel. -s defines the step the wall should make. -n is number of repeated runs to test. 
We used "-d3 -l5 -r2 -x-0.5 -t0.2 -D2 -s5 -n2 -g" for our tests.

To run the examples on a cluster you can use the scripts provided in the jobs directory.

