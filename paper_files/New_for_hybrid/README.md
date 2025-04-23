# README

## Comparison of the Partitioning Algorithm

To compare the timings between the introduction of the New-for-hybrid algorithm link with
- t8code v4.0.0 for the "old" performance
- t8code add tag for the new performance

## What is benchmarked?

In benchmark.cxx we create a cmesh, either from a file or from our examples. We use a hybrid mesh from the example, consisting of a hexahedron, a tetrahedron, a prism and a pyramid. For the published run-times we used a cmesh created from xy.msh.

The mesh is then adaptively refined, coarsened and repartitioned for multiple timesteps. In each timestep elements inside a wall are refined. If an element is outside the wall it is coarsened up until a minimum level. In each timestep the wall is repositioned, enforcing a repartitioning of the forest in each timestep.
