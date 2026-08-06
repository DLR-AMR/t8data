Compute Pi example

With https://github.com/DLR-AMR/t8code/pull/1800 we added an example that approximates Pi to t8code.
On Pi approximation day 2025, i.e. July 22nd, we posted a video to LinkedIn showcasing this example.
We approximate Pi by computing the area of the unit circle. To do so we iteratively refine at the boundary of the unit circle and 
compute the area of all triangles inside of it.

<p align="center">
  <img width="1000px" src=ComputePi.jpeg>
</p>

This data repo consists of:

- The screenshot above: https://github.com/DLR-AMR/t8data/blob/t8ComputePi/pictures/compute_pi/ComputePi.jpeg
- The video that we posted: https://github.com/DLR-AMR/t8data/tree/t8ComputePi/pictures/compute_pi/video
- A paraview state from which we created the screenshot and video: https://github.com/DLR-AMR/t8data/blob/t8ComputePi/pictures/compute_pi/compute_pi_pvstate.pvsm
- The paraview "pvtu" and "vtu" files that we used for the video. The paraview state uses these files: https://github.com/DLR-AMR/t8data/tree/t8ComputePi/pictures/compute_pi/vtu
