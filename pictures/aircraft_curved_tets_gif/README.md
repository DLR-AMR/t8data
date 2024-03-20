This folder contain all the data used for the generation of the aircraft_curved_tets_gif.

The vtk files are generated with our tutorials/features/t8_features_curved_meshes.cxx file with the following input:

./t8_features_curved_meshes -f aircraft -d3 -l0 -p -x-5 -X40 -r2 -t10 -n30 -o

The aircraft.brep and aircraft.msh files contain the geometry and the mesh.
The visualization state in Paraview is saved in the aircraft_gif_state.pvsm file.
Each frame can be exported with the "Export Animation" option in Paraview (png export works best for the gif generation).
The quality of the export can be set manually.
With the make_gif_pillow.py script it is possible to generate a GIF out of the exported pictures.
Usage of make_gif_pillow.py:

python make_gif_pillow.py <string included in all image filenames> !<string which should be omited>