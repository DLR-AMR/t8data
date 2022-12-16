The christmas card 2022 was created using png2mesh, commit 3cc10d79ff0a37a818dc98be2f2b778e16af0e92, linked against t8code with subelements, using the following call:


mpirun -np 2 ./png2mesh_demo -f t8codeWeihnachtsgruesse2022.png -m12 -e3


By loading the paraview state t8codeWeihnachtsgruesse2022.pvsm you can recreate the screenshot.
Note, that this requires you to build the pvtu and vtu files using the above png2mesh call.
