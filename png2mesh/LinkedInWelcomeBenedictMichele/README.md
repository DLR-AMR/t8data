To replicate the welcome video:

split_mp4_into_png.scp WelcomeBenedictMichele.mp4
run_png2mesh_on_frames.scp Welcome 216

Open paraview and load the state "Welcome_paraview.pvsm"
Export to Animation.

Paraview might not be able to export it as video, but only as individual png files.
In that case, we merge them again using:

merge_png_to_mp4.scp Welcome_tri Welcome_tri.mp4
./optimize_mp4_for_linkedin.scp Welcome_tri.mp4 Welcome_tri_linkedin.mp4
