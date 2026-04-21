<img width="1024" height="1001" alt="grafik" src="https://github.com/user-attachments/assets/9ee7a2bf-888c-477b-969d-1f3832363826" />

These files can be used to reconstruct our welcome post for Benedict and Michele, which is a video of an adapting mesh around their names.
We created it with Powerpoint and png2mesh.

You find the scripts to replicate the welcome video in this folder and in the [scripts folder](https://github.com/DLR-AMR/t8data/tree/linkein_welcome_post_2026/scripts).

```bash
split_mp4_into_png.scp WelcomeBenedictMichele.mp4
run_png2mesh_on_frames.scp Welcome 216
```

- Open paraview and load the state "Welcome_paraview.pvsm"
- Export to Animation.

Paraview might not be able to export it as video, but only as individual png files.
In that case, we merge them again using:

```bash
merge_png_to_mp4.scp Welcome_tri Welcome_tri.mp4
./optimize_mp4_for_linkedin.scp Welcome_tri.mp4 Welcome_tri_linkedin.mp4
```
