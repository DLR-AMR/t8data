# t8code README gallery helper files

The t8code README (https://github.com/DLR-AMR/t8code/blob/ddc86e4e890b627e41da7187c8f0feca39bec855/README.md) showcases a .gif with application results from t8code.
In this folder we collect the raw input images and videos and scripts to generate the .gif file.

## build_slide.py

This script can be used to create a single slide of the .gif.
Given an input file and a caption, this creates an image or video file of the input together with the caption in the correct format.

Usage:
```bash
python3 build_slide.py [-h] [--subtitle SUBTITLE] [--fps FPS] [--quality QUALITY] input output title
```

## generate_readme_gallery.py

This script builds the final .gif from all present slides.
It scans the folder for image content and puts all together in a .gif.

Usage:
```bash
python3 generate_readme_gallery.py
```
