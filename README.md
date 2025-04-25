# t8data

This repository manages large data associated to [t8code](https://github.com/holke/t8code) and is currenctly under development.
The data is managed via [git lfs](https://git-lfs.github.com/), please read the documentation of git lfs befor contributing. 

If you don't want to download all data when you clone the repository (especially if you want to use a benchmark on a cluster and don't want do download all pictures, meshfiles, gifs, videos, etc ...) you can run
`export GIT_LFS_SKIP_SMUDGE=1` befor cloning. This only sets an environment variable. If you want to ignor LFS files permanently you have to set the variable in you git config. 

You will find examples of large meshes managed by t8code, the raw data from benchmark runs used for publishing and much more data used by or comming from t8code. 
