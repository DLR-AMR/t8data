# Dockerfiles

This folder contains dockerfiles for t8code configurations.

# Building

To build a container (for example after changing the dockerfile):

```
docker build -t [BUILD-NAME] /path/to/dockerfile
```

# Running

To run a container interactively use:

```
docker run -it [BUILD-NAME]
```



# Currently available containers:

## ubuntu-22_04

Contains a debug and release configuration of t8code on Ubuntu 22.04.
See `/exec/t8code_debug`  and `/exec/t8code_release`.
