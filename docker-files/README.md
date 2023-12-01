# Dockerfiles

This folder contains dockerfiles for t8code configurations.
There are also [Docker images](https://hub.docker.com/r/dlramr/t8code-ubuntu/tags) listed in the Docker hub.
Those Docker files can be accessed in the [t8code-docker-images](https://github.com/DLR-AMR/t8code-docker-images) repository.

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
