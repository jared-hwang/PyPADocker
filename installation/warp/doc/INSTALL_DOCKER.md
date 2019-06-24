# Running Warp with Docker

This document describes how to run Warp within a Docker container. [Docker](https://www.docker.com/) is a framework that allows to run software in a portable manner on different plateforms.

## Downloading and installing Docker

Docker can be installed from [here](https://www.docker.com/products/overview).

## Preparing the Docker image of Warp

The Docker image of Warp (which contains a compiled version of Warp and all of its dependencies), can be obtained by downloading it from [DockerHub](https://hub.docker.com/).

To download the image, type
```
docker pull rlehe/warp
```

### Building the Warp image locally

To build the image locally, clone the Warp repository and build the image using the provided `Dockerfile`:
```
git clone https://bitbucket.org/berkeleylab/warp.git
cd warp
docker build -t warp .
```

## Running simulations

In order to run a simulation, create a new directory,
copy your Warp input script to this directory, and rename this script
to `warp_script.py`. (The folder `scripts/examples/` of the
[Warp repository](https://bitbucket.org/berkeleylab/warp/src) contains
several examples of input scripts.)

Then `cd` into this directory, and type the following command to enter the Docker container.
```
docker run -it -v $PWD:/home/warp_user/run rlehe/warp
```

Then launch the simulation by typing either (for serial simulations)
```
python warp_script.py
```
or (for e.g. a parallel simulation with a 2x3x2 domain decomposition in 3D)
```
mpirun -np 12 python warp_script.py -p 2 3 2
```
