The [Dockerfile](https://bitbucket.org/berkeleylab/warp/src/master/Dockerfile) is used as a starting point to generate a Singularity recipe (`warp_python2.def and warp_python3.def`).
Warp is installed from [source](https://bitbucket.org/berkeleylab/warp/src/master) within the image.

## Dockerfile and Singularity recipe

The main differences in the Singularity recipe (apart from syntactic differences) w.r.t. the Dockerfile are that:
1. The Warp source is cloned at build time and Warp is installed in a system location (`/usr/local/warp`).
1. No user is created.
1. Both a `/home_mnt` and a `/lustre` directory are created in the image to allow the user to
   1. bind mount his/her $HOME (without conflicting with /home, problem cause by old kernel not supporting OverlayFS) and
   1. bind mount the Lustre FS from the host into the container

   respectively.

## Usage
One would do this on a node with a command line of this sort:

``singularity exec  -H $HOME:/home_mnt/ warp.simg command command_arguments``

Since Warp is installed in /usr/local/, to run a Warp example the user would execute it like this:

``singularity exec  -H $HOME:/home_mnt/ warp.simg python /usr/local/warp/examples/Solenoid_transport.py``

and the output files would be created on the current directory, which therefore needs to be bind mounted (so, in the example,
either within /lustre or $HOME or /tmp.

``singularity inspect warp.simg``

to get some information about the image.
