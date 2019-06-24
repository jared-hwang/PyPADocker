# Installing Warp with mpi4py on Edison

This document describes how to install Warp on the Edison cluster.

## Setting up the environnement

Create a directory that will contain all the software related to Warp:
```
mkdir $SCRATCH/warp_install
```
NB: This directory is on the `$SCRATCH` directory, for fast reading
access. Also, the `$SCRATCH` directory is specific to each cluster
(it is a different folder for Hopper, Edison and Cori) and is thus
well-adapted for system-specific installation, whereas the
`$HOME` directory is shared for all the clusters.

Then enter the following lines in the file `.bashrc.ext` that resides in the `$HOME` directory. (This automatically loads a few modules used by Warp, whenever a simulation is run.)
```
if [ "$NERSC_HOST" == "edison" ] 
then
   module swap PrgEnv-intel PrgEnv-gnu
   module load python/2.7-anaconda
   module load h5py-parallel
   export PATH=$SCRATCH/warp_install/bin:$PATH
   export PYTHONPATH=$SCRATCH/warp_install/lib/python:$PYTHONPATH
fi
```

Then run `source .bashrc` in the `$HOME`directory (or alternatively
log out of Edison, and then log in again).

## Installing Forthon

As explained on the Warp website, in order to install Forthon :

- Type `git clone https://github.com/dpgrote/Forthon.git`
- Then `cd` into the directory `Forthon` and run `python setup.py install --home=$SCRATCH/warp_install`

## Installing Warp itself

### In order to install warp with gfortran

- Download the source by entering `git clone https://bitbucket.org/berkeleylab/warp.git`

- Go to the `warp/pywarp90` directory and create a file called `Makefile.local.pympi`. Enter the following lines in this file :
```
FCOMP = -F gfortran --fargs "-fPIC" --cargs "-fPIC"
FCOMPEXEC =  --fcompexec ftn
INSTALLOPTIONS = --home=$(SCRATCH)/warp_install
```
  
- Then, in the directory `warp/pywarp90`, enter `make pinstall`. The compilation then lasts for a few minutes.

### In order to install warp with the intel compiler

- Download the source by entering `git clone https://bitbucket.org/berkeleylab/warp.git`

- Go to the `warp/pywarp90` directory and create a file called `Makefile.local.pympi`. Enter the following lines in this file :
```
FCOMP = -F intel --fargs "-fPIC" --cargs "-fPIC"
FCOMPEXEC =  --fcompexec ftn
INSTALLOPTIONS = --home=$(SCRATCH)/warp_install
```
  
- Be sure you are using the intel compiler: `module load PrgEnv-intel`  
  
- Then, in the directory `warp/pywarp90`, enter `make -f Makefile.Forthon.pympi install`. The compilation then lasts for a few minutes.

## Running simulations

In order to run a simulation, create a new directory,
copy your Warp input script to this directory, and rename this script
to `warp_script.py`. (The folder `scripts/examples/` of the
[Warp repository](https://bitbucket.org/berkeleylab/warp/src) contains
several examples of input scripts.)

Then create a submission script named `submission_script`. Here is an
example of a typical submission script.
```
#!/bin/bash
#SBATCH --job-name=test_simulation
#SBATCH --time=00:30:00
#SBATCH -n 32
#SBATCH --partition=debug
#SBATCH -e test_simulation.err
#SBATCH -o test_simulation.out

export mydir="$SCRATCH/test_simulation"
rm -fr $mydir
mkdir -p $mydir

cd $SLURM_SUBMIT_DIR

cp ./* $mydir/.
cd $mydir

srun -n 32 python-mpi -i warp_script.py -p 2 1 16
```

Then submit the simulation by typing `sbatch submission_script`.  The
progress of the simulation can be seen by typing ```squeue -u `whoami` ```. 
