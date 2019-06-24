# Installing Warp with mpi4py on Cori

This document describes how to install Warp on the Cori cluster.
The installation can be done by compiling the sources (recommended).
It can also be done by using Shifter (see the corresponding section).

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
if [ "$NERSC_HOST" == "cori" ]
then
   module swap PrgEnv-intel PrgEnv-gnu
   module load python/2.7-anaconda
   module load h5py-parallel
   WARP=$SCRATCH/warp_install/
   export PATH=$WARP/bin:$PATH
   export PYTHONPATH=$WARP/lib/python:$PYTHONPATH
fi
```

Then run `source .bashrc` in the `$HOME`directory (or alternatively
log out of Cori, and then log in again).

## Installing Forthon

As explained on the Warp website, in order to install Forthon :

- Type `git clone https://github.com/dpgrote/Forthon.git`
- Then `cd` into the directory `Forthon` and run `python setup.py install --home=$WARP`

## Installing Warp itself

### In order to install warp with gfortran

- Download the source by entering `git clone https://bitbucket.org/berkeleylab/warp.git`

- Go to the `warp/pywarp90` directory and create a file called `Makefile.local.pympi`. Enter the following lines in this file :
```
FCOMP = -F gfortran --fargs "-fPIC" --cargs "-fPIC"
FCOMPEXEC =  --fcompexec ftn
INSTALLOPTIONS = --home=$(SCRATCH)/warp_install/
```

- Then, in the directory `warp/pywarp90`, enter `make pinstall`. The compilation then lasts for a few minutes.

### In order to install Warp with the Intel compiler

- Download the source by entering `git clone https://bitbucket.org/berkeleylab/warp.git`

- Go to the `warp/pywarp90` directory and create a file called `Makefile.local.pympi`.
If you want to compile for Haswell Architecture, enter the following lines in this file :

```
FCOMP = -F intel --fargs "-fPIC -O3 -xCORE-AVX2" --cargs "-fPIC"
FCOMPEXEC =  --fcompexec ftn
INSTALLOPTIONS = --home=$(SCRATCH)/warp_install/
```

For MIC architecture, such as Intel Xeon Phi KNL, replace `-xCORE-AVX2` by `-xMIC-AVX512`.
This activates the use of AVX512 instructions for best performance on KNL.
However, note that KNL supports previous Intel instructions and your code will
work even compiled for Haswell (Cori phase 1) or Ivy Bridge (Edison) architectures.

```
FCOMP = -F intel --fargs "-fPIC -O3 -xMIC-AVX512" --cargs "-fPIC"
FCOMPEXEC =  --fcompexec ftn
INSTALLOPTIONS = --home=$(SCRATCH)/warp_install/
```

- Be sure you are using the intel compiler: `module load PrgEnv-intel`  

- Then, in the directory `warp/pywarp90`, enter `make -f Makefile.Forthon.pympi install`.
The compilation then lasts for a few minutes.

### Running simulations

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
#SBATCH -C haswell
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

## Using Shifter

You can also use Shifter
[Shifter](http://www.nersc.gov/research-and-development/user-defined-images/)
to run Warp on Cori. Shifter handles Linux container (similar to
Docker), and allows to easily port codes from one
architecture to another.

Shifter should **not** use the installation of Warp which is on the
system itself. Therefore, if you previously compiled and installed
Warp in `$SCRATCH/warp_install`, please remove this directory to avoid
errors.
```
rm -rf $SCRATCH/warp_install
```

In addition, make sure that your `$PATH` and `$PYTHONPATH` variables
are unmodified, and that no module is loaded.
In particular, if you used to modify these variables
in `.bashrc.ext`, please remove or comment out the corresponding
lines, as in the following example (note that the lines are commented out):

```
# if [ "$NERSC_HOST" == "cori" ]
# then
#  module swap PrgEnv-intel PrgEnv-gnu
#	module load python/2.7-anaconda
#	module load mpi4py
#  export PATH=$SCRATCH/warp_install/bin:$PATH
#	export PYTHONPATH=$SCRATCH/warp_install/lib/python:$PYTHONPATH
# fi
```


### Running simulations with Shifter

In order to run a simulation, create a new directory,
copy your Warp input script to this directory, and rename this script
to `warp_script.py`. (The folder `scripts/examples/` of the
[Warp repository](https://bitbucket.org/berkeleylab/warp/src) contains
several examples of input scripts.)

Then create a submission script named `submission_script`. Here is an
example of a typical submission script.
```
#!/bin/csh
#SBATCH --job-name=test_cori_shifter
#SBATCH --time=00:10:00
#SBATCH -N 1
#SBATCH --partition=debug
#SBATCH -e test_cori_shifter.err
#SBATCH -o test_cori_shifter.out
#SBATCH -C haswell
#SBATCH --image=docker:rlehe/warp:latest
#SBATCH --volume=<your$SCRATCH>:/home/warp_user/run

setenv mydir "$SCRATCH/test_simulation"
rm -fr $mydir
mkdir -p $mydir

cd $SLURM_SUBMIT_DIR

cp ./* $mydir/.
cd $mydir

setenv OMP_NUM_THREADS 1
srun -n 32 -c 2 shifter python warp_script.py -p 4 1 8
```
Note that the options `--image=docker:rlehe/warp:latest` and `--
volume=<your$SCRATCH>:/home/warp_user/run` are essential
and should be copied exactly (**do not** replace `warp_user` or
`rlehe` by your username), with the exception of `<your$SCRATCH>`,
which should be replaced by the full path to your SCRATCH directory.

Then submit the simulation by typing `sbatch submission_script`.  The
progress of the simulation can be seen by typing ```squeue -u `whoami` ```.
