# Installing Warp on Mac OSX

## Installing MPI

MPI can be installed **either** with MacPorts or HomeBrew.

### With MacPorts

If you do not have MacPorts installed, install it from [this page](https://www.macports.org/install.php). Then run:

```
sudo port selfupdate
sudo port install gcc48 +gfortran
sudo port install openmpi-gcc48
```

If MacPort suggests to make `openmpi-gcc48` the default mpi, please follow this advice (typically by doing ```sudo port select --set mpi openmpi-gcc48-fortran```).

Open a new terminal and check that the default `mpif90` is the MacPorts version by running `mpif90 -v`. (It should return some text finishing by a line of the form `gcc version 4.8.5 (MacPorts gcc48 4.8.5_0)`). If you get an error of the form `The Open MPI wrapper compiler was unable to find the specified compiler gfortran in your PATH.`, this can be fixed by running `sudo ln -s /opt/local/bin/gfortran-mp-4.8 /opt/local/bin/gfortran`.

Check also that the default `gcc` is the MacPorts version by running
`gcc -v` (It should return some text finishing by a line of the form
`gcc version 4.8.5 (MacPorts gcc48 4.8.5_0)`). If it is not the case,
run `sudo port select gcc mp-gcc48`

### With Homebrew

If you do not have HomeBrew installed, install it from [this page](http://brew.sh/). Then run:

```brew install gcc``` (This process is very
long, the installation hangs at `make bootstrap` and takes 30 min to
1 hour to complete. Please be patient.)

```brew install openmpi```

Open a new terminal and check that the default `mpif90` is the Homebrew version by running
`mpif90 -v`. (This should return some text finishing by a line of the form `gcc version 4.8.4 (Homebrew gcc48 4.8.4)`.)

## Installing Anaconda

If you do not have Anaconda installed, install it from [this page](https://www.continuum.io/downloads#_macosx) (Choose the Python2.7 version).

Open a new terminal and make sure that the default `python` is the Anaconda version. (Running `python` should print a message which starts with a line of the form `Python 2.7.10 |Anaconda 2.3.0`)

## Installing mpi4py

Make sure that the Anaconda version of `mpi4py` is not installed (it is buggy):

```conda uninstall mpi4py```

Install `mpi4py` with pip:

`pip install mpi4py`

## Installing Forthon

Run `pip install Forthon`

## Installing Warp

Run ```git clone https://bitbucket.org/berkeleylab/warp.git```

then `cd` into the repository `warp/pywarp90` and create two files:

- A file named `Makefile.local.pympi` which contains the following text:

```FCOMP= -F gfortran```

- A file named `setup.local.py` which contains the following text:

    If you used MacPorts:
```python
if parallel:
	library_dirs += ['/opt/local/lib/openmpi-gcc48']
	libraries = fcompiler.libs + ['mpi_usempi', 'mpi_mpifh', 'mpi']
```

    If you used HomeBrew:
```python
if parallel:
    library_dirs += ['/usr/local/Cellar/open-mpi/1.10.1_1/lib']
    libraries = fcompiler.libs + ['mpi_usempif08', 'mpi_usempi_ignore_tkr', 'mpi_mpifh', 'mpi']
```

Then install Warp by running:
```
make install
make pinstall
```

## Running a simulation

In order to run a simulation, create a new directory,
copy your Warp input script to this directory, and rename this script
to `warp_script.py`. (The folder `scripts/examples/` of the
[Warp repository](https://bitbucket.org/berkeleylab/warp/src) contains
several examples of input scripts.)

Then launch the simulation by typing either (for serial simulations)
```
python warp_script.py
```
or (for e.g. a parallel simulation with a 2x3x2 domain decomposition in 3D)
```
mpirun -np 12 python warp_script.py -p 2 3 2
```

