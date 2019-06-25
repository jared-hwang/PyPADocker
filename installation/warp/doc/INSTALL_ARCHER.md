# Warp on ARCHER

**This setup document explains how to install Warp on
  [ARCHER](http://www.archer.ac.uk/), a supercomputer located in the UK.**

ARCHER provides many of the necessary dependencies as modules, and these should be loaded first:
```
module load python-compute pc-numpy pc-scipy cray-hdf5
```
Installation has only been tested using the gnu compiler set, for which the environment must be swapped:
```
module swap PrgEnv-cray PrgEnv-gnu
```

However, installing warp and some of its other dependencies into a python distribution requires the user to create a python environment. First of all, note that this *must* be done on the `/work` filesystem, otherwise it will not be accessible by the backend during runtime. When in a suitable location, create the environment by entering:
```
virtualenv venv
```
where `venv` is a name for your new python distribution. To use this, enter
```
source venv/bin/activate
```

## Dependencies
Packages not provided by the supplied python distribution are Forthon, h5py, and python-dateutil. After setting up the compiler wrappers, these can be installed using pip:
```
export CC=cc
export FC=ftn
pip install python-dateutil Forthon h5py
```
If there is an error regarding cython, the location should be added to your `PATH` and then try again. At the time of writing, this is done with
```
PATH="/fs2/y07/y07/cse/python/modules/cython/0.21.1/bin/:$PATH"
```

### pygist
To install pygist:
```
git clone https://bitbucket.org/dpgrote/pygist.git
cd pygist
python setup.py config
python setup.py install
cd ..
```

## Warp
We should now be ready to install **warp**.
```
git clone https://bitbucket.org/berkeleylab/warp.git
cd warp/pywarp90
```
The compiler wrappers should sort out all the linking, so all we need to do is tell Fothon which compiler to use:
```
echo "FORTHON = Forthon -F gfortran --fcompexec ftn" > Makefile.local.pympi
```
We also need to tell the system how to build the final module:
```
export LDSHARED="ftn -shared"
```
Finally, **warp** can be built and installed with:
```
make pinstall
```

## Running
Parallel simulations can be run using `aprun`:
```
aprun -n X python script.py
```
where `X` is the number of processes required and `script.py` is your **warp** script. At the time of writing, it is necessary to use a different version of `mpi4py` than default, in which case this should be performed *before* running:
```
PYTHONPATH=/work/y07/y07/cse/mpi4py/2.0.0/lib/python2.7/site-packages:$PYTHONPATH
```
