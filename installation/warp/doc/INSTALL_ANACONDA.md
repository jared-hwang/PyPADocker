# Installing Warp using Anaconda

Here are the steps needed to install Warp using Anaconda version 5 or higher. These are the versions that can include a full compiler toolset. It is recommended to use a version with Python3.

## Setting up Anaconda

This assumes that Anaconda has already been installed.
The appropriate C and Fortran compiler needs to be installed. Details can be found at the Anaconda website, https://www.anaconda.com/utilizing-the-new-compilers-in-anaconda-distribution-5/.

Install the compilers for your architecture.

Linux:
```
conda install gcc_linux-64
conda install gfortran_linux-64
```

MacOS:
```
conda install clang_osx-64
conda install gfortran_osx-64
```

In order to use the compilers, they must be activated, for example:
```
source activate root
```

This sets up the environment in the shell, for example setting the appropriate execute path. The compilation must be done in this activated shell.

## Installing Forthon

Forthon provides the glue linking Fortran and python. It must be installed first before Warp can be built. It can be installed using pip:

```
pip install Forthon
```

## Installing Serial Warp

Download Warp from bitbucket, https://bitbucket.org/berkeleylab/warp/downloads. Download and extract the tar ball from the most recent release (not the "Download repository" link). Update it by cd-ing into the warp directory and running git.
```
git pull
```

Warp needs to be configured to use the Anaconda compilers. Go into the warp/pywarp90 directory and create the file Makefile.local3.

```
echo 'FCOMP = -F gfortran' >> Makefile.local3
echo 'FCOMPEXEC = --fcompexec x86_64-apple-darwin13.4.0-gfortran' >> Makefile.local3
```

Replace the fcompexec argument with the appropriate name of the gfortran compiler for your system. (Don't forget that if you upgrade your system, this name will need to be updated also.)

Serial Warp can now be compiled and installed.
```
make -j install3
```

## Installing Parallel Warp

To use parallel Warp, the mpi4py package needs to be installed. I recommend using the conda forge version (which has more versions available).
```
conda config --add channels conda-forge
conda install mpi4py
```

As with the serial version, the appropriate compiler must be specified. Create the file Makefile.local3.pympi.
```
echo 'FCOMP = -F gfortran' >> Makefile.local3.pympi
echo 'FCOMPEXEC = --fcompexec mpifort' >> Makefile.local3.pympi
```

Then build and install it.
```
make -j pinstall3
```

## Installing pygist

To install pygist, the plotting package:
```
git clone https://bitbucket.org/dpgrote/pygist.git
cd pygist
python setup.py config
python setup.py install
cd ..
```

