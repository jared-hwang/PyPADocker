# Warp on quartz.llnl.gov

**This setup document explains how to install Warp on
  [QUARTZ](quart.llnl.gov), a supercomputer located at LLNL.**

The best version of python to use is the one from the modules.
```
module load python
```

Installation has only been tested using the gnu compiler set, for which the environment must be swapped:
```
module swap intel gcc
```

## Installation location
For convenience, the user option is used when installing. The installation directory can be controlled by
setting the environment variable PYTHONUSERBASE. If unset, the installation directory defaults to ~/.local.
If setting PYTHONUSERBASE, it needs to be set both during installation and running, so it should be set
in the shell resource file.

## Dependencies
Forthon is not provided by the supplied python distribution. It must be downloaded (the pip version is broken on this machine).
The path of Python in the Forthon executable is wrong when using pip.
```
git clone https://github.com/dpgrote/Forthon.git
cd Forthon
python setup.py install --user
cd ..
```

## pygist
To install pygist:
```
git clone https://bitbucket.org/dpgrote/pygist.git
cd pygist
python setup.py config
python setup.py install --user
cd ..
```

## Path
The execute path needs to be setup to include the new bin directory where Forthon and the gist viewer are installed.

The following commands are with tcsh. Note that the third line adds the path modification to the csh resource file
so that you will have the modified path every login. If using PYTHONUSERBASE, "~/.local" should be replaced
with $PYTHONUSERBASE.
```
set path=(~/.local/bin $path)
rehash
echo 'set path = (~/.local/bin $path)' >> ~/.cshrc
```

## Warp
We should now be ready to install **warp**.
```
git clone https://bitbucket.org/berkeleylab/warp.git
cd warp/pywarp90
```
The compiler wrappers should sort out all the linking, so all we need to do is tell Fothon which compiler to use:
```
echo "FCOMP = -F gfortran" > Makefile.local
echo "INSTALLOPTIONS = --user" >> Makefile.local
echo "FCOMPEXEC =  --fcompexec mpif90" > Makefile.local.pympi
echo "" >> Makefile.local
echo "INSTALLOPTIONS = --user" >> Makefile.local.pympi
```

For the parallel version, the linking needs to be told about the parallel fortran libraries. Do the following command
to create the setup.local.py file which has that information.

```
echo "if parallel:\
    library_dirs += ['/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-4.9.3/lib']\
    libraries += ['mpich', 'mpifort']" > setup.local.py
```

Finally, **warp** can be built and installed with (build both serial and parallel versions):
```
make -j install
make -j pinstall
```

