#!/usr/bin/env python
# To use:
#       python setup.py install
#
import sys
import os
import subprocess
from Forthon.compilers import FCompiler
import getopt

try:
    import distutils
    from distutils.core import setup, Extension
    from distutils.dist import Distribution
    from distutils.command.build import build
except:
    raise SystemExit('Distutils problem')

optlist, args = getopt.getopt(sys.argv[1:], 'gt:F:',
                              ['parallel', 'fargs=', 'cargs=', 'fcompexec='])

machine = sys.platform
debug   = 0
fcomp   = None
parallel = 0
fcompexec = None
for o in optlist:
    if   o[0] == '-g': debug = 1
    elif o[0] == '-t': machine = o[1]
    elif o[0] == '-F': fcomp = o[1]
    elif o[0] == '--parallel': parallel = 1
    elif o[0] == '--fcompexec': fcompexec = o[1]

sys.argv = ['setup.py'] + args
fcompiler = FCompiler(machine=machine,
                      debug=debug,
                      fcompname=fcomp,
                      fcompexec=fcompexec)

dummydist = Distribution()
dummydist.parse_command_line()
dummybuild = dummydist.get_command_obj('build')
dummybuild.finalize_options()
builddir = dummybuild.build_temp

if dummydist.commands[-1] == 'install':
    # --- During an install, remove the build/lib directory, since distutils
    # --- doesn't update an older warpC.so even if there were changes.
    os.system('rm -rf %s'%dummybuild.build_platlib)

warppkgs = ['top', 'env', 'w3d', 'f3d', 'wxy', 'fxy', 'wrz', 'frz', 'her', 'cir', 'cho', 'em3d']

def makeobjects(pkg):
    return [pkg + '.o', pkg + '_p.o', pkg + 'pymodule.o']

warpobjects = []
if sys.hexversion < 0x03000000:
    # --- With Python2, everything is put into the one warpC.so file.
    for pkg in warppkgs:
        warpobjects = warpobjects + makeobjects(pkg)

    warpobjects = warpobjects + ['top_lattice.o', 'top_fsl.o', 'dtop.o',
                                 'dw3d.o', 'w3d_injection.o', 'w3d_interp.o',
                                 'w3d_collisions.o', 'w3d_utilities.o', 'w3d_load.o',
                                 'f3d_mgrid.o', 'f3d_ImplicitES.o', 'f3d_mgrid_be.o',
                                 'f3d_conductors.o', 'f3d_bfield.o',
                                 'fft.o', 'util.o',
                                 'fxy_mgrid.o',
                                 'dwrz.o',
                                 'frz_mgrid.o', 'frz_mgrid_be.o', 'frz_ImplicitES.o',
                                 #'em2d_apml.o', 'em2d_apml_cummer.o', 'em2d_maxwell.o',
                                 'em3d_maxwell.o']
    if parallel:
        warpobjects = warpobjects + ['f3dslave.o', 'frzslave.o', 'topslave.o',
                                     'w3dslave.o']

    warpobjects = map(lambda p:os.path.join(builddir, p), warpobjects)

    # --- With Python3, each packages has it's own .so file and warpC.so is mostly a dummy package.

library_dirs = fcompiler.libdirs
libraries = fcompiler.libs
extra_link_args = ['-g'] + fcompiler.extra_link_args
include_dirs = [builddir]

# --- Add the numpy include paths
import numpy
define_macros = []
include_dirs.append(numpy.get_include())

if parallel:
    # --- This is only needed by warpC_Forthon.c
    define_macros += [('MPIPARALLEL', None)]

datestring = os.popen('git log --branches=master --remotes=origin -n 1 --pretty=%aD').read().strip()
define_macros += [('GITORIGINDATE', '"' + datestring + '"')]
datestring = os.popen('git log -n 1 --pretty=%aD').read().strip()
define_macros += [('GITLOCALDATE', '"' + datestring + '"')]
commithash = os.popen('git log -n 1 --pretty=%h').read().strip()
define_macros += [('GITCOMMITHASH', '"' + commithash + '"')]

define_macros.append(('FORTHON_PKGNAME', '"warp"'))

if parallel:
    name = 'warpCparallel'
else:
    name = 'warpC'

# --- The behavior of distutils changed from 2.2 to 2.3. In 2.3, the object
# --- files are always put in a build/temp directory relative to where the
# --- source file is, rather than relative to the main build directory.
# --- This tells distutils to put the objects in the same directory
# --- as the source files.
if sys.hexversion >= 0x020300f0 and dummydist.commands[-1] == 'build':
    sys.argv += ['--build-temp', '']

if machine == 'darwin':
# --- Machines running csh/tcsh seem to have MACHTYPE defined and this is the safest way to set -arch.
    try:
        machtype = os.environ['MACHTYPE']
        if machtype == 'i386':
            os.environ['ARCHFLAGS'] = '-arch i386'
        elif machtype == 'x86_64':
            os.environ['ARCHFLAGS'] = '-arch x86_64'
        elif machtype == 'powerpc':
            os.environ['ARCHFLAGS'] = '-arch ppc'
    except KeyError:
# ---  If the shell is bash, MACHTYPE is undefined.  So get what we can from uname. We will assume that if
# ---  we are running Snow Leopard we are -arch x86-64 and if running Leopard on intel we are -arch i386.
# ---  This can be over-ridden by defining MACHTYPE.
        archtype = os.uname()[-1]
        if archtype in ['Power Macintosh', 'ppc']:
            os.environ['ARCHFLAGS'] = '-arch ppc'
        elif archtype in ['i386', 'x86_64']:
            kernel_major = eval(os.uname()[2].split('.')[0])
            if kernel_major < 10 :
                os.environ['ARCHFLAGS'] = '-arch i386'  # Leopard or earlier
            else:
                os.environ['ARCHFLAGS'] = '-arch x86_64'  # Snow Leopard

# --- Check if there is a file setup.local.py holding local definitions
# --- that might be needed to build Warp. Note that execfile is used so
# --- that everything defined up to this point is available, and anything
# --- can be redefined.
if os.access('setup.local.py', os.F_OK):
    if sys.hexversion < 0x03000000:
        execfile('setup.local.py')
    else:
        exec(compile(open('setup.local.py').read(), 'setup.local.py', 'exec'))

elif parallel:
    if fcompexec is None or fcompexec == 'mpifort' or fcompexec == 'mpif90':
        # --- If parallel, try the "mpifort -show" to get the mpi libraries
        try:
            show = subprocess.check_output([fcompexec, '-show'], universal_newlines=True, stderr=subprocess.STDOUT)
            for s in show.split():
                if s.startswith('-L'):
                    library_dirs += [s[2:]]
                elif s.startswith('-l'):
                    libraries += [s[2:]]
        except (OSError, subprocess.CalledProcessError):
            pass

setup (name = 'warp',
       version = '4.5',
       author = 'David P. Grote, Jean-Luc Vay, et. al.',
       author_email = 'dpgrote@lbl.gov',
       description = 'Warp PIC accelerator code',
       long_description = """
Warp is a PIC code designed to model particle accelerators and similar
machines that are space-charge dominated.""",
       url = 'http://warp.lbl.gov',
       platforms = 'Linux, Unix, Windows (bash), Mac OSX',
       ext_modules = [Extension('warp.' + name,
                                ['warpC_Forthon.c',
                                 os.path.join(builddir, 'Forthon.c'),
                                 'pmath_rng.c', 'ranf.c', 'ranffortran.c'],
                                include_dirs=include_dirs,
                                library_dirs=library_dirs,
                                libraries=libraries,
                                define_macros=define_macros,
                                extra_objects=warpobjects,
                                extra_link_args=extra_link_args,
                                extra_compile_args=fcompiler.extra_compile_args
                               )]

       )
