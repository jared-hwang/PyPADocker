#!/usr/bin/env python
# To use:
#       python setup.py install
#

import os
import sys
import glob

try:
    import distutils
    from distutils.command.install import INSTALL_SCHEMES
    # --- For installing into site-packages, with python setup.py install
    from distutils.core import setup
    # --- For creating an egg file, with python setup.py bdist_egg
    #from setuptools import setup
except:
    raise SystemExit('Distutils problem')

try:
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
    from distutils.command.build_py import build_py

# --- Get around a "bug" in disutils on 64 bit systems. When there is no
# --- extension to be installed, distutils will put the scripts in
# --- /usr/lib/... instead of /usr/lib64. This fix will force the scripts
# --- to be installed in the same place as the .so (if an install was
# --- done in pywarp90 - though that shouldn't be done anymore).
if distutils.sysconfig.get_config_vars()["LIBDEST"].find('lib64') != -1:
    for scheme in INSTALL_SCHEMES.values():
        scheme['purelib'] = scheme['platlib']

# --- Write out git versioning information
with open('__version__.py', 'w') as ff:
    ff.write('__origindate__ = "%s"\n'%os.popen('git log --branches=master --remotes=origin -n 1 --pretty=%aD').read().strip())
    ff.write('__localdate__ = "%s"\n'%os.popen('git log -n 1 --pretty=%aD').read().strip())
    ff.write('__hash__ = "%s"\n'%os.popen('git log -n 1 --pretty=%h').read().strip())
    ff.write('__fullhash__ = "%s"\n'%os.popen('git log -n 1 --pretty=%H').read().strip())

setup(name='warp',
      version='4.5',
      author='David P. Grote, Jean-Luc Vay, et. al.',
      author_email='dpgrote@lbl.gov',
      description='Warp PIC accelerator code',
      long_description="""
Warp is a PIC code designed to model particle accelerators and similar
machines that are space-charge dominated.""",
      url='http://warp.lbl.gov',
      platforms='Linux, Unix, Windows (cygwin), Mac OSX',
      packages=['warp', 'warpoptions', 'warp_parallel',
                'warp.attic',
                'warp.data_dumping',
                'warp.data_dumping.openpmd_diag',
                'warp.diagnostics',
                'warp.diagnostics.palettes',
                'warp.envelope',
                'warp.field_solvers',
                'warp.field_solvers.laser',
                'warp.GUI',
                'warp.init_tools',
                'warp.lattice',
                'warp.particles',
                'warp.run_modes',
                'warp.utils'],
      package_dir={'warp': '.'},
      package_data={'warp': ['diagnostics/palettes/*.gs',
                             'diagnostics/palettes/*.gp',
                             'particles/aladdin_8.txt']},
      cmdclass={'build_py': build_py}
      )
