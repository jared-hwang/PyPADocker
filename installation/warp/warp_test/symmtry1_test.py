from warp import *
import unittest

# --- Set four-character run id, comment lines, user's name.
top.pline2   = "Test of 3d field symmetries"
top.pline1   = "test 1"
top.runmaker = "Arun Persaud"

# --- Turn off plot output
top.lprntpara = false
top.lpsplots = false

# --- Setup the grid size
w3d.nx = 4
w3d.ny = 4
w3d.nz = 4

w3d.xmmin = -1.
w3d.xmmax =  1.
w3d.ymmin = -1.
w3d.ymmax =  1.
w3d.zmmin = -1.
w3d.zmmax =  1.

ions = Species(type=Hydrogen, charge_state=1, name='H')
ions.sw = 1.e8

# --- for field solve
w3d.bound0 = dirichlet
w3d.boundnz = dirichlet
w3d.boundxy = dirichlet

w3d.l4symtry = false
w3d.l2symtry = true
if w3d.l4symtry:
    w3d.xmmin = 0.
    w3d.nx /= 2
if w3d.l2symtry or w3d.l4symtry:
    w3d.ymmin = 0.
    w3d.ny /= 2

@callfromparticleloader
def loadparticles():
    ions.addparticles(x=[-0.1,0.1], y=0., z=0.)

package('w3d')
generate()

class TestSymmetryFunctions(unittest.TestCase):
    def test_xsymmetry1(self):
        self.assertAlmostEqual(-ions.ex[0], ions.ex[1], places=12)

    def test_xsymmetry2(self):
        self.assertLess(ions.ex[0], 0.)

    def test_xsymmetry3(self):
        self.assertGreater(ions.ex[1], 0.)

    def test_ysymmetry1(self):
        self.assertAlmostEqual(ions.ey[0], 0., places=12)

    def test_ysymmetry2(self):
        self.assertAlmostEqual(ions.ey[1], 0., places=12)

    def test_zsymmetry1(self):
        self.assertAlmostEqual(ions.ez[0], 0., places=12)

    def test_zsymmetry2(self):
        self.assertAlmostEqual(ions.ez[1], 0., places=12)

if __name__ == '__main__':
    unittest.main()

