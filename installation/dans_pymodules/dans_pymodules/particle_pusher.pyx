import sys
import numpy as np
cimport numpy as np
cimport cython

if sys.version_info.major == 3:
    from .particles import IonSpecies
    from .field import Field
else:
    from particles import IonSpecies
    from field import Field

DTYPE1 = np.float64
DTYPE2 = np.int
ctypedef np.float64_t DTYPE1_t
ctypedef np.int_t DTYPE2_t

class ParticlePusher(object):

    def __init__(self, ion, algorithm="boris"):

        self._alg_switch = {
            "boris": self._v_boris,
            "leapfrog": self._v_leapfrog,
            "tajima_implicit": self._v_tajima_implicit
        }

        assert isinstance(ion, IonSpecies), "'ion' needs to be an IonSpecies object!"

        self._ion = ion

        algorithm = algorithm.lower()

        assert algorithm in self._alg_switch.keys(), "'{}' is not a known algorithm!".format(algorithm)

        self._algorithm = algorithm

        self._v = self._alg_switch[self._algorithm]

        self._efield = None
        self._bfield = None
        self._bounding_electrodes = None

    def set_efield(self, efield):
        self._efield = efield

    def set_bfield(self, bfield):
        self._bfield = bfield

    def set_bds(self, electrodes):
        self._bounding_electrodes = electrodes

    @cython.boundscheck(False)  # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def track(self,
              np.ndarray[DTYPE1_t, ndim=1] r0,
              np.ndarray[DTYPE1_t, ndim=1] v0,
              DTYPE2_t nsteps,
              DTYPE1_t dt):

        assert self._efield is not None or self._efield is not None, "Either efield or bfield (or both) have to be set!"

        cdef np.ndarray[DTYPE1_t, ndim=2] r = np.zeros([nsteps + 1, 3])
        cdef np.ndarray[DTYPE1_t, ndim=2] v = np.zeros([nsteps + 1, 3])

        r[0] = r0
        v[0] = v0

        # initialize the velocity half a step back:
        cdef np.ndarray[DTYPE1_t, ndim=1] ef = self._efield(r[0])
        cdef np.ndarray[DTYPE1_t, ndim=1] bf = self._bfield(r[0])

        _, v[0] = self.push(r[0], v[0], ef, bf, -0.5 * dt)

        cdef int i = 0

        if self._bounding_electrodes is None:
            # Track without electrode boundary check
            for i in range(nsteps):
                ef = self._efield(r[i])
                bf = self._bfield(r[i])

                r[i + 1], v[i + 1] = self.push(r[i], v[i], ef, bf, dt)

            # Push velocity one more half-step forward
            # TODO: Think about where to grab the fields and which r to use...
            ef = self._efield(r[nsteps])
            bf = self._bfield(r[nsteps])
            _, v[nsteps] = self.push(r[nsteps], v[nsteps], ef, bf, 0.5 * dt)

        else:
            # Track with electrode boundary check
            for i in range(nsteps):
                ef = self._efield(r[i])
                bf = self._bfield(r[i])

                r[i + 1], v[i + 1] = self.push(r[i], v[i], ef, bf, dt)

                if self._bounding_electrodes.points_inside(r[i + 1]):

                    r = r[:i + 1]
                    v = v[:i + 1]

                    break

        return r, v

    def algorithm(self, algorithm=None):

        if algorithm in self._alg_switch.keys():
            self._algorithm = algorithm
            self._v = self._alg_switch[self._algorithm]
        else:
            print("'{}' is not a known algorithm, keeping current algorithm ('{}')".format(algorithm, self._algorithm))

        return self._algorithm

    @cython.boundscheck(False)  # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def _v_boris(self,
                 np.ndarray[DTYPE1_t, ndim=1] _v,
                 np.ndarray[DTYPE1_t, ndim=1] _efield,
                 np.ndarray[DTYPE1_t, ndim=1] _bfield,
                 DTYPE1_t _dt):

        cdef np.ndarray[DTYPE1_t, ndim=1] t = 0.5 * self._ion.q_over_m() * _bfield * _dt
        cdef np.ndarray[DTYPE1_t, ndim=1] s = 2.0 * t / (1.0 + np.linalg.norm(t) ** 2.0)
        cdef np.ndarray[DTYPE1_t, ndim=1] v_minus = _v + 0.5 * self._ion.q_over_m() * _efield * _dt
        cdef np.ndarray[DTYPE1_t, ndim=1] v_prime = v_minus + np.cross(v_minus, t)
        cdef np.ndarray[DTYPE1_t, ndim=1] v_plus = v_minus + np.cross(v_prime, s)
        cdef np.ndarray[DTYPE1_t, ndim=1] result = v_plus + 0.5 * self._ion.q_over_m() * _efield * _dt

        return result

    @cython.boundscheck(False)  # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def _v_leapfrog(self,
                    np.ndarray[DTYPE1_t, ndim=1] _v,
                    np.ndarray[DTYPE1_t, ndim=1] _efield,
                    np.ndarray[DTYPE1_t, ndim=1] _bfield,
                    DTYPE1_t _dt):

        cdef np.ndarray[DTYPE1_t, ndim=1] result = _v + self._ion.q_over_m() * (_efield + np.cross(_v, _bfield)) * _dt

        return result

    @cython.boundscheck(False)  # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def _v_tajima_implicit(self,
                           np.ndarray[DTYPE1_t, ndim=1] _v,
                           np.ndarray[DTYPE1_t, ndim=1] _efield,
                           np.ndarray[DTYPE1_t, ndim=1] _bfield,
                           DTYPE1_t _dt):

        cdef np.ndarray[DTYPE1_t, ndim=1] b_norm = np.linalg.norm(_bfield)

        cdef np.ndarray[DTYPE1_t, ndim=2] rot_mat = (1.0 / b_norm) * np.array([[0.0, _bfield[2], -_bfield[1]],
                                                                               [-_bfield[2], 0.0, _bfield[0]],
                                                                               [-_bfield[1], -_bfield[0], 0.0]])

        cdef np.ndarray[DTYPE1_t, ndim=1] epsilon = 0.5 * self._ion.q_over_m() * b_norm * _dt
        cdef np.ndarray[DTYPE1_t, ndim=2] mat_p = np.eye(3) + epsilon * rot_mat
        cdef np.ndarray[DTYPE1_t, ndim=2] mat_m_inv = np.linalg.inv(np.eye(3) - epsilon * rot_mat)
        _v = np.matmul(np.matmul(mat_m_inv, mat_p), _v) + \
                                               np.matmul(mat_m_inv, _efield) * self._ion.q_over_m() * _dt

        return _v

    @cython.boundscheck(False)  # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def push(self,
             np.ndarray[DTYPE1_t, ndim=1] _r,
             np.ndarray[DTYPE1_t, ndim=1] _v,
             np.ndarray[DTYPE1_t, ndim=1] _efield,
             np.ndarray[DTYPE1_t, ndim=1] _bfield,
             DTYPE1_t _dt):

        _v = self._v(_v, _efield, _bfield, _dt)  # Call the velocity function determined by the algorithm

        return  _r + _v * _dt, _v

if __name__ == "__main__":
    h2p = IonSpecies(name="H2_1+", energy_mev=1.0, label="$\mathrm{H}_2^+$")
    h2p.calculate_from_energy_mev(0.07 / h2p.a())
    print("Cyclotron radius should be {} m".format(h2p.b_rho()))

    pusher = ParticlePusher(h2p, "boris")
    efield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": 0.0})
    bfield1 = Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.0})

    nsteps = 10000
    dt = 1e-10  # (1 ps --> s)

    r = np.zeros([nsteps + 1, 3])
    v = np.zeros([nsteps + 1, 3])
    r[0] = [h2p.b_rho(), 0.0, 0.0]
    v[0] = [0.0, h2p.v_m_per_s(), 0.0]

    # initialize the velocity half a step back:
    ef = efield1(r[0])
    bf = bfield1(r[0])
    _, v[0] = pusher.push(r[0], v[0], ef, bf, -0.5 * dt)

    for i in range(nsteps):

        ef = efield1(r[i])
        bf = bfield1(r[i])

        r[i + 1], v[i + 1] = pusher.push(r[i], v[i], ef, bf, dt)

    import matplotlib.pyplot as plt

    plt.plot(r[:, 0], r[:, 1])
    plt.gca().set_aspect('equal')
    plt.show()
