import numpy as np

__author__ = "Daniel Winklehner"
__doc__ = "Simple matrix based coordinate transformation in 3D"


class TransformationMatrix(object):

    def __init__(self, matrix=None):
        """
        A 3x3 transformation matrix that can be multiplied with lists or numpy arrays of shape (3,) or (3, N)
        :param matrix: a 3x3 list or numpy array. If None, identity transformation will be initialized.
        """

        if matrix is not None:
            self.matrix = matrix
        else:
            self.matrix = np.eye(3, 3)

        self.epsilon = 1.0e-10  # Small value to test for identities

    @property
    def epsilon(self):
        return self._epsilon

    @epsilon.setter
    def epsilon(self, epsilon):
        self._epsilon = epsilon

    def matrix_from_vectors(self, _u, _v, _w, normalized=False):
        """
        Creates a matrix from three orthogonal vectors. These could be 
        e.g. the local coordinate vectors (or non-normalized directions) in the global frame.
        :param _u: local 
        :param _v: 
        :param _w: 
        :param normalized: 
        :return: 
        """

        assert np.dot(_u, _v) <= self.epsilon and np.dot(_v, _w) <= self.epsilon, "Input vectors are not orthogonal!"

        if not normalized:
            _u /= np.linalg.norm(_u)
            _v /= np.linalg.norm(_v)
            _w /= np.linalg.norm(_w)

        self.matrix = np.array([_u,
                                _v,
                                _w]).T

        pass

    def __mul__(self, other):

        other = np.array(other)

        assert other.ndim <= 2, "Error: Cannot multiply 3x3 matrix with object of dimensions {}".format(other.ndim)

        if other.ndim == 1:

            assert other.shape[0] == 3, "Error: Point has to have three entries, got {}".format(other.shape[0])

            return np.dot(self.matrix, other)

        elif other.ndim == 2:

            ncols, nrows = other.shape

            assert nrows == 3, "Error: Points have to have three entries, got {}".format(nrows)

            return np.array([np.dot(self.matrix, point) for point in other])


if __name__ == '__main__':

    mat = TransformationMatrix()

    u = [1, 1, 0]
    v = [-1, 1, 0]
    w = [0, 0, 1]

    mat.matrix_from_vectors(u, v, w)

    v2 = [1, 2, 3]

    # noinspection PyTypeChecker
    print(mat * v2)
