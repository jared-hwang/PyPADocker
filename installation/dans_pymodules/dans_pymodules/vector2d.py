import numpy as np

__author__ = "Daniel Winklehner"
__doc__ = "Simple 2D vector class covering basic vector operations."


class Vector2D(object):
    """
    Simple 2D vector class covering basic vector operations
    """

    def __init__(self, vector=np.ones(2, 'd'), p0=None, p1=None):
        """
        Constructor
        """
        if p0 is not None and p1 is not None:
            assert len(p0) == 2 and len(p1) == 2, "Both points have to have dimension 2 for 2D Vector!"
            self.vector = np.array(p1, 'd') - np.array(p0, 'd')
        else:
            assert len(vector) == 2, "Length has to be 2 for 2D Vector!"
            self.vector = np.array(vector, 'd')

        self.length = self.get_length()

    def __getitem__(self, item):
        return self.vector[item]

    def __str__(self):
        return str(self.vector)

    def __mul__(self, other):
        try:

            return self[0] * other.vector[0] + self[1] * other.vector[1]

        except AttributeError:

            return Vector2D(vector=(self.vector * other))

    def __rmul__(self, other):
        try:

            return self[0] * other.vector[0] + self[1] * other.vector[1]

        except AttributeError:

            return Vector2D(vector=(self.vector * other))

    def __div__(self, other):
        return Vector2D(vector=(self.vector / other))

    def __rdiv__(self, other):
        return Vector2D(vector=(other / self.vector))

    def __add__(self, other):
        try:

            return Vector2D(vector=(self.vector + other.vector))

        except AttributeError:

            print("Can only add and subtract vectors from vectors!")

    def __radd__(self, other):
        try:

            return Vector2D(vector=(self.vector + other.vector))

        except AttributeError:

            print("Can only add and subtract vectors from vectors!")

    def __sub__(self, other):
        try:

            return Vector2D(vector=(self.vector - other.vector))

        except AttributeError:

            print("Can only add and subtract vectors from vectors!")

    def __rsub__(self, other):
        try:

            return Vector2D(vector=(other.vector - self.vector))

        except AttributeError:

            print("Can only add and subtract vectors from vectors!")

    def get_length(self):
        return np.sqrt(self.vector[0] ** 2.0 + self.vector[1] ** 2.0)

    def normalize(self):

        return Vector2D(vector=(self.vector / self.length))

    def rotate_ccw(self):

        return Vector2D(vector=np.array([-self.vector[1], self.vector[0]]))

    def rotate_cw(self):

        return Vector2D(vector=np.array([self.vector[1], -self.vector[0]]))

    def angle(self, other):
        """
        Calculate the angle from the present vector to the given second_vector
        """
        return np.arccos(self * other / self.get_length() / other.get_length())


if __name__ == '__main__':
    # Test the vector class
    v1 = Vector2D(np.array([1, 2]))
    v2 = Vector2D(np.array([-1, 2]))

    print(v1 + v2)
    print(v1 - v2)
    print(v1 * v2)
    print(v1.get_length())
