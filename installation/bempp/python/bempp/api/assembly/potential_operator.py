"""Definition of potential operators."""


class PotentialOperator:
    """Provides an interface to potential operators.

    This class is not supposed to be instantiated directly.
    """

    def __init__(self, op, component_count, space, evaluation_points):
        """Constructor. Should not be called by the user."""

        self._op = op
        self._component_count = component_count
        self._space = space
        self._evaluation_points = evaluation_points

    def evaluate(self, grid_fun):
        """
        Apply the potential operator to a grid function.

        Parameters
        ----------
        grid_fun : bempp.api.GridFunction
            A GridFunction object that represents the boundary density to
            which the potential is applied to.

        """

        res = self._op * grid_fun.coefficients
        return res.reshape(self._component_count, -1, order='F')

    def __is_compatible(self, other):
        import numpy as np

        return (self.component_count == other.component_count and
                np.linalg.norm(self.evaluation_points -
                               other.evaluation_points, ord=np.inf) == 0 and
                self.space.is_compatible(other.space))

    def __add__(self, other):

        if not self.__is_compatible(other):
            raise ValueError("Potential operators not compatible.")

        return PotentialOperator(
            self.discrete_operator + other.discrete_operator,
            self.component_count,
            self.space, self.evaluation_points)

    def __mul__(self, obj):
        import numpy as np
        from bempp.api import GridFunction

        if not isinstance(self, PotentialOperator):
            return obj * self

        if np.isscalar(obj):
            return PotentialOperator(obj * self.discrete_operator,
                                     self.component_count,
                                     self.space, self.evaluation_points, )
        elif isinstance(obj, GridFunction):
            return self.evaluate(obj)
        else:
            raise NotImplementedError(
                "Cannot multiply with object of type %s", str(type(obj)))

    def __neg__(self):

        return self.__mul__(-1.0)

    def __sub__(self, other):

        return self.__add__(-other)

    @property
    def space(self):
        """Return the underlying function space."""
        return self._space

    @property
    def component_count(self):
        """Number of components of the potential (1 for scalar potentials)."""
        return self._component_count

    @property
    def evaluation_points(self):
        """Return the evaluation points."""
        return self._evaluation_points

    @property
    def discrete_operator(self):
        """Return the underlying discrete operator."""
        return self._op

    @property
    def dtype(self):
        """Data type of the potential."""
        return self._op.dtype
