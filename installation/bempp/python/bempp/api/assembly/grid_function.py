# pylint: disable-msg=too-many-arguments
"""Definition of Grid functions in BEM++"""


class GridFunction(object):
    """
    Representation of functions on a grid.

    Attributes
    ----------
    coefficients : np.ndarray
        Return or set the vector of coefficients.
    component_count : int
        Return the number of components of the grid
        function values.
    space : bemp.api.space.Space
        Return the space over which the GridFunction is defined.
    grid : bempp.api.grid.Grid
        Return the underlying grid.
    parameters : bempp.api.ParameterList
        Return the set of parameters.
    representation : string
        Return 'primal' if the coefficients of the Gridfunction
        are known. Return 'dual' if only the coefficients in the
        dual space are known.

    """

    def __init__(self, space, dual_space=None, fun=None, coefficients=None,
                 projections=None, parameters=None):
        """
        Construct a grid function.

        A grid function can be initialized in three different ways.

        1. By providing a Python callable. Any Python callable of the
           following form is valid.::

                callable(x,n,domain_index,result)

           Here, x, n, and result are all numpy arrays. x contains the current
           evaluation point, n the associated outward normal direction and
           result is a numpy array that will store the result of the Python
           callable. The variable domain_index stores the index of the
           subdomain on which x lies (default 0). This makes it possible to
           define different functions for different subdomains.

           The following example defines input data that is the inner product
           of the coordinate x with the normal direction n.::

                fun(x,n,domain_index,result):
                    result[0] =  np.dot(x,n)

        2. By providing a vector of coefficients at the nodes. This is
           preferable if the coefficients of the data are coming from an
           external code.

        3. By providing a vector of projection data and a corresponding
           dual space.

        Parameters
        ----------
        space : bempp.api.space.Space
            The space over which the GridFunction is defined.
        dual_space : bempp.api.Space
            A representation of the dual space. If not specified
            then space == dual_space is assumed (optional).
        fun : callable
            A Python function from which the GridFunction is constructed
            (optional).
        coefficients : np.ndarray
            A 1-dimensional array with the coefficients of the GridFunction
            at the interpolatoin points of the space (optional).
        projections : np.ndarray
            A 1-dimensional array with the projections of the GridFunction
            onto a dual space (optional).
        parameters : bempp.api.ParameterList
            A ParameterList object used for the assembly of
            the GridFunction (optional).

        Notes
        -----
        * Only one of projections, coefficients, or fun is allowed as parameter.

        Examples
        --------
        To create a GridFunction from a Python callable my_fun use

        >>> grid_function = GridFunction(space, fun=my_fun)

        To create a GridFunction from a vector of coefficients coeffs use

        >>> grid_function = GridFunction(space,coefficients=coeffs)

        To create a GridFunction from a vector of projections proj use

        >>> grid_function = GridFunction(
                space,dual_space=dual_space, projections=proj)

        """
        import bempp.api

        if space is None:
            raise ValueError("space must not be None.")

        if parameters is None:
            parameters = bempp.api.global_parameters

        if (sum([1 for a in [fun, coefficients, projections] if a is not None])
                != 1):
            raise ValueError(
                "Exactly one of 'fun', 'coefficients' or 'projections' must " +
                "be given.")

        self._coefficients = None
        self._dual_coefficients = None
        self._space = space
        self._parameters = parameters

        if dual_space is not None:
            self._dual_space = dual_space
        else:
            self._dual_space = space

        if coefficients is not None:
            self.coefficients = coefficients
            self._representation = 'primal'

        if fun is not None:
            #pylint: disable=no-name-in-module
            from bempp.core.assembly.function_projector import \
                calculate_projection

            #pylint: disable=protected-access
            self._dual_coefficients = calculate_projection(
                parameters, fun, self._dual_space._impl)
            self._representation = 'dual'

        if projections is not None:
            self._dual_coefficients = projections
            self._representation = 'dual'

    def _compute_coefficients(self, projections, dual_space):
        """Compute coefficients from projections."""
        import bempp.api
        import numpy as np
        np_proj = 1.0 * np.asarray(projections).squeeze()
        if np_proj.ndim > 1:
            raise ValueError("'projections' must be a 1-d array.")

        from bempp.api.assembly import InverseSparseDiscreteBoundaryOperator
        ident = bempp.api.operators.boundary.sparse.identity(
            self.space, self.space, dual_space).weak_form()

        inv_ident = InverseSparseDiscreteBoundaryOperator(ident)

        return inv_ident * projections

    def plot(self, mode='element', transformation='real'):
        """
        Plot the grid function.
        
        Attributes
        ----------
        mode : string
            One of 'element' or 'node'. If 'element' is chosen
            the color is determined by the mid-point of the faces
            of the grid. For 'vertices' the vertex values are
            chosen (default: 'element')
        transformation : string or object
            One of 'real', 'imag', 'abs', 'log_abs' or
            'abs_squared' or a callable object. 
            Describes the data transformation
            before plotting. For functions with vector values
            only 'abs', 'log_abs' or 'abs_squared' are allowed.
            If a callable object is given this is applied instead.
            It is important that the callable returns numpy arrays
            with the same number of dimensions as before.
        
        """
        import bempp.api
        from bempp.api.external.viewers import visualize
        visualize(self, mode, transformation)

    def projections(self, dual_space=None):
        """

        Compute the vector of projections onto the
        given dual space.

        Parameters
        ----------
        dual_space : bempp.api.space.Space
            A representation of the dual space. If not specified
            then fun.dual_space is used.

        Returns
        -------
        out : np.ndarray
            A vector of projections onto the dual space.

        """
        import bempp.api

        if dual_space is None:
            dual_space = self._dual_space

        if (dual_space == self._dual_space and
                self._dual_coefficients is not None):
            return self._dual_coefficients

        ident = bempp.api.operators.boundary.sparse.identity(
            self.space, self.space, dual_space).weak_form()
        self._dual_space = dual_space
        self._dual_coefficients = ident * self.coefficients

        return self._dual_coefficients

    def evaluate(self, element, local_coordinates):
        """Evaluate grid function on a single element."""
        import numpy as np
        coefficients = self.coefficients
        # Get global dof ids and weights
        global_dofs, weights = self.space.get_global_dofs(
            element, dof_weights=True)
        dof_values = np.asarray(
            [coefficients[dof] if dof >= 0 else 0 for dof in global_dofs]) * \
                np.asarray(weights)
        return self.space.evaluate_local_basis(
            element, local_coordinates, dof_values)

    def evaluate_on_element_centers(self):
        """Evaluate the grid function on all element centers."""
        import numpy as np

        grid = self.space.grid
        local_coordinates = np.array([[1./3],[1./3]])
        index_set = grid.leaf_view.index_set()

        values = np.zeros((self.component_count, grid.leaf_view.entity_count(0)),
                dtype=self.dtype)
        for element in grid.leaf_view.entity_iterator(0):
            index = index_set.entity_index(element)
            local_values = self.evaluate(element, local_coordinates)
            values[:, index] = local_values.flat
        return values

    def evaluate_on_vertices(self):
        """Evaluate the grid function on all vertices."""
        import numpy as np

        grid = self.space.grid
        local_coordinates = np.array([[0, 1, 0],[0, 0, 1]])
        index_set = grid.leaf_view.index_set()

        values = np.zeros((self.component_count, grid.leaf_view.entity_count(2)),
                dtype=self.dtype)

        for element in grid.leaf_view.entity_iterator(0):
            local_data = self.evaluate(element, local_coordinates)
            for i in range(3):
                index = index_set.sub_entity_index(element, i, 2)
                values[:, index] = local_data[:, i]
        return values

    def evaluate_surface_gradient(self, element, local_coordinates):
        """Evaluate surface gradient of grid function for scalar spaces."""
        import numpy as np
        coefficients = self.coefficients
        global_dofs, weights = self.space.get_global_dofs(
            element, dof_weights=True)
        dof_values = np.asarray(
            [coefficients[dof] if dof >= 0 else 0 for dof in global_dofs]) * \
                np.asarray(weights)
        return self.space.evaluate_surface_gradient(
            element, local_coordinates, dof_values)

    def integrate(self, element=None):
        """Integrate the function over the grid or a single element."""
        from bempp.api.integration import gauss_triangle_points_and_weights
        import numpy as np

        components = self.component_count
        res = np.zeros((components, 1), dtype='float64')
        accuracy_order = self.parameters.quadrature.far.single_order
        points, weights = gauss_triangle_points_and_weights(accuracy_order)

        element_list = [element] if element is not None else list(
            self.grid.leaf_view.entity_iterator(0))

        for element in element_list:
            integration_elements = element.geometry.integration_elements(
                points)
            res += np.sum(
                self.evaluate(element, points) * weights * integration_elements,
                axis=1)
        return res

    def surface_grad_norm(self, element=None):
        """Norm of the surface gradient on a single element or on the mesh."""
        from bempp.api.integration import gauss_triangle_points_and_weights
        import numpy as np

        res = 0
        accuracy_order = self.parameters.quadrature.far.single_order
        points, weights = gauss_triangle_points_and_weights(accuracy_order)

        element_list = [element] if element is not None else list(
            self.grid.leaf_view.entity_iterator(0))

        for element in element_list:
            integration_elements = element.geometry.integration_elements(
                points)
            abs_surface_gradient_square = np.sum(
                np.abs(
                    self.evaluate_surface_gradient(element, points))**2, axis=0)
            res += np.sum(abs_surface_gradient_square *
                          weights * integration_elements)

        return np.sqrt(res)

    def l2_norm(self, element=None):
        """L^2 norm of the function on a single element or on the mesh."""
        import bempp.api
        from bempp.api.integration import gauss_triangle_points_and_weights
        import numpy as np

        # L2-Norm on the whole space
        if element is None:
            mass = bempp.api.operators.boundary.sparse.identity(
                self.space, self.space, self.space,
                parameters=self.parameters).weak_form().sparse_operator
            vec = self.coefficients
            return np.sqrt(np.abs(vec.conjugate().T.dot(mass.dot(vec))))

        # L2-Norm on a single element
        res = 0
        accuracy_order = self.parameters.quadrature.far.single_order
        points, weights = gauss_triangle_points_and_weights(accuracy_order)

        element_list = [element] if element is not None else list(
            self.grid.leaf_view.entity_iterator(0))

        for element in element_list:
            integration_elements = element.geometry.integration_elements(
                points)
            abs_surface_value_squared = np.sum(
                np.abs(self.evaluate(element, points))**2, axis=0)
            res += np.sum(abs_surface_value_squared *
                          weights * integration_elements)

        return np.sqrt(res)

    #pylint: disable=too-many-locals
    def relative_error(self, fun, element=None):
        """
        Relative L^2 error compared to a given analytic function.

        Parameters
        ----------
        fun : callable
            A python callable of the form f(p), where p is a point with
            the three space components x, y, z at which to evaluate the
            function.
        element : bempp.api.grid.entity.Entity
            An entity of codimension 0.

        """
        from bempp.api.integration import gauss_triangle_points_and_weights
        import numpy as np

        global_diff = 0
        fun_l2_norm = 0
        accuracy_order = self.parameters.quadrature.far.single_order
        points, weights = gauss_triangle_points_and_weights(accuracy_order)
        npoints = points.shape[1]

        element_list = [element] if element is not None else list(
            self.grid.leaf_view.entity_iterator(0))

        for element in element_list:
            integration_elements = element.geometry.integration_elements(
                points)
            global_dofs = element.geometry.local2global(points)
            fun_vals = np.zeros(
                (self.component_count, npoints), dtype=self.dtype)

            for j in range(npoints):
                fun_vals[:, j] = fun(global_dofs[:, j])

            diff = np.sum(np.abs(self.evaluate(
                element, points) - fun_vals)**2, axis=0)
            global_diff += np.sum(diff * integration_elements * weights)
            abs_fun_squared = np.sum(np.abs(fun_vals)**2, axis=0)
            fun_l2_norm += np.sum(abs_fun_squared *
                                  integration_elements * weights)

        return np.sqrt(global_diff / fun_l2_norm)

    def __add__(self, other):
        """Add two grid functions."""

        if self.space != other.space:
            raise ValueError("Spaces are not identical.")

        if self.representation == 'dual' and other.representation == 'dual':
            if self.dual_space == other.dual_space:
                return GridFunction(
                    self.space,
                    projections=self.projections() + other.projections(),
                    dual_space=self.dual_space)

        return GridFunction(self.space,
                            coefficients=self.coefficients + other.coefficients)

    def __mul__(self, alpha):
        import numpy as np

        if np.isscalar(alpha):
            if self.representation == 'dual':
                return GridFunction(self.space,
                                    projections=alpha * self._dual_coefficients,
                                    dual_space=self.dual_space,
                                    parameters=self.parameters)
            else:
                return GridFunction(self.space,
                                    coefficients=alpha * self.coefficients,
                                    parameters=self.parameters)
        else:
            return NotImplemented

    def __rmul__(self, alpha):
        import numpy as np

        if np.isscalar(alpha):
            return self * alpha
        else:
            return NotImplemented

    def __div__(self, alpha):
        if not isinstance(self, GridFunction):
            return (1. / alpha) * self

        return self * (1. / alpha)

    def __truediv__(self, alpha):
        return self.__div__(alpha)

    def __neg__(self):
        return self.__mul__(-1.0)

    def __sub__(self, other):
        if self.space != other.space:
            raise ValueError("Spaces are not identical.")

        return self + (-other)

    @property
    def space(self):
        """Return the Space object."""
        return self._space

    @property
    def dual_space(self):
        """Return the dual space."""
        return self._dual_space

    @property
    def representation(self):
        """Return 'dual' or 'primal'."""
        return self._representation

    @property
    def grid(self):
        """Return the underlying grid."""
        return self.space.grid

    @property
    def parameters(self):
        """Return the parameters."""
        return self._parameters

    @property
    def coefficients(self):
        """Return the function coefficients."""
        if self._coefficients is None:
            self._coefficients = self._compute_coefficients(
                self._dual_coefficients, self.dual_space)
            self._representation = 'primal'
        return self._coefficients

    @coefficients.setter
    def coefficients(self, value):
        """Set the coefficients of the grid function."""
        import numpy as np
        np_coeffs = 1.0 * np.asarray(value)
        if np_coeffs.ndim > 1:
            raise ValueError("'coefficients' must be a 1-d array.")
        self._coefficients = np_coeffs
        self._representation = 'primal'
        self._dual_coefficients = None

    @property
    def component_count(self):
        """Return the number of components of the grid function values."""
        return self.space.codomain_dimension

    @property
    def dtype(self):
        """Return the dtype."""
        if self.representation == 'primal':
            return self._coefficients.dtype
        else:
            return self._dual_coefficients.dtype

    @property
    def real(self):
        """Return a new grid function consisting of the real part of this function."""
        import numpy as np
        import bempp.api
        if self.representation == 'primal':
            return bempp.api.GridFunction(
                space=self.space, dual_space=self.dual_space,
                coefficients=np.real(self.coefficients))
        else:
            return bempp.api.GridFunction(
                space=self.space, dual_space=self.dual_space,
                projections=np.real(self.projections()))

    @property
    def imag(self):
        """Return a new grid function consisting of the imaginary part of this function."""
        import numpy as np
        import bempp.api
        if self.representation == 'primal':
            return bempp.api.GridFunction(
                space=self.space, dual_space=self.dual_space,
                coefficients=np.imag(self.coefficients))
        else:
            return bempp.api.GridFunction(
                space=self.space, dual_space=self.dual_space,
                projections=np.imag(self.projections()))

        
    @classmethod
    def from_random(cls, space):
        """Create a random grid function normalized to unit norm. """
        from numpy.random import randn
        ndofs = space.global_dof_count
        fun = cls(space, coefficients=randn(ndofs))
        return fun / fun.l2_norm()

    @classmethod
    def from_ones(cls, space):
        """Create a grid function with all coefficients set to one. """
        from numpy import ones
        ndofs = space.global_dof_count

        return cls(space, coefficients=ones(ndofs))

    @classmethod
    def from_zeros(cls, space):
        """Create a grid function with all coefficients set to one. """
        from numpy import zeros
        ndofs = space.global_dof_count

        return cls(space, coefficients=zeros(ndofs))

