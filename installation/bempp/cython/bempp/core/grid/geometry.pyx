from bempp.core.utils cimport Matrix
from bempp.core.utils cimport RowVector
from bempp.core.utils.eigen cimport eigen_matrix_to_np_float64
from bempp.core.utils.eigen cimport eigen_row_vector_to_np_float64
from bempp.core.utils.eigen cimport np_to_eigen_matrix_float64
from cython.operator cimport dereference as deref
import numpy as _np
cimport numpy as _np


cdef class Geometry0:

    def integration_elements(self, local_coordinates):
        """Return the integration elements for the given local coordinates."""

        cdef Matrix[double] local_vector = np_to_eigen_matrix_float64(local_coordinates)
        cdef RowVector[double] int_elements

        self.impl_.getIntegrationElements(local_vector, int_elements)

        return eigen_row_vector_to_np_float64(int_elements)

    def normals(self, local_coordinates):
        """Return normal directions for the given local coordinates."""

        cdef Matrix[double] local_vector = np_to_eigen_matrix_float64(local_coordinates)
        cdef Matrix[double] normal_directions

        self.impl_.getNormals(local_vector, normal_directions)

        return eigen_matrix_to_np_float64(normal_directions)

    def local2global(self, local_coordinates):
        """Return the global coordinates associated with the given local coordinates."""

        cdef Matrix[double] local_vector = np_to_eigen_matrix_float64(local_coordinates)
        cdef Matrix[double] global_coords

        self.impl_.local2global(local_vector, global_coords)

        return eigen_matrix_to_np_float64(global_coords)

    def jacobians_transposed(self, local_coordinates):
        """Return the transposed Jacobians associated with the given local coordinates."""

        cdef Matrix[double] local_vector = np_to_eigen_matrix_float64(local_coordinates)
        cdef vector[Matrix[double]] jacobians_t

        self.impl_.getJacobiansTransposed(local_vector, jacobians_t)

        return [eigen_matrix_to_np_float64(mat) for mat in jacobians_t]

    def jacobian_inverses_transposed(self, local_coordinates):
        """Return the inverse transposed Jacobians associated with the given local coordinates."""

        cdef Matrix[double] local_vector = np_to_eigen_matrix_float64(local_coordinates)
        cdef vector[Matrix[double]] jacobians_t

        self.impl_.getJacobianInversesTransposed(local_vector, jacobians_t)

        return [eigen_matrix_to_np_float64(mat) for mat in jacobians_t]

    property corners:
        """Corners of entity"""
        def __get__(self):

            cdef Matrix[double]* c = new Matrix[double](self.dim_world,self.corner_count)
            cdef _np.ndarray res 
            self.impl_.getCorners(deref(c))
            res = eigen_matrix_to_np_float64(deref(c))
            del c
            return res

    property corner_count:
        """Number of corners of element"""
        def __get__(self):
            return self.impl_.cornerCount()

    property affine:
        """Return if element is affine"""
        def __get__(self):
            return self.impl_.affine()

    property dim:
        """" Dimension of the entity. """
        def __get__(self):
            return self.impl_.dim()

    property dim_world:
        """ Dimension of the space containing the entity. """
        def __get__(self):
            return self.impl_.dimWorld()

    property volume:
        """ Return volume of the entity. """
        def __get__(self):
            return self.impl_.volume()

cdef class Geometry1:

    def integration_elements(self, local_coordinates):
        """Return the integration elements for the given local coordinates."""

        cdef Matrix[double] local_vector = np_to_eigen_matrix_float64(local_coordinates)
        cdef RowVector[double] int_elements

        self.impl_.getIntegrationElements(local_vector, int_elements)

        return eigen_row_vector_to_np_float64(int_elements)

    property corners:
        """Corners of entity"""
        def __get__(self):

            cdef Matrix[double]* c = new Matrix[double](self.dim_world,self.corner_count)
            cdef _np.ndarray res 
            self.impl_.getCorners(deref(c))
            res = eigen_matrix_to_np_float64(deref(c))
            del c
            return res

    property corner_count:
        """Number of corners of element"""
        def __get__(self):
            return self.impl_.cornerCount()

    property affine:
        """Return if element is affine"""
        def __get__(self):
            return self.impl_.affine()

    property dim:
        """" Dimension of the entity. """
        def __get__(self):
            return self.impl_.dim()

    property dim_world:
        """ Dimension of the space containing the entity. """
        def __get__(self):
            return self.impl_.dimWorld()

    property volume:
        """ Return volume of the entity. """
        def __get__(self):
            return self.impl_.volume()

cdef class Geometry2:
    
    def integration_elements(self, local_coordinates):
        """Return the integration elements for the given local coordinates."""

        cdef Matrix[double] local_vector = np_to_eigen_matrix_float64(local_coordinates)
        cdef RowVector[double] int_elements

        self.impl_.getIntegrationElements(local_vector, int_elements)

        return eigen_row_vector_to_np_float64(int_elements)

    property corners:
        """Corners of entity"""
        def __get__(self):

            cdef Matrix[double]* c = new Matrix[double](self.dim_world,self.corner_count)
            cdef _np.ndarray res 
            self.impl_.getCorners(deref(c))
            res = eigen_matrix_to_np_float64(deref(c))
            del c
            return res

    property corner_count:
        """Number of corners of element"""
        def __get__(self):
            return self.impl_.cornerCount()

    property affine:
        """Return if element is affine"""
        def __get__(self):
            return self.impl_.affine()

    property dim:
        """" Dimension of the entity. """
        def __get__(self):
            return self.impl_.dim()

    property dim_world:
        """ Dimension of the space containing the entity. """
        def __get__(self):
            return self.impl_.dimWorld()

    property volume:
        """ Return volume of the entity. """
        def __get__(self):
            return self.impl_.volume()
