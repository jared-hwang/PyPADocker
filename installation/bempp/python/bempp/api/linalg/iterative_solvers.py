"""Iterative solver interfaces."""

import scipy.sparse.linalg
import numpy as np
from bempp.api.assembly import GridFunction
from bempp.api.assembly import BoundaryOperator
from bempp.api.assembly.blocked_operator import BlockedOperatorBase

#pylint: disable=invalid-name
#pylint: disable=too-many-arguments
#pylint: disable=too-many-locals

class _it_counter(object):
    """Iteration Counter class."""

    def __init__(self, store_residuals,
                 iteration_is_cg=False, operator=None, rhs=None):
        self._count = 0
        self._store_residuals = store_residuals
        self._residuals = []
        self._iteration_is_cg = iteration_is_cg
        self._operator = operator
        self._rhs = rhs

    def __call__(self, x):
        self._count += 1
        if self._store_residuals:
            if self._iteration_is_cg:
                res = self._rhs - self._operator * x
            else:
                res = x
            self._residuals.append(np.linalg.norm(res))

    @property
    def count(self):
        """Return the number of iterations."""
        return self._count

    @property
    def residuals(self):
        """Return the vector of residuals."""
        return self._residuals


def gmres(A, b, tol=1E-5,
          restart=None, maxiter=None,
          use_strong_form=False, return_residuals=False,
          return_iteration_count=False):
    """Interface to the scipy.sparse.linalg.gmres function.

    This function behaves like the scipy.sparse.linalg.gmres function. But
    instead of a linear operator and a vector b it takes a boundary operator
    and a grid function or a blocked operator and a list of grid functions.
    The result is returned as a grid function or as a list of grid functions
    in the correct spaces.

    """
    if isinstance(A, BoundaryOperator):
        return _gmres_single_op_imp(
            A, b, tol, restart, maxiter, use_strong_form, return_residuals,
            return_iteration_count)

    if isinstance(A, BlockedOperatorBase):
        return _gmres_block_op_imp(
            A, b, tol, restart, maxiter, use_strong_form, return_residuals,
            return_iteration_count)

    raise ValueError(
        "A must be a BoundaryOperator or BlockedBoundaryOperator")



def cg(A, b, tol=1E-5, maxiter=None,
       use_strong_form=False, return_residuals=False,
       return_iteration_count=False):
    """Interface to the scipy.sparse.linalg.cg function.

    This function behaves like the scipy.sparse.linalg.cg function. But
    instead of a linear operator and a vector b it takes a boundary operator
    and a grid function. The result is returned as a grid function in the
    correct space.

    """
    import bempp.api
    import time

    if not isinstance(A, BoundaryOperator):
        raise ValueError("A must be of type BoundaryOperator")

    if not isinstance(b, GridFunction):
        raise ValueError("b must be of type GridFunction")

    if use_strong_form:
        if not A.range.is_compatible(b.space):
            raise ValueError(
                "The range of A and the domain of A must " +
                "have the same number of unknowns if the strong form is used.")
        A_op = A.strong_form()
        b_vec = b.coefficients
    else:
        A_op = A.weak_form()
        b_vec = b.projections(A.dual_to_range)

    callback = _it_counter(return_residuals, True, A_op, b_vec)
    bempp.api.log("Starting CG iteration")
    start_time = time.time()
    x, info = scipy.sparse.linalg.cg(
        A_op, b_vec,
        tol=tol, maxiter=maxiter, callback=callback)
    end_time = time.time()
    bempp.api.log(
        "CG finished in %i iterations and took %.2E sec." %(
        callback.count, end_time - start_time))

    res_fun = GridFunction(A.domain, coefficients=x.ravel())

    if return_residuals and return_iteration_count:
        return res_fun, info, callback.residuals, callback.count

    if return_residuals:
        return res_fun, info, callback.residuals

    if return_iteration_count:
        return res_fun, info, callback.count

    return res_fun, info

def _gmres_single_op_imp(
        A, b, tol=1E-5,
        restart=None, maxiter=None,
        use_strong_form=False, return_residuals=False,
        return_iteration_count=False):
    """Implementation for single operators."""

    import bempp.api
    import time

    if not isinstance(b, GridFunction):
        raise ValueError("b must be of type GridFunction")

    # Assemble weak form before the logging messages

    if use_strong_form:
        if not A.range.is_compatible(b.space):
            raise ValueError(
                "The range of A and the domain of A must have" +
                "the same number of unknowns if the strong form is used.")
        A_op = A.strong_form()
        b_vec = b.coefficients
    else:
        A_op = A.weak_form()
        b_vec = b.projections(A.dual_to_range)

    callback = _it_counter(return_residuals)

    bempp.api.log("Starting GMRES iteration")
    start_time = time.time()
    x, info = scipy.sparse.linalg.gmres(
        A_op, b_vec,
        tol=tol, restart=restart, maxiter=maxiter, callback=callback)
    end_time = time.time()
    bempp.api.log(
        "GMRES finished in %i iterations and took %.2E sec." %(
        callback.count, end_time - start_time))

    res_fun = GridFunction(A.domain, coefficients=x.ravel())

    if return_residuals and return_iteration_count:
        return res_fun, info, callback.residuals, callback.count

    if return_residuals:
        return res_fun, info, callback.residuals

    if return_iteration_count:
        return res_fun, info, callback.count

    return res_fun, info

def _gmres_block_op_imp(
        A, b, tol=1E-5,
        restart=None, maxiter=None,
        use_strong_form=False, return_residuals=False,
        return_iteration_count=False):
    """Implementation for blocked operators."""

    import bempp.api
    import time
    from bempp.api.assembly.blocked_operator import \
        coefficients_of_grid_function_list, \
        projections_of_grid_function_list, \
        grid_function_list_from_coefficients

    # Assemble weak form before the logging messages

    if use_strong_form:
        b_vec = coefficients_of_grid_function_list(b)
        A_op = A.strong_form()
    else:
        A_op = A.weak_form()
        b_vec = projections_of_grid_function_list(
            b, A.dual_to_range_spaces)

    callback = _it_counter(return_residuals)

    bempp.api.log("Starting GMRES iteration")
    start_time = time.time()
    x, info = scipy.sparse.linalg.gmres(
        A_op, b_vec,
        tol=tol, restart=restart, maxiter=maxiter, callback=callback)
    end_time = time.time()
    bempp.api.log(
        "GMRES finished in %i iterations and took %.2E sec." %(
        callback.count, end_time - start_time))

    res_fun = grid_function_list_from_coefficients(
        x.ravel(), A.domain_spaces)

    if return_residuals and return_iteration_count:
        return res_fun, info, callback.residuals, callback.count

    if return_residuals:
        return res_fun, info, callback.residuals

    if return_iteration_count:
        return res_fun, info, callback.count

    return res_fun, info
