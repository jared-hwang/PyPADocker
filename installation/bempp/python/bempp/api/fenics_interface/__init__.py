"""Routines for the BEM++ FEniCS interface."""

import bempp

#pylint: disable=import-error

if bempp.api.HAVE_DOLFIN:
    from .fenics_operator import FenicsOperator
    from .coupling import boundary_grid_from_fenics_mesh, \
        fenics_to_bempp_trace_data, fenics_space_info
