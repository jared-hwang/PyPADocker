"""Definition of the electric and magnetic far field operators for Maxwell."""
from bempp.api.operators.potential import _common

#pylint: disable=protected-access

@_common.potential_logger
def electric_field(space, evaluation_points, wave_number, parameters=None):
    """Return the Maxwell electric far field operator

    Parameters
    ----------
    space : bempp.api.space.Space
        The function space over which to assemble the potential.
    evaluation_points : numpy.ndarray
        A (3 x N) array of N evaluation points, where each column corresponds to
        the coordinates of one evaluation point.
    wave_number : complex
        Wavenumber of the operator.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.

    """

    import bempp
    from bempp.api.assembly.potential_operator import PotentialOperator
    from bempp.api.assembly.discrete_boundary_operator import \
        GeneralNonlocalDiscreteBoundaryOperator
    from bempp.core.operators.far_field.maxwell import electric_field_ext

    if parameters is None:
        parameters = bempp.api.global_parameters

    return PotentialOperator(GeneralNonlocalDiscreteBoundaryOperator(
        electric_field_ext(space._impl, evaluation_points,
                           wave_number, parameters)),
                             3, space, evaluation_points)


@_common.potential_logger
def magnetic_field(space, evaluation_points, wave_number, parameters=None):
    """Return the Maxwell magnetic far field operator

    Parameters
    ----------
    space : bempp.api.space.Space
        The function space over which to assemble the potential.
    evaluation_points : numpy.ndarray
        A (3 x N) array of N evaluation points, where each column corresponds to
        the coordinates of one evaluation point.
    wave_number : complex
        Wavenumber of the operator.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.

    """

    import bempp
    from bempp.api.assembly.potential_operator import PotentialOperator
    from bempp.api.assembly.discrete_boundary_operator import \
        GeneralNonlocalDiscreteBoundaryOperator
    from bempp.core.operators.far_field.maxwell import magnetic_field_ext

    if parameters is None:
        parameters = bempp.api.global_parameters

    return PotentialOperator(GeneralNonlocalDiscreteBoundaryOperator(
        magnetic_field_ext(space._impl, evaluation_points, wave_number,
                           parameters)),
                             3, space, evaluation_points)
