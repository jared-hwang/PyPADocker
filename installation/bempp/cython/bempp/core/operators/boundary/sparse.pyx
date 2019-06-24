from bempp.core.assembly.abstract_boundary_operator cimport ElementaryLocalOperator
from bempp.core.assembly.abstract_boundary_operator cimport c_ElementaryLocalOperator
from bempp.core.space.space cimport c_Space, Space
from bempp.core.utils cimport shared_ptr
from bempp.core.utils import _convert_to_bytes
from bempp.core.utils.enum_types cimport SymmetryMode, symmetry_mode
from bempp.core.utils cimport ParameterList, c_ParameterList
from libcpp.string cimport string
from cython.operator cimport dereference as deref


cdef extern from "bempp/operators/sparse_operators.hpp" namespace "Bempp":
    shared_ptr[const c_ElementaryLocalOperator] identity_operator "Bempp::identityOperator<double,double>"(
            const c_ParameterList&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            string label, SymmetryMode symmetry)
    shared_ptr[const c_ElementaryLocalOperator] maxwell_identity_operator "Bempp::maxwellIdentityOperator<double,double>"(
            const c_ParameterList&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            string label, SymmetryMode symmetry)
    shared_ptr[const c_ElementaryLocalOperator] laplace_beltrami_operator "Bempp::laplaceBeltramiOperator<double,double>"(
            const c_ParameterList&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            string label, SymmetryMode symmetry)

cdef extern from "bempp/core/operators/boundary/support_operators.hpp" namespace "Bempp":
    shared_ptr[const c_ElementaryLocalOperator] curl_value_local_operator(
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            int component)
    shared_ptr[const c_ElementaryLocalOperator] value_times_normal_local_operator(
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            int component)
    shared_ptr[const c_ElementaryLocalOperator] vector_value_times_scalar_local_operator(
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            int component)
    shared_ptr[const c_ElementaryLocalOperator] div_times_scalar_local_operator(
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&)
    shared_ptr[const c_ElementaryLocalOperator] div_times_div_local_operator(
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&)
    shared_ptr[const c_ElementaryLocalOperator] grad_times_hcurl_value_local_operator(
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&)
    shared_ptr[const c_ElementaryLocalOperator] hcurl_times_hcurl_value_local_operator(
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&)
    shared_ptr[const c_ElementaryLocalOperator] hcurl_curl_times_curl_local_operator(
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&)


def identity_ext(
        ParameterList parameters,
        Space domain,
        Space range,
        Space dual_to_range,
        object label='', object symmetry='no_symmetry'):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()
    op.impl_.assign(identity_operator(
        deref(parameters.impl_),domain.impl_, range.impl_, dual_to_range.impl_,
        _convert_to_bytes(label), symmetry_mode(_convert_to_bytes(symmetry))))
    return op

def maxwell_identity_ext(
        ParameterList parameters,
        Space domain,
        Space range,
        Space dual_to_range,
        object label='', object symmetry='no_symmetry'):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()
    op.impl_.assign(maxwell_identity_operator(
        deref(parameters.impl_),domain.impl_, range.impl_, dual_to_range.impl_,
        _convert_to_bytes(label), symmetry_mode(_convert_to_bytes(symmetry))))
    return op

def div_times_div_ext(
        Space domain,
        Space range,
        Space dual_to_range):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()
    op.impl_.assign(div_times_div_local_operator(domain.impl_, range.impl_, dual_to_range.impl_))
    return op


def laplace_beltrami_ext(
        ParameterList parameters,
        Space domain,
        Space range,
        Space dual_to_range,
        object label='', object symmetry='no_symmetry'):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()
    op.impl_.assign(laplace_beltrami_operator(
        deref(parameters.impl_),domain.impl_, range.impl_, dual_to_range.impl_,
        _convert_to_bytes(label), symmetry_mode(_convert_to_bytes(symmetry))))
    return op

def curl_value_ext(Space domain, Space range, Space dual_to_range,
                                  int component):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()
    op.impl_.assign(curl_value_local_operator(domain.impl_, range.impl_, dual_to_range.impl_,
                    component))
    return op

def value_times_normal_ext(Space domain, Space range, Space dual_to_range,
                           int component):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()
    op.impl_.assign(value_times_normal_local_operator(domain.impl_, range.impl_, dual_to_range.impl_,
                    component))
    return op

def vector_value_times_scalar_ext(Space domain, Space range, Space dual_to_range,
                           int component):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()
    op.impl_.assign(vector_value_times_scalar_local_operator(domain.impl_, range.impl_, dual_to_range.impl_,
                                                      component))
    return op

def div_times_scalar_ext(Space domain, Space range, Space dual_to_range):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()
    op.impl_.assign(div_times_scalar_local_operator(domain.impl_, range.impl_, dual_to_range.impl_))
    return op

def grad_times_hcurl_value_ext(Space domain, Space range, Space dual_to_range):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()
    op.impl_.assign(grad_times_hcurl_value_local_operator(domain.impl_, range.impl_, dual_to_range.impl_))
    return op

def hcurl_times_hcurl_value_ext(Space domain, Space range, Space dual_to_range):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()
    op.impl_.assign(hcurl_times_hcurl_value_local_operator(domain.impl_, range.impl_, dual_to_range.impl_))
    return op

def hcurl_curl_times_curl_ext(Space domain, Space range, Space dual_to_range):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()
    op.impl_.assign(hcurl_curl_times_curl_local_operator(domain.impl_, range.impl_, dual_to_range.impl_))
    return op
