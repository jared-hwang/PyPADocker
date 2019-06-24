// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "abstract_maxwell_identity_operator.hpp"
#include "../assembly/context.hpp"

#include "../common/boost_make_shared_fwd.hpp"
#include "../fiber/default_test_trial_integral_imp.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/maxwell_3d_test_trial_integrand_functor.hpp"

#include <boost/type_traits/is_complex.hpp>

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
AbstractMaxwellIdentityOperator<BasisFunctionType, ResultType>::
    AbstractMaxwellIdentityOperator(
        const shared_ptr<const Space<BasisFunctionType>> &domain,
        const shared_ptr<const Space<BasisFunctionType>> &range,
        const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
        const std::string &label, int symmetry)
    : Base(domain, range, dualToRange, label, symmetry) {
  typedef Fiber::Maxwell3dTestTrialIntegrandFunctor<
      BasisFunctionType, ResultType> IntegrandFunctor;
  typedef Fiber::DefaultTestTrialIntegral<IntegrandFunctor> Integral;
  m_integral.reset(new Integral((IntegrandFunctor())));
}

template <typename BasisFunctionType, typename ResultType>
AbstractMaxwellIdentityOperator<
    BasisFunctionType, ResultType>::~AbstractMaxwellIdentityOperator() {}

template <typename BasisFunctionType, typename ResultType>
const typename AbstractMaxwellIdentityOperator<
    BasisFunctionType, ResultType>::CollectionOfBasisTransformations &
AbstractMaxwellIdentityOperator<BasisFunctionType,
                                ResultType>::testTransformations() const {
  return this->dualToRange()->basisFunctionValue();
}

template <typename BasisFunctionType, typename ResultType>
const typename AbstractMaxwellIdentityOperator<
    BasisFunctionType, ResultType>::CollectionOfBasisTransformations &
AbstractMaxwellIdentityOperator<BasisFunctionType,
                                ResultType>::trialTransformations() const {
  return this->domain()->basisFunctionValue();
}

template <typename BasisFunctionType, typename ResultType>
const typename AbstractMaxwellIdentityOperator<BasisFunctionType,
                                               ResultType>::TestTrialIntegral &
AbstractMaxwellIdentityOperator<BasisFunctionType, ResultType>::integral()
    const {
  return *m_integral;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    AbstractMaxwellIdentityOperator);

} // namespace Bempp
