// Copyright (C) 2011-2012 by the Bem++ Authors
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

#ifndef fiber_default_local_assembler_for_local_operators_on_surfaces_hpp
#define fiber_default_local_assembler_for_local_operators_on_surfaces_hpp

#include "../common/common.hpp"

#include "local_assembler_for_local_operators.hpp"

#include "numerical_quadrature.hpp"
#include "numerical_test_trial_integrator.hpp"
#include "shared_ptr.hpp"
#include "types.hpp"

#include "../common/boost_ptr_map_fwd.hpp"
#include <tbb/concurrent_unordered_map.h>
#include <tbb/mutex.h>
#include <cstring>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace Fiber {

/** \cond FORWARD_DECL */
class OpenClHandler;

template <typename CoordinateType>
class QuadratureDescriptorSelectorForLocalOperators;
template <typename CoordinateType> class SingleQuadratureRuleFamily;
/** \endcond */

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
class DefaultLocalAssemblerForLocalOperatorsOnSurfaces
    : public LocalAssemblerForLocalOperators<ResultType> {
public:
  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  DefaultLocalAssemblerForLocalOperatorsOnSurfaces(
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
          &testShapesets,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
          &trialShapesets,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const TestTrialIntegral<BasisFunctionType, ResultType>>
          &integral,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const shared_ptr<const QuadratureDescriptorSelectorForLocalOperators<
          CoordinateType>> &quadDescSelector,
      const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>>
          &quadRuleFamily);
  ~DefaultLocalAssemblerForLocalOperatorsOnSurfaces();

  virtual void evaluateLocalWeakForms(const std::vector<int> &elementIndices,
                                      std::vector<Matrix<ResultType>> &result);

private:
  /** \cond PRIVATE */
  const TestTrialIntegrator<BasisFunctionType, ResultType> &
  selectIntegrator(int elementIndex);

  const TestTrialIntegrator<BasisFunctionType, ResultType> &
  getIntegrator(const SingleQuadratureDescriptor &desc);

private:
  // typedef boost::ptr_map<SingleQuadratureDescriptor,
  //                       TestTrialIntegrator<BasisFunctionType, ResultType>>
  //    IntegratorMap;
  typedef tbb::concurrent_unordered_map<
      SingleQuadratureDescriptor,
      TestTrialIntegrator<BasisFunctionType, ResultType> *> IntegratorMap;
  typedef DefaultLocalAssemblerForOperatorsOnSurfacesUtilities<
      BasisFunctionType> Utilities;

private:
  shared_ptr<const GeometryFactory> m_geometryFactory;
  shared_ptr<const RawGridGeometry<CoordinateType>> m_rawGeometry;
  shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
      m_testShapesets;
  shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
      m_trialShapesets;
  shared_ptr<const CollectionOfShapesetTransformations<CoordinateType>>
      m_testTransformations;
  shared_ptr<const CollectionOfShapesetTransformations<CoordinateType>>
      m_trialTransformations;
  shared_ptr<const TestTrialIntegral<BasisFunctionType, ResultType>> m_integral;
  shared_ptr<const OpenClHandler> m_openClHandler;
  shared_ptr<const QuadratureDescriptorSelectorForLocalOperators<
      CoordinateType>> m_quadDescSelector;
  shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>> m_quadRuleFamily;

  IntegratorMap m_testTrialIntegrators;
  mutable tbb::mutex m_integratorCreationMutex;
  /** \endcond */
};

} // namespace Fiber

#include "default_local_assembler_for_local_operators_on_surfaces_imp.hpp"

#endif
