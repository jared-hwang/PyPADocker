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

#ifndef fiber_default_local_assembler_for_potential_operators_on_surfaces_hpp
#define fiber_default_local_assembler_for_potential_operators_on_surfaces_hpp

#include "../common/common.hpp"

#include "local_assembler_for_potential_operators.hpp"

#include "_2d_array.hpp"
#include "accuracy_options.hpp"
#include "default_local_assembler_for_operators_on_surfaces_utilities.hpp"
#include "element_pair_topology.hpp"
#include "numerical_quadrature.hpp"
#include "parallelization_options.hpp"
#include "shared_ptr.hpp"
#include "kernel_trial_integrator.hpp"
#include "verbosity_level.hpp"

#include <boost/static_assert.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <tbb/concurrent_unordered_map.h>
#include <cstring>
#include <climits>
#include <set>
#include <utility>
#include <vector>

namespace Fiber {

/** \cond FORWARD_DECL */
class OpenClHandler;
template <typename CoordinateType> class CollectionOfShapesetTransformations;
template <typename ValueType> class CollectionOfKernels;
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class KernelTrialIntegral;
template <typename CoordinateType> class RawGridGeometry;
template <typename BasisFunctionType>
class QuadratureDescriptorSelectorForPotentialOperators;
template <typename CoordinateType> class SingleQuadratureRuleFamily;
/** \endcond */

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
class DefaultLocalAssemblerForPotentialOperatorsOnSurfaces
    : public LocalAssemblerForPotentialOperators<ResultType> {
public:
  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  DefaultLocalAssemblerForPotentialOperatorsOnSurfaces(
      const Matrix<CoordinateType> &points,
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
          &trialShapesets,
      const shared_ptr<const CollectionOfKernels<KernelType>> &kernels,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const KernelTrialIntegral<BasisFunctionType, KernelType,
                                                 ResultType>> &integral,
      const ParallelizationOptions &parallelizationOptions,
      VerbosityLevel::Level verbosityLevel,
      const shared_ptr<const QuadratureDescriptorSelectorForPotentialOperators<
          BasisFunctionType>> &quadDescSelector,
      const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>>
          &quadRuleFamily);
  virtual ~DefaultLocalAssemblerForPotentialOperatorsOnSurfaces();

  virtual void evaluateLocalContributions(
      const std::vector<int> &pointIndices, int trialElementIndex,
      LocalDofIndex localTrialDofIndex, std::vector<Matrix<ResultType>> &result,
      CoordinateType nominalDistance = -1.);

  virtual void
  evaluateLocalContributions(int pointIndex, int componentIndex,
                             const std::vector<int> &trialElementIndices,
                             std::vector<Matrix<ResultType>> &result,
                             CoordinateType nominalDistance = -1.);

  virtual void
  evaluateLocalContributions(const std::vector<int> &pointIndices,
                             const std::vector<int> &trialElementIndices,
                             Fiber::_2dArray<Matrix<ResultType>> &result,
                             CoordinateType nominalDistance = -1.);

  virtual int resultDimension() const;

  virtual CoordinateType estimateRelativeScale(CoordinateType minDist) const;

private:
  /** \cond PRIVATE */
  typedef KernelTrialIntegrator<BasisFunctionType, KernelType, ResultType>
      Integrator;
  typedef DefaultLocalAssemblerForOperatorsOnSurfacesUtilities<
      BasisFunctionType> Utilities;

  const Integrator &selectIntegrator(int testElementIndex,
                                     int trialElementIndex,
                                     CoordinateType nominalDistance = -1.);

  int order(int pointIndex, int trialElementIndex,
            CoordinateType nominalDistance) const;

  const Integrator &getIntegrator(const SingleQuadratureDescriptor &index);

  CoordinateType pointElementDistanceSquared(int pointIndex,
                                             int trialElementIndex) const;

  void precalculateElementSizesAndCenters();

private:
  Matrix<CoordinateType> m_points;
  shared_ptr<const GeometryFactory> m_geometryFactory;
  shared_ptr<const RawGridGeometry<CoordinateType>> m_rawGeometry;
  shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
      m_trialShapesets;
  shared_ptr<const CollectionOfKernels<KernelType>> m_kernels;
  shared_ptr<const CollectionOfShapesetTransformations<CoordinateType>>
      m_trialTransformations;
  shared_ptr<const KernelTrialIntegral<BasisFunctionType, KernelType,
                                       ResultType>> m_integral;
  ParallelizationOptions m_parallelizationOptions;
  VerbosityLevel::Level m_verbosityLevel;
  shared_ptr<const QuadratureDescriptorSelectorForPotentialOperators<
      BasisFunctionType>> m_quadDescSelector;
  shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>> m_quadRuleFamily;

  typedef tbb::concurrent_unordered_map<SingleQuadratureDescriptor,
                                        Integrator *> IntegratorMap;
  IntegratorMap m_kernelTrialIntegrators;

  enum { INVALID_INDEX = INT_MAX };
  /** \endcond */
};

} // namespace Fiber

#include "default_local_assembler_for_potential_operators_on_surfaces_imp.hpp"

#endif
