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

#ifndef fiber_modified_maxwell_3d_single_layer_potential_operator_integrand_functor_hpp
#define fiber_modified_maxwell_3d_single_layer_potential_operator_integrand_functor_hpp

#include "../common/common.hpp"

#include "collection_of_2d_arrays.hpp"
#include "collection_of_4d_arrays.hpp"
#include "geometrical_data.hpp"
#include "conjugate.hpp"

#include <cassert>
#include <vector>

namespace Fiber {

template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class ModifiedMaxwell3dSingleLayerPotentialOperatorIntegrandFunctor {
public:
  typedef BasisFunctionType_ BasisFunctionType;
  typedef KernelType_ KernelType;
  typedef ResultType_ ResultType;
  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  void addGeometricalDependencies(size_t &trialGeomDeps) const {
    // Do nothing
  }

  int resultDimension() const { return 3; }

  template <template <typename T> class CollectionOf1dSlicesOfConstNdArrays,
            typename TrialValueType>
  void evaluate(
      const ConstGeometricalDataSlice<CoordinateType> & /* trialGeomData */,
      const CollectionOf2dSlicesOfConst4dArrays<KernelType> &kernelValues,
      const CollectionOf1dSlicesOfConstNdArrays<TrialValueType>
          &weightedTransformedTrialValues,
      std::vector<ResultType> &result) const {
    const int dimWorld = 3;

    // Assert that there are at least two kernels with correct dimensions
    assert(kernelValues.size() >= 2);
    assert(kernelValues[0].extent(0) == 1);
    assert(kernelValues[0].extent(1) == 1);
    assert(kernelValues[1].extent(0) == 3);
    assert(kernelValues[1].extent(1) == 1);

    // Assert that there are at least two test and trial transformations
    // (function value and surface div) of correct dimensions
    assert(weightedTransformedTrialValues.size() >= 2);
    typedef
        typename CollectionOf1dSlicesOfConstNdArrays<TrialValueType>::ConstSlice
            _1dSliceOfConstNdArray;

    _1dSliceOfConstNdArray trialValues = weightedTransformedTrialValues[0];
    _1dSliceOfConstNdArray trialSurfaceDivs = weightedTransformedTrialValues[1];
    assert(trialValues.extent(0) == 3);
    assert(trialSurfaceDivs.extent(0) == 1);

    // Assert that the result vector is three-dimensional
    assert(result.size() == dimWorld);

    for (int dim = 0; dim < dimWorld; ++dim)
      result[dim] = -kernelValues[0](0, 0) * trialValues(dim) +
                    kernelValues[1](dim, 0) * trialSurfaceDivs(0);
  }
};

} // namespace Fiber

#endif
