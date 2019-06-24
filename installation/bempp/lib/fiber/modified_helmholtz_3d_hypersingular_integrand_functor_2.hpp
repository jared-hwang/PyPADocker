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

#ifndef fiber_modified_helmholtz_3d_hypersingular_integrand_functor_2_hpp
#define fiber_modified_helmholtz_3d_hypersingular_integrand_functor_2_hpp

#include "../common/common.hpp"

#include <cassert>
#include "collection_of_3d_arrays.hpp"
#include "geometrical_data.hpp"
#include "conjugate.hpp"

namespace Fiber {

/** \brief Functor evaluating the integrand of the hypersingular operator for
 *  the modified Helmholtz equation in 3D. */
template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class ModifiedHelmholtz3dHypersingularIntegrandFunctor2 {
public:
  typedef BasisFunctionType_ BasisFunctionType;
  typedef KernelType_ KernelType;
  typedef ResultType_ ResultType;
  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  void addGeometricalDependencies(size_t &testGeomDeps,
                                  size_t &trialGeomDeps) const {
    testGeomDeps |= NORMALS;
    trialGeomDeps |= NORMALS;
  }

  template <template <typename T> class CollectionOf2dSlicesOfConstNdArrays>
  ResultType
  evaluate(const ConstGeometricalDataSlice<CoordinateType> &testGeomData,
           const ConstGeometricalDataSlice<CoordinateType> &trialGeomData,
           const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>
               &testTransfValues,
           const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>
               &trialTransfValues,
           const CollectionOf2dSlicesOfConstNdArrays<KernelType> &kernelValues)
      const {
    const int dimWorld = 3;

    // Assert that there are at least two scalar-valued kernels
    assert(kernelValues.size() >= 2);
    assert(kernelValues[0].extent(0) == 1);
    assert(kernelValues[0].extent(1) == 1);
    assert(kernelValues[1].extent(0) == 1);
    assert(kernelValues[1].extent(1) == 1);

    // Assert that there are at least two test and trial transformations
    // (function value and surface curl) of correct dimensions
    assert(testTransfValues.size() >= 2);
    assert(trialTransfValues.size() >= 2);
    _1dSliceOfConst3dArray<BasisFunctionType> testValues = testTransfValues[0];
    _1dSliceOfConst3dArray<BasisFunctionType> trialValues =
        trialTransfValues[0];
    _1dSliceOfConst3dArray<BasisFunctionType> testSurfaceCurls =
        testTransfValues[1];
    _1dSliceOfConst3dArray<BasisFunctionType> trialSurfaceCurls =
        trialTransfValues[1];
    assert(testValues.extent(0) == 1);
    assert(trialValues.extent(0) == 1);
    assert((int)testSurfaceCurls.extent(0) == dimWorld);
    assert((int)trialSurfaceCurls.extent(0) == dimWorld);

    // Let K_0(x, y) = K(x, y) and K_1(x, y) = kappa^2 K(x, y).
    // Return
    // K_0(x, y) curl u*(x) . curl v(y) + K_1(x, y) u*(x) v(y) n(x) . n(y)

    BasisFunctionType dotProduct0 = 0.;
    for (int dim = 0; dim < dimWorld; ++dim)
      dotProduct0 += conjugate(testSurfaceCurls(dim)) * trialSurfaceCurls(dim);
    ResultType term0 = kernelValues[0](0, 0) * dotProduct0;
    CoordinateType dotProduct1 = 0.;
    for (int dim = 0; dim < dimWorld; ++dim)
      dotProduct1 += testGeomData.normal(dim) * trialGeomData.normal(dim);
    ResultType term1 =
        kernelValues[1](0, 0) *
        (dotProduct1 * conjugate(testValues(0)) * trialValues(0));
    return term0 + term1;
  }
};

} // namespace Fiber

#endif
