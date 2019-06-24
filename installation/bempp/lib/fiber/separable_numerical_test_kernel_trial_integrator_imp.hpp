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

#include "../common/common.hpp"

#include "bempp/common/config_opencl.hpp"

#include "separable_numerical_test_kernel_trial_integrator.hpp" // To keep IDEs happy

#include "_2d_array.hpp"
#include "_3d_array.hpp"
#include "_4d_array.hpp"

#include "shapeset.hpp"
#include "basis_data.hpp"
#include "conjugate.hpp"
#include "collection_of_shapeset_transformations.hpp"
#include "geometrical_data.hpp"
#include "collection_of_kernels.hpp"
#include "opencl_handler.hpp"
#include "raw_grid_geometry.hpp"
#include "test_kernel_trial_integral.hpp"
#include "types.hpp"
#include "CL/separable_numerical_double_integrator.cl.str"

#include "../common/auto_timer.hpp"

#include <cassert>
#include <memory>

namespace Fiber {

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
SeparableNumericalTestKernelTrialIntegrator<BasisFunctionType, KernelType,
                                            ResultType, GeometryFactory>::
    SeparableNumericalTestKernelTrialIntegrator(
        const Matrix<CoordinateType> &localTestQuadPoints,
        const Matrix<CoordinateType> &localTrialQuadPoints,
        const std::vector<CoordinateType> &testQuadWeights,
        const std::vector<CoordinateType> &trialQuadWeights,
        const GeometryFactory &testGeometryFactory,
        const GeometryFactory &trialGeometryFactory,
        const RawGridGeometry<CoordinateType> &testRawGeometry,
        const RawGridGeometry<CoordinateType> &trialRawGeometry,
        const CollectionOfShapesetTransformations<CoordinateType>
            &testTransformations,
        const CollectionOfKernels<KernelType> &kernels,
        const CollectionOfShapesetTransformations<CoordinateType>
            &trialTransformations,
        const TestKernelTrialIntegral<BasisFunctionType, KernelType, ResultType>
            &integral,
        const OpenClHandler &openClHandler, bool cacheGeometricalData)
    : m_localTestQuadPoints(localTestQuadPoints),
      m_localTrialQuadPoints(localTrialQuadPoints),
      m_testQuadWeights(testQuadWeights), m_trialQuadWeights(trialQuadWeights),
      m_testGeometryFactory(testGeometryFactory),
      m_trialGeometryFactory(trialGeometryFactory),
      m_testRawGeometry(testRawGeometry), m_trialRawGeometry(trialRawGeometry),
      m_testTransformations(testTransformations), m_kernels(kernels),
      m_trialTransformations(trialTransformations), m_integral(integral),
      m_openClHandler(openClHandler),
      m_cacheGeometricalData(cacheGeometricalData) {
  if (localTestQuadPoints.cols() != testQuadWeights.size())
    throw std::invalid_argument(
        "SeparableNumericalTestKernelTrialIntegrator::"
        "SeparableNumericalTestKernelTrialIntegrator(): "
        "numbers of test points and weights do not match");
  if (localTrialQuadPoints.cols() != trialQuadWeights.size())
    throw std::invalid_argument(
        "SeparableNumericalTestKernelTrialIntegrator::"
        "SeparableNumericalTestKernelTrialIntegrator(): "
        "numbers of trial points and weights do not match");

#ifdef WITH_OPENCL
  if (openClHandler.UseOpenCl()) {
    // push integration points to CL device
    clTestQuadPoints =
        openClHandler.pushMatrix<CoordinateType>(localTestQuadPoints);
    clTrialQuadPoints =
        openClHandler.pushMatrix<CoordinateType>(localTrialQuadPoints);
    clTestQuadWeights =
        openClHandler.pushVector<CoordinateType>(testQuadWeights);
    clTrialQuadWeights =
        openClHandler.pushVector<CoordinateType>(trialQuadWeights);
  }
#endif

  if (cacheGeometricalData)
    precalculateGeometricalData();
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
SeparableNumericalTestKernelTrialIntegrator<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::~SeparableNumericalTestKernelTrialIntegrator() {
#ifdef WITH_OPENCL
  if (m_openClHandler.UseOpenCl()) {
    delete clTestQuadPoints;
    delete clTrialQuadPoints;
    delete clTestQuadWeights;
    delete clTrialQuadWeights;
  }
#endif
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void SeparableNumericalTestKernelTrialIntegrator<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::precalculateGeometricalData() {
  size_t testBasisDeps = 0, trialBasisDeps = 0; // ignored in this function
  size_t testGeomDeps = 0, trialGeomDeps = 0;

  m_testTransformations.addDependencies(testBasisDeps, testGeomDeps);
  m_trialTransformations.addDependencies(trialBasisDeps, trialGeomDeps);
  m_kernels.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
  m_integral.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

  precalculateGeometricalDataOnSingleGrid(
      m_localTestQuadPoints, m_testGeometryFactory, m_testRawGeometry,
      testGeomDeps, m_cachedTestGeomData);
  precalculateGeometricalDataOnSingleGrid(
      m_localTrialQuadPoints, m_trialGeometryFactory, m_trialRawGeometry,
      trialGeomDeps, m_cachedTrialGeomData);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void SeparableNumericalTestKernelTrialIntegrator<BasisFunctionType, KernelType,
                                                 ResultType, GeometryFactory>::
    precalculateGeometricalDataOnSingleGrid(
        const Matrix<CoordinateType> &localQuadPoints,
        const GeometryFactory &geometryFactory,
        const RawGridGeometry<CoordinateType> &rawGeometry, size_t geomDeps,
        std::vector<GeometricalData<CoordinateType>> &geomData) {
  geomData.resize(rawGeometry.elementCount());

  typedef typename GeometryFactory::Geometry Geometry;
  std::unique_ptr<Geometry> geometry = geometryFactory.make();
  for (size_t e = 0; e < rawGeometry.elementCount(); ++e) {
    rawGeometry.setupGeometry(e, *geometry);
    geometry->getData(geomDeps, localQuadPoints, geomData[e]);
    if (geomDeps & DOMAIN_INDEX)
      geomData[e].domainIndex = rawGeometry.domainIndex(e);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void SeparableNumericalTestKernelTrialIntegrator<BasisFunctionType, KernelType,
                                                 ResultType, GeometryFactory>::
    integrate(CallVariant callVariant, const std::vector<int> &elementIndicesA,
              int elementIndexB, const Shapeset<BasisFunctionType> &basisA,
              const Shapeset<BasisFunctionType> &basisB,
              LocalDofIndex localDofIndexB,
              const std::vector<Matrix<ResultType> *> &result) const {
  if (m_openClHandler.UseOpenCl()) {
    integrateCl(callVariant, elementIndicesA, elementIndexB, basisA, basisB,
                localDofIndexB, result);
  } else {
    integrateCpu(callVariant, elementIndicesA, elementIndexB, basisA, basisB,
                 localDofIndexB, result);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void SeparableNumericalTestKernelTrialIntegrator<BasisFunctionType, KernelType,
                                                 ResultType, GeometryFactory>::
    integrateCpu(CallVariant callVariant,
                 const std::vector<int> &elementIndicesA, int elementIndexB,
                 const Shapeset<BasisFunctionType> &basisA,
                 const Shapeset<BasisFunctionType> &basisB,
                 LocalDofIndex localDofIndexB,
                 const std::vector<Matrix<ResultType> *> &result) const {
  const int testPointCount = m_localTestQuadPoints.cols();
  const int trialPointCount = m_localTrialQuadPoints.cols();
  const int elementACount = elementIndicesA.size();

  if (result.size() != elementIndicesA.size())
    throw std::invalid_argument(
        "SeparableNumericalTestKernelTrialIntegrator::integrate(): "
        "arrays 'result' and 'elementIndicesA' must have the same number "
        "of elements");
  if (testPointCount == 0 || trialPointCount == 0 || elementACount == 0)
    return;
  // TODO: in the (pathological) case that pointCount == 0 but
  // geometryCount != 0, set elements of result to 0.

  // Evaluate constants

  const int dofCountA = basisA.size();
  const int dofCountB = localDofIndexB == ALL_DOFS ? basisB.size() : 1;
  const int testDofCount = callVariant == TEST_TRIAL ? dofCountA : dofCountB;
  const int trialDofCount = callVariant == TEST_TRIAL ? dofCountB : dofCountA;

  BasisData<BasisFunctionType> testBasisData, trialBasisData;
  GeometricalData<CoordinateType> *testGeomData = &m_testGeomData.local();
  GeometricalData<CoordinateType> *trialGeomData = &m_trialGeomData.local();
  const GeometricalData<CoordinateType> *constTestGeomData = testGeomData;
  const GeometricalData<CoordinateType> *constTrialGeomData = trialGeomData;

  size_t testBasisDeps = 0, trialBasisDeps = 0;
  size_t testGeomDeps = 0, trialGeomDeps = 0;

  m_testTransformations.addDependencies(testBasisDeps, testGeomDeps);
  m_trialTransformations.addDependencies(trialBasisDeps, trialGeomDeps);
  m_kernels.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
  m_integral.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

  typedef typename GeometryFactory::Geometry Geometry;
  std::unique_ptr<Geometry> geometryA, geometryB;
  const RawGridGeometry<CoordinateType> *rawGeometryA = 0, *rawGeometryB = 0;
  if (!m_cacheGeometricalData) {
    if (callVariant == TEST_TRIAL) {
      geometryA = m_testGeometryFactory.make();
      geometryB = m_trialGeometryFactory.make();
      rawGeometryA = &m_testRawGeometry;
      rawGeometryB = &m_trialRawGeometry;
    } else {
      geometryA = m_trialGeometryFactory.make();
      geometryB = m_testGeometryFactory.make();
      rawGeometryA = &m_trialRawGeometry;
      rawGeometryB = &m_testRawGeometry;
    }
  }

  CollectionOf3dArrays<BasisFunctionType> testValues, trialValues;
  CollectionOf4dArrays<KernelType> kernelValues;

  for (size_t i = 0; i < result.size(); ++i) {
    assert(result[i]);
    result[i]->resize(testDofCount, trialDofCount);
  }

  if (!m_cacheGeometricalData)
    rawGeometryB->setupGeometry(elementIndexB, *geometryB);
  if (callVariant == TEST_TRIAL) {
    basisA.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS,
                    testBasisData);
    basisB.evaluate(trialBasisDeps, m_localTrialQuadPoints, localDofIndexB,
                    trialBasisData);
    if (m_cacheGeometricalData)
      constTrialGeomData = &m_cachedTrialGeomData[elementIndexB];
    else {
      geometryB->getData(trialGeomDeps, m_localTrialQuadPoints, *trialGeomData);
      if (trialGeomDeps & DOMAIN_INDEX)
        trialGeomData->domainIndex = rawGeometryB->domainIndex(elementIndexB);
    }
    m_trialTransformations.evaluate(trialBasisData, *constTrialGeomData,
                                    trialValues);
  } else {
    basisA.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS,
                    trialBasisData);
    basisB.evaluate(testBasisDeps, m_localTestQuadPoints, localDofIndexB,
                    testBasisData);
    if (m_cacheGeometricalData)
      constTestGeomData = &m_cachedTestGeomData[elementIndexB];
    else {
      geometryB->getData(testGeomDeps, m_localTestQuadPoints, *testGeomData);
      if (testGeomDeps & DOMAIN_INDEX)
        testGeomData->domainIndex = rawGeometryB->domainIndex(elementIndexB);
    }
    m_testTransformations.evaluate(testBasisData, *constTestGeomData,
                                   testValues);
  }

  // Iterate over the elements
  for (int indexA = 0; indexA < elementACount; ++indexA) {
    const int elementIndexA = elementIndicesA[indexA];
    if (!m_cacheGeometricalData)
      rawGeometryA->setupGeometry(elementIndexA, *geometryA);
    if (callVariant == TEST_TRIAL) {
      if (m_cacheGeometricalData)
        constTestGeomData = &m_cachedTestGeomData[elementIndicesA[indexA]];
      else {
        geometryA->getData(testGeomDeps, m_localTestQuadPoints, *testGeomData);
        if (testGeomDeps & DOMAIN_INDEX)
          testGeomData->domainIndex = rawGeometryA->domainIndex(elementIndexA);
      }
      m_testTransformations.evaluate(testBasisData, *constTestGeomData,
                                     testValues);
    } else {
      if (m_cacheGeometricalData)
        constTrialGeomData = &m_cachedTrialGeomData[elementIndicesA[indexA]];
      else {
        geometryA->getData(trialGeomDeps, m_localTrialQuadPoints,
                           *trialGeomData);
        if (trialGeomDeps & DOMAIN_INDEX)
          trialGeomData->domainIndex = rawGeometryA->domainIndex(elementIndexA);
      }
      m_trialTransformations.evaluate(trialBasisData, *constTrialGeomData,
                                      trialValues);
    }

    m_kernels.evaluateOnGrid(*constTestGeomData, *constTrialGeomData,
                             kernelValues);
    m_integral.evaluateWithTensorQuadratureRule(
        *constTestGeomData, *constTrialGeomData, testValues, trialValues,
        kernelValues, m_testQuadWeights, m_trialQuadWeights, *result[indexA]);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void SeparableNumericalTestKernelTrialIntegrator<BasisFunctionType, KernelType,
                                                 ResultType, GeometryFactory>::
    integrateCl(CallVariant callVariant,
                const std::vector<int> &elementIndicesA, int elementIndexB,
                const Shapeset<BasisFunctionType> &basisA,
                const Shapeset<BasisFunctionType> &basisB,
                LocalDofIndex localDofIndexB,
                const std::vector<Matrix<ResultType> *> &result) const {
  //#ifdef WITH_OPENCL
  //  // DEBUG: test latency
  //  {
  //    int bufsize = 10;
  //    std::vector<ResultType> buf(bufsize);
  //    tbb::tick_count t0 = tbb::tick_count::now();
  //    cl::Buffer *clbuf = m_openClHandler.pushVector<ResultType> (buf);
  //    tbb::tick_count::interval_t dt1 = tbb::tick_count::now()-t0;

  //    t0 = tbb::tick_count::now();
  //    m_openClHandler.pullVector<ResultType> (*clbuf, buf, bufsize);
  //    tbb::tick_count::interval_t dt2 = tbb::tick_count::now()-t0;

  //    static int callcount = 0;
  //    std::cout << callcount++ << "\tpush: " << dt1.seconds() << "\tpull: " <<
  // dt2.seconds() << std::endl;
  //  }

  //  //    tbb::tick_count t_start, t_end, t0, t1;
  //  //    tbb::tick_count::interval_t dt_buf, dt_kern, dt_prog, dt_pull;

  //  //    t_start = tbb::tick_count::now();

  //    // Temporary code. TODO: add support for multiple-term expressions.
  //    if (!m_testTransformation.isTrivial() ||
  // !m_trialTransformation.isTrivial())
  //        throw
  // std::runtime_error("SeparableNumericalTestKernelTrialIntegrator::"
  //                                 "integrateCl(): multiple-term expression
  // lists are "
  //                                 "not supported in the OpenCL mode yet");
  //    const Expression<CoordinateType>& m_testExpression =
  // m_testTransformation.term(0);
  //    const Expression<CoordinateType>& m_trialExpression =
  // m_trialTransformation.term(0);

  //    const int testPointCount = m_localTestQuadPoints.cols();
  //    const int trialPointCount = m_localTrialQuadPoints.cols();
  //    const int elementACount = elementIndicesA.size();
  //    const int pointDim = m_localTestQuadPoints.rows();
  //    const int meshDim = m_openClHandler.meshGeom().size.dim;

  //    if (testPointCount == 0 || trialPointCount == 0 || elementACount == 0)
  //        return;
  //    // TODO: in the (pathological) case that pointCount == 0 but
  //    // geometryCount != 0, set elements of result to 0.

  //    // Evaluate constants
  //    const int testComponentCount = m_testExpression.codomainDimension();
  //    const int trialComponentCount = m_trialExpression.codomainDimension();
  //    const int dofCountA = basisA.size();
  //    const int dofCountB = localDofIndexB == ALL_DOFS ? basisB.size() : 1;
  //    const int testDofCount = callVariant == TEST_TRIAL ? dofCountA :
  // dofCountB;
  //    const int trialDofCount = callVariant == TEST_TRIAL ? dofCountB :
  // dofCountA;

  //    const int kernelRowCount = m_kernels.codomainDimension();
  //    const int kernelColCount = m_kernels.domainDimension();

  //    // Assert that the kernel tensor dimensions are compatible
  //    // with the number of components of the functions

  //    // TODO: This will need to be modified once we allow scalar-valued
  // kernels
  //    // (treated as if they were multiplied by the unit tensor) with
  //    // vector-valued functions
  //    assert(testComponentCount == kernelRowCount);
  //    assert(kernelColCount == trialComponentCount);

  //    int argIdx;
  //    int testBasisDeps = 0, trialBasisDeps = 0;
  //    int testGeomDeps = INTEGRATION_ELEMENTS;
  //    int trialGeomDeps = INTEGRATION_ELEMENTS;

  //    cl::Buffer *clElementIndicesA;
  //    cl::Buffer *clGlobalTrialPoints;
  //    cl::Buffer *clGlobalTestPoints;
  //    cl::Buffer *clGlobalTrialNormals;
  //    cl::Buffer *clTrialValues;
  //    cl::Buffer *clTestValues;
  //    cl::Buffer *clTrialIntegrationElements;
  //    cl::Buffer *clTestIntegrationElements;
  //    cl::Buffer *clResult;

  //    m_testExpression.addDependencies(testBasisDeps, testGeomDeps);
  //    m_trialExpression.addDependencies(trialBasisDeps, trialGeomDeps);
  //    m_kernels.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

  //    result.set_size(testDofCount, trialDofCount, elementACount);

  //    //    t0 = tbb::tick_count::now();
  //    clElementIndicesA = m_openClHandler.pushVector<int> (elementIndicesA);
  //    clResult = m_openClHandler.createBuffer<ResultType>
  // (testDofCount*trialDofCount*elementACount,
  //							 CL_MEM_WRITE_ONLY);
  //    //    dt_buf += tbb::tick_count::now()-t0;

  //    // Build the OpenCL program
  //    //    t0 = tbb::tick_count::now();
  //    cl::Program::Sources sources;
  //    sources.push_back (m_openClHandler.initStr());
  //    sources.push_back (basisA.clCodeString(true));
  //    sources.push_back (basisB.clCodeString(false));
  //    sources.push_back (m_kernels.evaluateClCode());
  //    sources.push_back (clStrIntegrateRowOrCol());
  //    m_openClHandler.loadProgramFromStringArray (sources);
  //    //    dt_prog += tbb::tick_count::now()-t0;

  //    // Call the CL kernels to map the trial and test quadrature points
  //    if (callVariant == TEST_TRIAL)
  //    {
  //        //        t0 = tbb::tick_count::now();
  //        clGlobalTestPoints = m_openClHandler.createBuffer<CoordinateType>(
  //            elementACount*testPointCount*meshDim, CL_MEM_READ_WRITE);
  //	clTestIntegrationElements =
  // m_openClHandler.createBuffer<CoordinateType>(
  //	    elementACount*testPointCount, CL_MEM_READ_WRITE);
  //	//	dt_buf += tbb::tick_count::now()-t0;
  //	cl::Kernel &clMapTest = m_openClHandler.setKernel
  //("clMapPointsToElements");
  //	argIdx = m_openClHandler.SetGeometryArgs (clMapTest, 0);
  //	clMapTest.setArg (argIdx++, *clTestQuadPoints);
  //	clMapTest.setArg (argIdx++, testPointCount);
  //	clMapTest.setArg (argIdx++, pointDim);
  //	clMapTest.setArg (argIdx++, *clElementIndicesA);
  //	clMapTest.setArg (argIdx++, elementACount);
  //	clMapTest.setArg (argIdx++, *clGlobalTestPoints);
  //	clMapTest.setArg (argIdx++, *clTestIntegrationElements);

  //	//        t0 = tbb::tick_count::now();
  //	m_openClHandler.enqueueKernel (cl::NDRange(elementACount,
  // testPointCount));
  //	//	dt_kern += tbb::tick_count::now()-t0;

  //	//        t0 = tbb::tick_count::now();
  //        clGlobalTrialPoints = m_openClHandler.createBuffer<CoordinateType> (
  //	    trialPointCount*meshDim, CL_MEM_READ_WRITE);
  //        clGlobalTrialNormals = m_openClHandler.createBuffer<CoordinateType>
  // (
  //	    trialPointCount*meshDim, CL_MEM_READ_WRITE);
  //	clTrialIntegrationElements =
  // m_openClHandler.createBuffer<CoordinateType>(
  //	    trialPointCount, CL_MEM_READ_WRITE);
  //	//	dt_buf += tbb::tick_count::now()-t0;
  //	cl::Kernel &clMapTrial = m_openClHandler.setKernel
  //("clMapPointsAndNormalsToElement");
  //	argIdx = m_openClHandler.SetGeometryArgs (clMapTrial, 0);
  //	clMapTrial.setArg (argIdx++, *clTestQuadPoints);
  //	clMapTrial.setArg (argIdx++, trialPointCount);
  //	clMapTrial.setArg (argIdx++, pointDim);
  //	clMapTrial.setArg (argIdx++, elementIndexB);
  //	clMapTrial.setArg (argIdx++, *clGlobalTrialPoints);
  //	clMapTrial.setArg (argIdx++, *clGlobalTrialNormals);
  //	clMapTrial.setArg (argIdx++, *clTrialIntegrationElements);
  //	//        t0 = tbb::tick_count::now();
  //	m_openClHandler.enqueueKernel (cl::NDRange(trialPointCount));
  //	//	dt_kern += tbb::tick_count::now()-t0;

  //	//        t0 = tbb::tick_count::now();
  //	clTestValues = m_openClHandler.createBuffer<BasisFunctionType> (
  //	    elementACount*testPointCount*testDofCount, CL_MEM_READ_WRITE);
  //	//	dt_buf += tbb::tick_count::now()-t0;
  //	cl::Kernel &clBasisTest = m_openClHandler.setKernel
  //("clBasisfAElements");
  //	argIdx = m_openClHandler.SetGeometryArgs (clBasisTest, 0);
  //	clBasisTest.setArg (argIdx++, *clElementIndicesA);
  //	clBasisTest.setArg (argIdx++, elementACount);
  //	clBasisTest.setArg (argIdx++, *clTestQuadPoints);
  //	clBasisTest.setArg (argIdx++, testPointCount);
  //	clBasisTest.setArg (argIdx++, pointDim);
  //	clBasisTest.setArg (argIdx++, testDofCount);
  //	clBasisTest.setArg (argIdx++, *clTestValues);
  //	//        t0 = tbb::tick_count::now();
  //	m_openClHandler.enqueueKernel (cl::NDRange(elementACount,
  // testPointCount));
  //	//	dt_kern += tbb::tick_count::now()-t0;

  //	//        t0 = tbb::tick_count::now();
  //	clTrialValues = m_openClHandler.createBuffer<BasisFunctionType> (
  //	    trialPointCount*trialDofCount, CL_MEM_READ_WRITE);
  //	//	dt_buf += tbb::tick_count::now()-t0;
  //	cl::Kernel &clBasisTrial = m_openClHandler.setKernel
  //("clBasisfBElement");
  //	argIdx = m_openClHandler.SetGeometryArgs (clBasisTrial, 0);
  //	clBasisTrial.setArg (argIdx++, elementIndexB);
  //	clBasisTrial.setArg (argIdx++, *clTrialQuadPoints);
  //	clBasisTrial.setArg (argIdx++, trialPointCount);
  //	clBasisTrial.setArg (argIdx++, pointDim);
  //	clBasisTrial.setArg (argIdx++, trialDofCount);
  //	clBasisTrial.setArg (argIdx++, localDofIndexB);
  //	clBasisTrial.setArg (argIdx++, *clTrialValues);
  //	//        t0 = tbb::tick_count::now();
  //	m_openClHandler.enqueueKernel (cl::NDRange(trialPointCount));
  //	//	dt_kern += tbb::tick_count::now()-t0;
  //    }
  //    else
  //    {
  //      //        t0 = tbb::tick_count::now();
  //        clGlobalTrialPoints = m_openClHandler.createBuffer<CoordinateType> (
  //	    elementACount*trialPointCount*meshDim, CL_MEM_READ_WRITE);
  //        clGlobalTrialNormals = m_openClHandler.createBuffer<CoordinateType>
  // (
  //	    elementACount*trialPointCount*meshDim, CL_MEM_READ_WRITE);
  //	clTrialIntegrationElements =
  // m_openClHandler.createBuffer<CoordinateType> (
  //	    elementACount*trialPointCount, CL_MEM_READ_WRITE);
  //	//	dt_buf += tbb::tick_count::now()-t0;
  //	cl::Kernel &clMapTrial = m_openClHandler.setKernel
  //("clMapPointsAndNormalsToElements");
  //	argIdx = m_openClHandler.SetGeometryArgs (clMapTrial, 0);
  //	clMapTrial.setArg (argIdx++, *clTrialQuadPoints);
  //	clMapTrial.setArg (argIdx++, trialPointCount);
  //	clMapTrial.setArg (argIdx++, pointDim);
  //	clMapTrial.setArg (argIdx++, *clElementIndicesA);
  //	clMapTrial.setArg (argIdx++, elementACount);
  //	clMapTrial.setArg (argIdx++, *clGlobalTrialPoints);
  //	clMapTrial.setArg (argIdx++, *clGlobalTrialNormals);
  //	clMapTrial.setArg (argIdx++, *clTrialIntegrationElements);
  //	//        t0 = tbb::tick_count::now();
  //	m_openClHandler.enqueueKernel (cl::NDRange(elementACount,
  // trialPointCount));
  //	//	dt_kern += tbb::tick_count::now()-t0;

  //	//        t0 = tbb::tick_count::now();
  //	clGlobalTestPoints = m_openClHandler.createBuffer<CoordinateType> (
  //            testPointCount*meshDim, CL_MEM_READ_WRITE);
  //	clTestIntegrationElements = m_openClHandler.createBuffer<CoordinateType>
  //(
  //	    testPointCount, CL_MEM_READ_WRITE);
  //	//	dt_buf += tbb::tick_count::now()-t0;
  //	cl::Kernel &clMapTest = m_openClHandler.setKernel
  //("clMapPointsToElement");
  //	argIdx = m_openClHandler.SetGeometryArgs (clMapTest, 0);
  //	clMapTest.setArg (argIdx++, *clTestQuadPoints);
  //	clMapTest.setArg (argIdx++, testPointCount);
  //	clMapTest.setArg (argIdx++, pointDim);
  //	clMapTest.setArg (argIdx++, elementIndexB);
  //	clMapTest.setArg (argIdx++, *clGlobalTestPoints);
  //	clMapTest.setArg (argIdx++, *clTestIntegrationElements);
  //	//        t0 = tbb::tick_count::now();
  //	m_openClHandler.enqueueKernel (cl::NDRange(testPointCount));
  //	//	dt_kern += tbb::tick_count::now()-t0;

  //	//        t0 = tbb::tick_count::now();
  //	clTrialValues = m_openClHandler.createBuffer<BasisFunctionType> (
  //	    elementACount*trialPointCount*trialDofCount, CL_MEM_READ_WRITE);
  //	//	dt_buf += tbb::tick_count::now()-t0;
  //	cl::Kernel &clBasisTrial = m_openClHandler.setKernel
  //("clBasisfAElements");
  //	argIdx = m_openClHandler.SetGeometryArgs (clBasisTrial, 0);
  //	clBasisTrial.setArg (argIdx++, *clElementIndicesA);
  //	clBasisTrial.setArg (argIdx++, elementACount);
  //	clBasisTrial.setArg (argIdx++, *clTrialQuadPoints);
  //	clBasisTrial.setArg (argIdx++, trialPointCount);
  //	clBasisTrial.setArg (argIdx++, pointDim);
  //	clBasisTrial.setArg (argIdx++, trialDofCount);
  //	clBasisTrial.setArg (argIdx++, *clTrialValues);
  //	//        t0 = tbb::tick_count::now();
  //	m_openClHandler.enqueueKernel (cl::NDRange(elementACount,
  // trialPointCount));
  //	//	dt_kern += tbb::tick_count::now()-t0;

  //	//        t0 = tbb::tick_count::now();
  //	clTestValues = m_openClHandler.createBuffer<BasisFunctionType> (
  //	    testPointCount*testDofCount, CL_MEM_READ_WRITE);
  //	//	dt_buf += tbb::tick_count::now()-t0;
  //	cl::Kernel &clBasisTest = m_openClHandler.setKernel
  //("clBasisfBElement");
  //	argIdx = m_openClHandler.SetGeometryArgs (clBasisTest, 0);
  //	clBasisTest.setArg (argIdx++, elementIndexB);
  //	clBasisTest.setArg (argIdx++, *clTestQuadPoints);
  //	clBasisTest.setArg (argIdx++, testPointCount);
  //	clBasisTest.setArg (argIdx++, pointDim);
  //	clBasisTest.setArg (argIdx++, testDofCount);
  //	clBasisTest.setArg (argIdx++, localDofIndexB);
  //	clBasisTest.setArg (argIdx++, *clTestValues);
  //	//        t0 = tbb::tick_count::now();
  //	m_openClHandler.enqueueKernel (cl::NDRange(testPointCount));
  //	//	dt_kern += tbb::tick_count::now()-t0;
  //    }

  //    // Build the OpenCL kernel
  //    cl::Kernel &clKernel = m_openClHandler.setKernel ("clIntegrate");

  //    // Set kernel arguments
  //    argIdx = m_openClHandler.SetGeometryArgs (clKernel, 0);
  //    clKernel.setArg (argIdx++, *clGlobalTrialPoints);
  //    clKernel.setArg (argIdx++, *clGlobalTestPoints);
  //    clKernel.setArg (argIdx++, *clGlobalTrialNormals);
  //    clKernel.setArg (argIdx++, *clTrialIntegrationElements);
  //    clKernel.setArg (argIdx++, *clTestIntegrationElements);
  //    clKernel.setArg (argIdx++, *clTrialValues);
  //    clKernel.setArg (argIdx++, *clTestValues);
  //    clKernel.setArg (argIdx++, *clTrialQuadWeights);
  //    clKernel.setArg (argIdx++, *clTestQuadWeights);
  //    clKernel.setArg (argIdx++, trialPointCount);
  //    clKernel.setArg (argIdx++, testPointCount);
  //    clKernel.setArg (argIdx++, trialComponentCount);
  //    clKernel.setArg (argIdx++, testComponentCount);
  //    clKernel.setArg (argIdx++, trialDofCount);
  //    clKernel.setArg (argIdx++, testDofCount);
  //    clKernel.setArg (argIdx++, elementACount);
  //    clKernel.setArg (argIdx++, callVariant == TEST_TRIAL ? 1:0);
  //    clKernel.setArg (argIdx++, *clElementIndicesA);
  //    clKernel.setArg (argIdx++, elementIndexB);
  //    clKernel.setArg (argIdx++, *clResult);

  //    // Run the CL kernel
  //    //    t0 = tbb::tick_count::now();
  //    m_openClHandler.enqueueKernel (cl::NDRange(elementACount));
  //    //    dt_kern += tbb::tick_count::now()-t0;

  //    // Copy results back
  //    //    t0 = tbb::tick_count::now();
  //    m_openClHandler.pullCube<ResultType> (*clResult, result);
  //    //  dt_pull += tbb::tick_count::now()-t0;

  //    // Clean up local device buffers
  //    delete clElementIndicesA;
  //    delete clGlobalTrialPoints;
  //    delete clGlobalTestPoints;
  //    delete clGlobalTrialNormals;
  //    delete clTestValues;
  //    delete clTrialValues;
  //    delete clTrialIntegrationElements;
  //    delete clTestIntegrationElements;
  //    delete clResult;

  //    //  t_end = tbb::tick_count::now();
  //    //  static int callcount = 0;
  //    //  std::cout << callcount++ << '\t' << (t_end-t_start).seconds()
  //    //	      << "\tprog=" << dt_prog.seconds()
  //    //	      << "\tbuf=" << dt_buf.seconds()
  //    //	      << "\tkern=" << dt_kern.seconds()
  //    //	      << "\tpull=" << dt_pull.seconds()
  //    //	      << std::endl;

  //#else
  //    throw std::runtime_error ("Trying to call OpenCL method without OpenCL
  // support");
  //#endif // WITH_OPENCL
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void SeparableNumericalTestKernelTrialIntegrator<BasisFunctionType, KernelType,
                                                 ResultType, GeometryFactory>::
    integrate(const std::vector<ElementIndexPair> &elementIndexPairs,
              const Shapeset<BasisFunctionType> &testShapeset,
              const Shapeset<BasisFunctionType> &trialShapeset,
              const std::vector<Matrix<ResultType> *> &result) const {
  if (m_openClHandler.UseOpenCl()) {
    integrateCl(elementIndexPairs, testShapeset, trialShapeset, result);
  } else {
    integrateCpu(elementIndexPairs, testShapeset, trialShapeset, result);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void SeparableNumericalTestKernelTrialIntegrator<BasisFunctionType, KernelType,
                                                 ResultType, GeometryFactory>::
    integrateCpu(const std::vector<ElementIndexPair> &elementIndexPairs,
                 const Shapeset<BasisFunctionType> &testShapeset,
                 const Shapeset<BasisFunctionType> &trialShapeset,
                 const std::vector<Matrix<ResultType> *> &result) const {
  const int testPointCount = m_localTestQuadPoints.cols();
  const int trialPointCount = m_localTrialQuadPoints.cols();
  const int geometryPairCount = elementIndexPairs.size();

  if (result.size() != elementIndexPairs.size())
    throw std::invalid_argument(
        "NonseparableNumericalTestKernelTrialIntegrator::integrate(): "
        "arrays 'result' and 'elementIndicesA' must have the same number "
        "of elements");
  if (testPointCount == 0 || trialPointCount == 0 || geometryPairCount == 0)
    return;
  // TODO: in the (pathological) case that pointCount == 0 but
  // geometryPairCount != 0, set elements of result to 0.

  // Evaluate constants

  const int testDofCount = testShapeset.size();
  const int trialDofCount = trialShapeset.size();

  BasisData<BasisFunctionType> testBasisData, trialBasisData;
  GeometricalData<CoordinateType> *testGeomData = &m_testGeomData.local();
  GeometricalData<CoordinateType> *trialGeomData = &m_trialGeomData.local();
  const GeometricalData<CoordinateType> *constTestGeomData = testGeomData;
  const GeometricalData<CoordinateType> *constTrialGeomData = trialGeomData;

  size_t testBasisDeps = 0, trialBasisDeps = 0;
  size_t testGeomDeps = 0, trialGeomDeps = 0;

  m_testTransformations.addDependencies(testBasisDeps, testGeomDeps);
  m_trialTransformations.addDependencies(trialBasisDeps, trialGeomDeps);
  m_kernels.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
  m_integral.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

  typedef typename GeometryFactory::Geometry Geometry;
  std::unique_ptr<Geometry> testGeometry;
  std::unique_ptr<Geometry> trialGeometry;
  if (!m_cacheGeometricalData) {
    testGeometry = m_testGeometryFactory.make();
    trialGeometry = m_trialGeometryFactory.make();
  }

  CollectionOf3dArrays<BasisFunctionType> testValues, trialValues;
  CollectionOf4dArrays<KernelType> kernelValues;

  for (size_t i = 0; i < result.size(); ++i) {
    assert(result[i]);
    result[i]->resize(testDofCount, trialDofCount);
  }

  testShapeset.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS,
                        testBasisData);
  trialShapeset.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS,
                         trialBasisData);

  // Iterate over the elements
  for (int pairIndex = 0; pairIndex < geometryPairCount; ++pairIndex) {
    const int testElementIndex = elementIndexPairs[pairIndex].first;
    const int trialElementIndex = elementIndexPairs[pairIndex].second;
    if (m_cacheGeometricalData) {
      constTestGeomData = &m_cachedTestGeomData[testElementIndex];
      constTrialGeomData = &m_cachedTrialGeomData[trialElementIndex];
    } else {
      m_testRawGeometry.setupGeometry(testElementIndex, *testGeometry);
      m_trialRawGeometry.setupGeometry(trialElementIndex, *trialGeometry);
      testGeometry->getData(testGeomDeps, m_localTestQuadPoints, *testGeomData);
      if (testGeomDeps & DOMAIN_INDEX)
        testGeomData->domainIndex =
            m_testRawGeometry.domainIndex(testElementIndex);
      trialGeometry->getData(trialGeomDeps, m_localTrialQuadPoints,
                             *trialGeomData);
      if (trialGeomDeps & DOMAIN_INDEX)
        trialGeomData->domainIndex =
            m_trialRawGeometry.domainIndex(trialElementIndex);
    }
    m_testTransformations.evaluate(testBasisData, *constTestGeomData,
                                   testValues);
    m_trialTransformations.evaluate(trialBasisData, *constTrialGeomData,
                                    trialValues);

    m_kernels.evaluateOnGrid(*constTestGeomData, *constTrialGeomData,
                             kernelValues);
    m_integral.evaluateWithTensorQuadratureRule(
        *constTestGeomData, *constTrialGeomData, testValues, trialValues,
        kernelValues, m_testQuadWeights, m_trialQuadWeights,
        *result[pairIndex]);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void SeparableNumericalTestKernelTrialIntegrator<BasisFunctionType, KernelType,
                                                 ResultType, GeometryFactory>::
    integrateCl(const std::vector<ElementIndexPair> &elementIndexPairs,
                const Shapeset<BasisFunctionType> &testShapeset,
                const Shapeset<BasisFunctionType> &trialShapeset,
                const std::vector<Matrix<ResultType> *> &result) const {
  //#ifdef WITH_OPENCL
  //    const int testPointCount = m_localTestQuadPoints.cols();
  //    const int trialPointCount = m_localTrialQuadPoints.cols();
  //    const int geometryPairCount = elementIndexPairs.size();
  //    const int pointDim = m_localTestQuadPoints.rows();
  //    const int meshDim = m_openClHandler.meshGeom().size.dim;

  //    if (testPointCount == 0 || trialPointCount == 0 || geometryPairCount ==
  // 0)
  //        return;
  //    // TODO: in the (pathological) case that pointCount == 0 but
  //    // geometryPairCount != 0, set elements of result to 0.

  //    // TODO: add support for multiple-term expressions.
  //    if (!m_testTransformation.isTrivial() ||
  // !m_trialTransformation.isTrivial())
  //        throw
  // std::runtime_error("SeparableNumericalTestKernelTrialIntegrator::"
  //                                 "integrateCl(): multiple-term expression
  // lists are "
  //                                 "not supported in the OpenCL mode yet");
  //    const Expression<CoordinateType>& m_testExpression =
  // m_testTransformation.term(0);
  //    const Expression<CoordinateType>& m_trialExpression =
  // m_trialTransformation.term(0);

  //    // Evaluate constants
  //    const int testComponentCount = m_testExpression.codomainDimension();
  //    const int trialComponentCount = m_trialExpression.codomainDimension();
  //    const int testDofCount = testShapeset.size();
  //    const int trialDofCount = trialShapeset.size();

  //    const int kernelRowCount = m_kernels.codomainDimension();
  //    const int kernelColCount = m_kernels.domainDimension();

  //    // Assert that the kernel tensor dimensions are compatible
  //    // with the number of components of the functions

  //    const bool scalarKernel = (kernelRowCount == 1 && kernelColCount == 1);
  //    if (scalarKernel)
  //        assert(testComponentCount == trialComponentCount);
  //    else
  //    {
  //        assert(testComponentCount == kernelRowCount);
  //        assert(kernelColCount == trialComponentCount);
  //    }

  //    result.set_size(testDofCount, trialDofCount, geometryPairCount);

  //    int argIdx;

  //    // define device buffers
  //    cl::Buffer *clElementIndexA;
  //    cl::Buffer *clElementIndexB;
  //    cl::Buffer *clGlobalTrialPoints;
  //    cl::Buffer *clGlobalTestPoints;
  //    cl::Buffer *clGlobalTrialNormals;
  //    cl::Buffer *clTrialIntegrationElements;
  //    cl::Buffer *clTestIntegrationElements;
  //    cl::Buffer *clTrialValues;
  //    cl::Buffer *clTestValues;
  //    cl::Buffer *clResult;

  //    // Build the OpenCL program
  //    //std::vector<std::string> sources;
  //    cl::Program::Sources sources;
  //    sources.push_back (m_openClHandler.initStr());
  //    sources.push_back (testShapeset.clCodeString(true));
  //    sources.push_back (trialShapeset.clCodeString(false));
  //    sources.push_back (m_kernels.evaluateClCode());
  //    sources.push_back (clStrIntegrateRowOrCol());
  //    m_openClHandler.loadProgramFromStringArray (sources);

  //    // we need to separate the two index lists
  //    std::vector<int> elementIndexA(geometryPairCount);
  //    std::vector<int> elementIndexB(geometryPairCount);
  //    for (int i = 0; i < geometryPairCount; i++) {
  //        elementIndexA[i] = elementIndexPairs[i].first;
  //	elementIndexB[i] = elementIndexPairs[i].second;
  //    }

  //    // push the element index pairs
  //    clElementIndexA = m_openClHandler.pushVector<int> (elementIndexA);
  //    clElementIndexB = m_openClHandler.pushVector<int> (elementIndexB);
  //    clGlobalTestPoints = m_openClHandler.createBuffer<CoordinateType> (
  //        geometryPairCount*testPointCount*meshDim, CL_MEM_READ_WRITE);
  //    clGlobalTrialPoints = m_openClHandler.createBuffer<CoordinateType> (
  //        geometryPairCount*trialPointCount*meshDim, CL_MEM_READ_WRITE);
  //    clGlobalTrialNormals = m_openClHandler.createBuffer<CoordinateType> (
  //	geometryPairCount*trialPointCount*meshDim, CL_MEM_READ_WRITE);
  //    clTestIntegrationElements = m_openClHandler.createBuffer<CoordinateType>
  // (
  //        geometryPairCount*testPointCount, CL_MEM_READ_WRITE);
  //    clTrialIntegrationElements =
  // m_openClHandler.createBuffer<CoordinateType> (
  //        geometryPairCount*trialPointCount, CL_MEM_READ_WRITE);
  //    clTestValues = m_openClHandler.createBuffer<BasisFunctionType> (
  //	geometryPairCount*testPointCount*testDofCount, CL_MEM_READ_WRITE);
  //    clTrialValues = m_openClHandler.createBuffer<BasisFunctionType> (
  //	geometryPairCount*trialPointCount*trialDofCount, CL_MEM_READ_WRITE);
  //    clResult = m_openClHandler.createBuffer<ResultType>
  // (testDofCount*trialDofCount*geometryPairCount,
  //        CL_MEM_WRITE_ONLY);

  //    // Call the CL kernels to map the trial and test quadrature points
  //    cl::Kernel &clMapTest = m_openClHandler.setKernel
  // ("clMapPointsToElements");
  //    argIdx = m_openClHandler.SetGeometryArgs (clMapTest, 0);
  //    clMapTest.setArg (argIdx++, *clTestQuadPoints);
  //    clMapTest.setArg (argIdx++, testPointCount);
  //    clMapTest.setArg (argIdx++, pointDim);
  //    clMapTest.setArg (argIdx++, *clElementIndexA);
  //    clMapTest.setArg (argIdx++, geometryPairCount);
  //    clMapTest.setArg (argIdx++, *clGlobalTestPoints);
  //    clMapTest.setArg (argIdx++, *clTestIntegrationElements);
  //    m_openClHandler.enqueueKernel (cl::NDRange(geometryPairCount,
  // testPointCount));

  //    cl::Kernel &clMapTrial = m_openClHandler.setKernel
  // ("clMapPointsAndNormalsToElements");
  //    argIdx = m_openClHandler.SetGeometryArgs (clMapTrial, 0);
  //    clMapTrial.setArg (argIdx++, *clTrialQuadPoints);
  //    clMapTrial.setArg (argIdx++, trialPointCount);
  //    clMapTrial.setArg (argIdx++, pointDim);
  //    clMapTrial.setArg (argIdx++, *clElementIndexB);
  //    clMapTrial.setArg (argIdx++, *clGlobalTrialPoints);
  //    clMapTrial.setArg (argIdx++, *clGlobalTrialNormals);
  //    clMapTrial.setArg (argIdx++, *clTrialIntegrationElements);
  //    m_openClHandler.enqueueKernel (cl::NDRange(geometryPairCount,
  // trialPointCount));

  //    cl::Kernel &clBasisTest = m_openClHandler.setKernel
  // ("clBasisAElements");
  //    argIdx = m_openClHandler.SetGeometryArgs (clBasisTest, 0);
  //    clBasisTest.setArg (argIdx++, *clElementIndexA);
  //    clBasisTest.setArg (argIdx++, geometryPairCount);
  //    clBasisTest.setArg (argIdx++, *clTestQuadPoints);
  //    clBasisTest.setArg (argIdx++, testPointCount);
  //    clBasisTest.setArg (argIdx++, pointDim);
  //    clBasisTest.setArg (argIdx++, testDofCount);
  //    clBasisTest.setArg (argIdx++, *clTestValues);
  //    m_openClHandler.enqueueKernel (cl::NDRange(geometryPairCount,
  // testPointCount));

  //    cl::Kernel &clBasisTrial = m_openClHandler.setKernel
  // ("clBasisBElements");
  //    argIdx = m_openClHandler.SetGeometryArgs (clBasisTrial, 0);
  //    clBasisTest.setArg (argIdx++, *clElementIndexB);
  //    clBasisTest.setArg (argIdx++, geometryPairCount);
  //    clBasisTest.setArg (argIdx++, *clTrialQuadPoints);
  //    clBasisTest.setArg (argIdx++, trialPointCount);
  //    clBasisTest.setArg (argIdx++, pointDim);
  //    clBasisTest.setArg (argIdx++, trialDofCount);
  //    clBasisTest.setArg (argIdx++, *clTrialValues);
  //    m_openClHandler.enqueueKernel (cl::NDRange(geometryPairCount,
  // trialPointCount));

  //    // Build the OpenCL kernel
  //    cl::Kernel &clKernel = m_openClHandler.setKernel (scalarKernel ?
  //        "clIntegratePairsScalar" : "clIntegratePairs");

  //    // Set kernel arguments
  //    argIdx = m_openClHandler.SetGeometryArgs (clKernel, 0);
  //    clKernel.setArg (argIdx++, *clGlobalTrialPoints);
  //    clKernel.setArg (argIdx++, *clGlobalTestPoints);
  //    clKernel.setArg (argIdx++, *clGlobalTrialNormals);
  //    clKernel.setArg (argIdx++, *clTrialIntegrationElements);
  //    clKernel.setArg (argIdx++, *clTestIntegrationElements);
  //    clKernel.setArg (argIdx++, *clTrialValues);
  //    clKernel.setArg (argIdx++, *clTestValues);
  //    clKernel.setArg (argIdx++, *clTrialQuadWeights);
  //    clKernel.setArg (argIdx++, *clTestQuadWeights);
  //    clKernel.setArg (argIdx++, trialPointCount);
  //    clKernel.setArg (argIdx++, testPointCount);
  //    clKernel.setArg (argIdx++, trialComponentCount);
  //    clKernel.setArg (argIdx++, testComponentCount);
  //    clKernel.setArg (argIdx++, trialDofCount);
  //    clKernel.setArg (argIdx++, testDofCount);
  //    clKernel.setArg (argIdx++, geometryPairCount);
  //    clKernel.setArg (argIdx++, *clElementIndexA);
  //    clKernel.setArg (argIdx++, *clElementIndexB);
  //    clKernel.setArg (argIdx++, *clResult);

  //    // Run the CL kernel
  //    m_openClHandler.enqueueKernel (cl::NDRange(geometryPairCount));

  //    // Copy results back
  //    m_openClHandler.pullCube<ResultType> (*clResult, result);

  //    // clean up local device buffers
  //    delete clElementIndexA;
  //    delete clElementIndexB;
  //    delete clGlobalTestPoints;
  //    delete clGlobalTrialPoints;
  //    delete clGlobalTrialNormals;
  //    delete clTrialIntegrationElements;
  //    delete clTestIntegrationElements;
  //    delete clTestValues;
  //    delete clTrialValues;
  //    delete clResult;

  //#else
  //    throw std::runtime_error ("Trying to call OpenCL method without OpenCL
  // support");
  //#endif // WITH_OPENCL
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
const std::pair<const char *, int> SeparableNumericalTestKernelTrialIntegrator<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::clStrIntegrateRowOrCol() const {
  return std::make_pair(separable_numerical_double_integrator_cl,
                        separable_numerical_double_integrator_cl_len);
}

} // namespace Fiber
