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

#include "numerical_quadrature.hpp"

#include "bempp/common/config_data_types.hpp"

// Hyena code
#include "quadrature/galerkinduffy.hpp"
#include "quadrature/quadrature.hpp"

namespace Fiber {

// Helper functions in anonymous namespace
namespace {

// Regular integration

template <ELEMENT_SHAPE SHAPE, typename ValueType>
inline void reallyFillPointsAndWeightsRegular(int order,
                                              Matrix<ValueType> &points,
                                              std::vector<ValueType> &weights) {
  const int elementDim = 2;
  if (SHAPE == QUADRANGLE)
    // QuadratureRule<QUADRANGLE, GAUSS> expects to be given the number of
    // integration points in each dimension, unlike
    // QuadratureRule<TRIANGLE, GAUSS>, which expects the exactness degree.
    order = (order + 1 + 1) / 2;
  order = std::max(order, 1); // Hyena does not accept order == 0 for triangles
  const QuadratureRule<SHAPE, GAUSS> &rule(order);
  const int pointCount = rule.getNumPoints();
  points.resize(elementDim, pointCount);
  for (int i = 0; i < pointCount; ++i) {
    const Point2 point = rule.getPoint(i);
    for (int dim = 0; dim < elementDim; ++dim)
      points(dim, i) = point[dim];
    if (SHAPE == TRIANGLE)
      // Map points from
      // Hyena's reference triangle   (0,0)--(1,0)--(1,1)
      // to Dune's reference triangle (0,0)--(1,0)--(0,1).
      points(0, i) -= points(1, i);
  }
  weights.resize(pointCount);
  for (int i = 0; i < pointCount; ++i)
    weights[i] = rule.getWeight(i);
}

// Singular integration

template <typename ValueType>
inline void remapPointsSharedVertexTriangle(int sharedVertex,
                                            Matrix<ValueType> &points) {
  // Map vertex 0 to vertex #sharedVertex
  if (sharedVertex == 0)
    ; // do nothing
  else if (sharedVertex == 1) {
    const int pointCount = points.cols();
    Matrix<ValueType> oldPoints(points);
    for (int i = 0; i < pointCount; ++i) {
      points(0, i) = -oldPoints(0, i) - oldPoints(1, i) + 1.;
      points(1, i) = oldPoints(1, i);
    }
  } else if (sharedVertex == 2) {
    const int pointCount = points.cols();
    Matrix<ValueType> oldPoints(points);
    for (int i = 0; i < pointCount; ++i) {
      points(0, i) = oldPoints(0, i);
      points(1, i) = -oldPoints(0, i) - oldPoints(1, i) + 1.;
    }
  } else
    throw std::invalid_argument("remapPointsSharedVertexTriangle(): "
                                "sharedVertex must be 0, 1 or 2");
}

template <typename ValueType>
inline void remapPointsSharedVertexQuadrilateral(int sharedVertex,
                                                 Matrix<ValueType> &points) {
  throw std::runtime_error("remapPointsSharedVertexQuadrilateral(): "
                           "not implemented yet");
}

inline int nonsharedVertexTriangle(int sharedVertex0, int sharedVertex1) {
  //        if (sharedVertex0 + sharedVertex1 == 1)
  //            return 2;
  //        else if (sharedVertex0 + sharedVertex1 == 2)
  //            return 1;
  //        else if (sharedVertex0 + sharedVertex1 == 3)
  //            return 0;
  return 3 - sharedVertex0 - sharedVertex1;
}

template <typename ValueType>
inline void remapPointsSharedEdgeTriangle(int sharedVertex0, int sharedVertex1,
                                          Matrix<ValueType> &points) {
  assert(0 <= sharedVertex0 && sharedVertex0 <= 2);
  assert(0 <= sharedVertex1 && sharedVertex1 <= 2);
  assert(sharedVertex0 != sharedVertex1);
  assert(points.rows() == 2);

  // Vertices of the reference triangle
  Matrix<ValueType> refVertices(2, 3);
  refVertices.fill(0.);
  refVertices(0, 1) = 1.;
  refVertices(1, 2) = 1.;

  Matrix<ValueType> newVertices(2, 3);
  newVertices.col(0) = refVertices.col(sharedVertex0);
  newVertices.col(1) = refVertices.col(sharedVertex1);
  newVertices.col(2) =
      refVertices.col(nonsharedVertexTriangle(sharedVertex0, sharedVertex1));

  Vector<ValueType> b(newVertices.col(0));
  Matrix<ValueType> A(2, 2);
  A.col(0) = newVertices.col(1) - b;
  A.col(1) = newVertices.col(2) - b;

  Matrix<ValueType> oldPoints(points);

  // points := A * oldPoints + b[extended to pointCount columns)
  for (size_t col = 0; col < points.cols(); ++col)
    for (int dim = 0; dim < 2; ++dim)
      points(dim, col) = A(dim, 0) * oldPoints(0, col) +
                         A(dim, 1) * oldPoints(1, col) + b(dim);
}

template <typename ValueType>
inline void remapPointsSharedEdgeQuadrilateral(int sharedVertex0,
                                               int sharedVertex1,
                                               Matrix<ValueType> &points) {
  throw std::runtime_error("remapPointsSharedEdgeQuadrilateral(): "
                           "not implemented yet");
}

template <ELEMENT_SHAPE SHAPE, typename ValueType>
inline void remapPointsSharedVertex(int sharedVertex,
                                    Matrix<ValueType> &points) {
  // compile-time dispatch; circumvent the rule that no specialisation
  // of member template functions is possible
  BOOST_STATIC_ASSERT(SHAPE == TRIANGLE || SHAPE == QUADRANGLE);
  if (SHAPE == TRIANGLE)
    remapPointsSharedVertexTriangle(sharedVertex, points);
  else
    remapPointsSharedVertexQuadrilateral(sharedVertex, points);
}

template <ELEMENT_SHAPE SHAPE, typename ValueType>
inline void remapPointsSharedEdge(int sharedVertex0, int sharedVertex1,
                                  Matrix<ValueType> &points) {
  // compile-time dispatch; circumvent the rule that no specialisation
  // of member template functions is possible
  BOOST_STATIC_ASSERT(SHAPE == TRIANGLE || SHAPE == QUADRANGLE);
  if (SHAPE == TRIANGLE)
    remapPointsSharedEdgeTriangle<ValueType>(sharedVertex0, sharedVertex1,
                                             points);
  else
    remapPointsSharedEdgeQuadrilateral<ValueType>(sharedVertex0, sharedVertex1,
                                                  points);
}

template <ELEMENT_SHAPE SHAPE, SING_INT SINGULARITY, typename ValueType>
inline void reallyFillPointsAndWeightsSingular(
    const DoubleQuadratureDescriptor &desc, Matrix<ValueType> &testPoints,
    Matrix<ValueType> &trialPoints, std::vector<ValueType> &weights) {
  const int elementDim = 2;
  // Is there a more efficient way?
  const int order = std::max(desc.testOrder, desc.trialOrder);
  const int numPointsIn1d = (order + 1 + 1) / 2;
  // quadrangle regardless of SHAPE
  const QuadratureRule<QUADRANGLE, GAUSS> rule(numPointsIn1d);
  const GalerkinDuffyExpression<SHAPE, SINGULARITY> transform(rule);
  const int pointCount = rule.getNumPoints();
  const int regionCount = transform.getNumRegions();
  const int totalPointCount = regionCount * pointCount * pointCount;

  Point2 point;

  // Quadrature points
  testPoints.resize(elementDim, totalPointCount);
  trialPoints.resize(elementDim, totalPointCount);
  for (int region = 0; region < regionCount; ++region)
    for (int testIndex = 0; testIndex < pointCount; ++testIndex)
      for (int trialIndex = 0; trialIndex < pointCount; ++trialIndex) {
        int col = testIndex + trialIndex * pointCount +
                  region * pointCount * pointCount;

        // Test point
        point = transform.getPointX(testIndex, trialIndex, region);
        for (int dim = 0; dim < elementDim; ++dim)
          testPoints(dim, col) = point[dim];
        if (SHAPE == TRIANGLE)
          // Map points from
          // Hyena's reference triangle   (0,0)--(1,0)--(1,1)
          // to Dune's reference triangle (0,0)--(1,0)--(0,1).
          testPoints(0, col) -= testPoints(1, col);

        // Trial point
        point = transform.getPointY(testIndex, trialIndex, region);
        for (int dim = 0; dim < elementDim; ++dim)
          trialPoints(dim, col) = point[dim];
        if (SHAPE == TRIANGLE)
          // Map points from
          // Hyena's reference triangle   (0,0)--(1,0)--(1,1)
          // to Dune's reference triangle (0,0)--(1,0)--(0,1).
          trialPoints(0, col) -= trialPoints(1, col);
      }

  // Weights
  weights.resize(totalPointCount);
  for (int region = 0; region < regionCount; ++region)
    for (int testIndex = 0; testIndex < pointCount; ++testIndex)
      for (int trialIndex = 0; trialIndex < pointCount; ++trialIndex) {
        int col = testIndex + trialIndex * pointCount +
                  region * pointCount * pointCount;
        weights[col] = transform.getWeight(testIndex, trialIndex, region) *
                       rule.getWeight(testIndex) * rule.getWeight(trialIndex);
      }

  if (SINGULARITY == VRTX_ADJACENT) {
    remapPointsSharedVertex<SHAPE>(desc.topology.testSharedVertex0, testPoints);
    remapPointsSharedVertex<SHAPE>(desc.topology.trialSharedVertex0,
                                   trialPoints);
  } else if (SINGULARITY == EDGE_ADJACENT) {
    remapPointsSharedEdge<SHAPE>(desc.topology.testSharedVertex0,
                                 desc.topology.testSharedVertex1, testPoints);
    remapPointsSharedEdge<SHAPE>(desc.topology.trialSharedVertex0,
                                 desc.topology.trialSharedVertex1, trialPoints);
  }
}

} // namespace

// User-callable functions

template <typename ValueType>
void fillSingleQuadraturePointsAndWeights(int elementCornerCount,
                                          int accuracyOrder,
                                          Matrix<ValueType> &points,
                                          std::vector<ValueType> &weights) {
  if (elementCornerCount == 3)
    reallyFillPointsAndWeightsRegular<TRIANGLE>(accuracyOrder, points, weights);
  else if (elementCornerCount == 4)
    return reallyFillPointsAndWeightsRegular<QUADRANGLE>(accuracyOrder, points,
                                                         weights);
  else
    throw std::invalid_argument("fillSingleQuadraturePointsAndWeights(): "
                                "elementCornerCount must be either 3 or 4");
}

template <typename ValueType>
void fillDoubleSingularQuadraturePointsAndWeights(
    const DoubleQuadratureDescriptor &desc, Matrix<ValueType> &testPoints,
    Matrix<ValueType> &trialPoints, std::vector<ValueType> &weights) {
  const ElementPairTopology &topology = desc.topology;
  if (topology.testVertexCount == 3 && topology.trialVertexCount == 3) {
    if (topology.type == ElementPairTopology::SharedVertex)
      reallyFillPointsAndWeightsSingular<TRIANGLE, VRTX_ADJACENT>(
          desc, testPoints, trialPoints, weights);
    else if (topology.type == ElementPairTopology::SharedEdge)
      reallyFillPointsAndWeightsSingular<TRIANGLE, EDGE_ADJACENT>(
          desc, testPoints, trialPoints, weights);
    else if (topology.type == ElementPairTopology::Coincident)
      reallyFillPointsAndWeightsSingular<TRIANGLE, COINCIDENT>(
          desc, testPoints, trialPoints, weights);
    else
      throw std::invalid_argument(
          "fillDoubleSingularQuadraturePointsAndWeights(): "
          "Invalid element configuration");
  } else if (topology.testVertexCount == 4 && topology.trialVertexCount == 4) {
    if (topology.type == ElementPairTopology::SharedVertex)
      reallyFillPointsAndWeightsSingular<QUADRANGLE, VRTX_ADJACENT>(
          desc, testPoints, trialPoints, weights);
    else if (topology.type == ElementPairTopology::SharedEdge)
      reallyFillPointsAndWeightsSingular<QUADRANGLE, EDGE_ADJACENT>(
          desc, testPoints, trialPoints, weights);
    else if (topology.type == ElementPairTopology::Coincident)
      reallyFillPointsAndWeightsSingular<QUADRANGLE, COINCIDENT>(
          desc, testPoints, trialPoints, weights);
    else
      throw std::invalid_argument(
          "fillDoubleSingularQuadraturePointsAndWeights(): "
          "Invalid element configuration");
  } else
    throw std::invalid_argument(
        "fillDoubleSingularQuadraturePointsAndWeights(): "
        "Singular quadrature rules for mixed "
        "meshes are not implemented yet.");
}

#ifdef ENABLE_SINGLE_PRECISION
template void fillSingleQuadraturePointsAndWeights<float>(
    int elementCornerCount, int accuracyOrder, Matrix<float> &points,
    std::vector<float> &weights);
template void fillDoubleSingularQuadraturePointsAndWeights<float>(
    const DoubleQuadratureDescriptor &desc, Matrix<float> &testPoints,
    Matrix<float> &trialPoints, std::vector<float> &weights);
#endif
#ifdef ENABLE_DOUBLE_PRECISION
template void fillSingleQuadraturePointsAndWeights<double>(
    int elementCornerCount, int accuracyOrder, Matrix<double> &points,
    std::vector<double> &weights);
template void fillDoubleSingularQuadraturePointsAndWeights<double>(
    const DoubleQuadratureDescriptor &desc, Matrix<double> &testPoints,
    Matrix<double> &trialPoints, std::vector<double> &weights);
#endif

} // namespace Fiber
