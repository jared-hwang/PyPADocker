// Copyright (C) 2011-2013 by the BEM++ Authors
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

#include "assembled_potential_operator.hpp"

#include "discrete_boundary_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../space/space.hpp"
#include "../common/eigen_support.hpp"

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
AssembledPotentialOperator<BasisFunctionType, ResultType>::
    AssembledPotentialOperator(
        const shared_ptr<const Space<BasisFunctionType>> &space_,
        const shared_ptr<const Matrix<CoordinateType>> &evaluationPoints_,
        const shared_ptr<const DiscreteBoundaryOperator<ResultType>> &op_,
        int componentCount_)
    : m_space(space_), m_evaluationPoints(evaluationPoints_), m_op(op_),
      m_componentCount(componentCount_) {
  if (!m_op)
    throw std::invalid_argument("AssembledPotentialOperator::"
                                "AssembledPotentialOperator(): "
                                "op_ must not be null");
  if (m_componentCount < 1)
    throw std::invalid_argument("AssembledPotentialOperator::"
                                "AssembledPotentialOperator(): "
                                "componentCount must be positive");
  if (m_op->rowCount() % m_componentCount != 0)
    throw std::invalid_argument("AssembledPotentialOperator::"
                                "AssembledPotentialOperator(): "
                                "incorrect number of rows in op_");
  if (m_evaluationPoints &&
      m_op->rowCount() != m_componentCount * m_evaluationPoints->cols())
    throw std::invalid_argument("AssembledPotentialOperator::"
                                "AssembledPotentialOperator(): "
                                "incorrect number of rows in op_");
  if (m_space && m_op->columnCount() != m_space->globalDofCount())
    throw std::invalid_argument(
        "AssembledPotentialOperator::AssembledPotentialOperator(): "
        "number of columns of op_ must match the number of global DOFs "
        "of space_");
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Space<BasisFunctionType>>
AssembledPotentialOperator<BasisFunctionType, ResultType>::space() const {
  return m_space;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Matrix<typename AssembledPotentialOperator<
    BasisFunctionType, ResultType>::CoordinateType>>
AssembledPotentialOperator<BasisFunctionType, ResultType>::evaluationPoints()
    const {
  return m_evaluationPoints;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>>
AssembledPotentialOperator<BasisFunctionType, ResultType>::discreteOperator()
    const {
  return m_op;
}

template <typename BasisFunctionType, typename ResultType>
int AssembledPotentialOperator<BasisFunctionType, ResultType>::componentCount()
    const {
  return m_componentCount;
}

template <typename BasisFunctionType, typename ResultType>
Matrix<ResultType>
AssembledPotentialOperator<BasisFunctionType, ResultType>::apply(
    const Vector<ResultType> &coefficients) const {
  if (m_space->globalDofCount() != coefficients.rows())
    throw std::invalid_argument(
        "AssembledPotentialOperator::apply(): "
        "space used to expand 'argument' does not "
        "match the one used during operator construction");
  Matrix<ResultType> result(m_op->rowCount(), 1);
  Matrix<ResultType> coeffsMatrix(coefficients.rows(), 1);
  coeffsMatrix.col(0) = coefficients;
  m_op->apply(NO_TRANSPOSE, coeffsMatrix, result, 1., 0.);
  assert(result.rows() % m_componentCount == 0);
  result.resize(m_componentCount, result.rows() / m_componentCount);
  return result;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    AssembledPotentialOperator);

} // namespace Bempp
