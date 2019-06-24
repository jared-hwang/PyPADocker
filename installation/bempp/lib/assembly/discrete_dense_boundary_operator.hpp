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

#ifndef bempp_discrete_dense_boundary_operator_hpp
#define bempp_discrete_dense_boundary_operator_hpp

#include "../common/common.hpp"
#include "../common/eigen_support.hpp"

#include "discrete_boundary_operator.hpp"

#include "../common/shared_ptr.hpp"

namespace Bempp {

/** \ingroup discrete_boundary_operators
 *  \brief Discrete boundary operator stored as a dense matrix. */
template <typename ValueType>
class DiscreteDenseBoundaryOperator
    : public DiscreteBoundaryOperator<ValueType> {
public:
  /** \brief Constructor.
   *
   *  Construct a discrete boundary operator represented by the matrix \p mat.
   */
  explicit DiscreteDenseBoundaryOperator(const Matrix<ValueType> &mat);

  virtual void dump() const;

  virtual Matrix<ValueType> asMatrix() const;

  virtual unsigned int rowCount() const;
  virtual unsigned int columnCount() const;

  virtual void addBlock(const std::vector<int> &rows,
                        const std::vector<int> &cols, const ValueType alpha,
                        Matrix<ValueType> &block) const;

  const Matrix<ValueType> &data() const;

private:
  virtual void applyBuiltInImpl(const TranspositionMode trans,
                                const Eigen::Ref<Vector<ValueType>> &x_in,
                                Eigen::Ref<Vector<ValueType>> y_inout,
                                const ValueType alpha,
                                const ValueType beta) const;

private:
  /** \cond PRIVATE */
  Matrix<ValueType> m_mat;
  /** \endcond */
};

/** \brief
 *  Return a shared pointer to newly constructed DiscreteDenseBoundaryOperator
 *  wrapping a specified matrix. */
template <typename ValueType>
shared_ptr<DiscreteDenseBoundaryOperator<ValueType>>
discreteDenseBoundaryOperator(const Matrix<ValueType> &mat);

} // namespace Bempp

#endif
