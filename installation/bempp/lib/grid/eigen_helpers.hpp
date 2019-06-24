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
#ifndef bempp_eigen_helpers_hpp
#define bempp_eigen_helpers_hpp

#include "../common/common.hpp"
#include "../common/eigen_support.hpp"

#include <dune/common/fvector.hh>

using namespace Bempp;

// Internal implementations for general Eigen objects
template <typename M, typename T, int size>
bool _eigen_fieldvector_compare(const M &x,
                                const Dune::FieldVector<T, size> &y) {
  if (x.rows() != size)
    return false;
  for (int i = 0; i < size; ++i)
    if (x(i) != y[i])
      return false;
  return true;
}

template <typename M, typename T, int rows, int cols>
bool _eigen_fieldmatrix_compare(const M &x,
                                const Dune::FieldMatrix<T, rows, cols> &y) {
  if (x.rows() != rows || x.cols() != cols)
    return false;
  for (int j = 0; j < cols; ++j)
    for (int i = 0; i < rows; ++i)
      if (x(i, j) != y[i][j])
        return false;
  return true;
}

template <typename T, int size>
bool operator==(const Vector<T> &x, const Dune::FieldVector<T, size> &y) {
  return _eigen_fieldvector_compare(x, y);
}

template <typename T, int size>
bool operator==(const Dune::FieldVector<T, size> &x, const Vector<T> &y) {
  return _eigen_fieldvector_compare(y, x);
}

template <typename T, int rows, int cols>
bool operator==(const Bempp::Matrix<T> &x,
                const Dune::FieldMatrix<T, rows, cols> &y) {
  return _eigen_fieldmatrix_compare(x, y);
}

template <typename T, int rows, int cols>
bool operator==(const Dune::FieldMatrix<T, rows, cols> &x,
                const Bempp::Matrix<T> &y) {
  return _eigen_fieldmatrix_compare(y, x);
}

// For BOOST_CHECK_EQUAL to work, we need to copy these comparison operators
// to the boost::test_tools:tt_detail namespace...

namespace boost {
namespace test_tools {
namespace tt_detail {

template <typename T, int size>
bool operator==(const Bempp::Vector<T> &x,
                const Dune::FieldVector<T, size> &y) {
  return ::operator==(x, y);
}

template <typename T, int size>
bool operator==(const Dune::FieldVector<T, size> &x,
                const Bempp::Vector<T> &y) {
  return ::operator==(x, y);
}

template <typename T, int rows, int cols>
bool operator==(const Bempp::Matrix<T> &x,
                const Dune::FieldMatrix<T, rows, cols> &y) {
  return ::operator==(x, y);
}

template <typename T, int rows, int cols>
bool operator==(const Dune::FieldMatrix<T, rows, cols> &x,
                const Bempp::Matrix<T> &y) {
  return ::operator==(x, y);
}

} // namespace tt_detail
} // namespace test_tools
} // namespace boost

#endif
