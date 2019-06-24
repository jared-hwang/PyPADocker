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

#ifndef fiber_local_assembler_for_grid_functions_hpp
#define fiber_local_assembler_for_grid_functions_hpp

#include "../common/common.hpp"

#include "_2d_array.hpp"
#include "types.hpp"

#include <vector>
#include <stdexcept>

namespace Fiber {

/** \brief Local assembler for grid functions.

  This assembler provides methods that evaluate local (element-by-element) weak
  forms of integrals occurring in boundary-element matrices of grid functions.
 */
template <typename ResultType> class LocalAssemblerForGridFunctions {
public:
  virtual ~LocalAssemblerForGridFunctions() {}

  /** \brief Assemble local weak forms of a source term on specified elements.
   */
  virtual void
  evaluateLocalWeakForms(const std::vector<int> &elementIndices,
                         std::vector<Vector<ResultType>> &result) = 0;
};

} // namespace Fiber

#endif
