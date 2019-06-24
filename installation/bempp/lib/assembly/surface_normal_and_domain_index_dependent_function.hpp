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

#ifndef bempp_surface_normal_and_domain_index_dependent_function_hpp
#define bempp_surface_normal_and_domain_index_dependent_function_hpp

#include "../common/common.hpp"
#include "../common/eigen_support.hpp"

#include "../fiber/surface_normal_and_domain_index_dependent_function.hpp"

namespace Bempp {
using Fiber::SurfaceNormalAndDomainIndexDependentFunction;

/** \ingroup assembly_functions

  \brief Use a functor taking a point's global coordinates, the local
  surface-normal unit vector and the domain index of the element containing the
  point to construct an instance of the Function class.

  The template parameter \p Functor should be a class implementing the following
  interface:

  \code
  class Functor
  {
  public:
      // Type of the function's values (e.g. float or std::complex<double>)
      typedef <implementiation-defined> ValueType;
      typedef ScalarTraits<ValueType>::RealType CoordinateType;

      // Number of components of the function's arguments ("point" and "normal")
      int argumentDimension() const;

      // Number of components of the function's result
      int resultDimension() const;

      // Evaluate the function at the point "point" lying on an element from
      // domain "domainIndex", with "normal" the local unit normal vector, and
      // store result in the array "result". All arrays will be preinitialised
      // to correct dimensions.
      void evaluate(const Vector<CoordinateType>& point,
                    const Vector<CoordinateType>& normal,
                    int domainIndex,
                    Vector<ValueType>& result) const;
  };
  \endcode

  The constructed Function object can subsequently be passed into a constructor
  of the GridFunction class. */
template <typename Functor>
inline SurfaceNormalAndDomainIndexDependentFunction<Functor>
surfaceNormalAndDomainIndexDependentFunction(const Functor &functor) {
  return SurfaceNormalAndDomainIndexDependentFunction<Functor>(functor);
}

} // namespace Bempp

#endif
