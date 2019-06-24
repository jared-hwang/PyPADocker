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

#ifndef bempp_elementary_integral_operator_base_hpp
#define bempp_elementary_integral_operator_base_hpp

#include "../common/common.hpp"
#include "../common/types.hpp"

#include "abstract_boundary_operator.hpp"

#include "../common/shared_ptr.hpp"

#include <stdexcept>
#include <vector>

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForIntegralOperators;
template <typename CoordinateType> class RawGridGeometry;
template <typename ValueType> class Basis;
class OpenClHandler;
/** \endcond */

} // namespace Fiber

namespace Bempp {

/** \ingroup abstract_boundary_operators
 *  \brief Base class of ElementaryIntegralOperator, containing functionality
 *  independent from \c KernelType.
 */
template <typename BasisFunctionType_, typename ResultType_>
class ElementaryIntegralOperatorBase
    : public AbstractBoundaryOperator<BasisFunctionType_, ResultType_> {
  typedef AbstractBoundaryOperator<BasisFunctionType_, ResultType_> Base;

public:
  /** \copydoc AbstractBoundaryOperator::BasisFunctionType */
  typedef typename Base::BasisFunctionType BasisFunctionType;
  /** \copydoc AbstractBoundaryOperator::ResultType */
  typedef typename Base::ResultType ResultType;
  /** \copydoc AbstractBoundaryOperator::CoordinateType */
  typedef typename Base::CoordinateType CoordinateType;
  /** \copydoc AbstractBoundaryOperator::QuadratureStrategy */
  typedef typename Base::QuadratureStrategy QuadratureStrategy;
  /** \brief Type of the appropriate instantiation of
   * Fiber::LocalAssemblerForOperators. */
  typedef Fiber::LocalAssemblerForIntegralOperators<ResultType> LocalAssembler;

  /** \copydoc AbstractBoundaryOperator::AbstractBoundaryOperator */
  ElementaryIntegralOperatorBase(
      const shared_ptr<const Space<BasisFunctionType>> &domain,
      const shared_ptr<const Space<BasisFunctionType>> &range,
      const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
      const std::string &label, int symmetry);

  /** \brief Destructor. */
  ~ElementaryIntegralOperatorBase();

  /** \brief Construct a local assembler suitable for this operator.
   *
   *  \param[in] quadStrategy  Quadrature strategy to be used to construct the
   *assembler.
   *
   *  (TODO: finish description of the other parameters.)
   */
  std::unique_ptr<LocalAssembler> makeAssembler(
      const QuadratureStrategy &quadStrategy,
      const shared_ptr<const GeometryFactory> &testGeometryFactory,
      const shared_ptr<const GeometryFactory> &trialGeometryFactory,
      const shared_ptr<const Fiber::RawGridGeometry<CoordinateType>>
          &testRawGeometry,
      const shared_ptr<const Fiber::RawGridGeometry<CoordinateType>>
          &trialRawGeometry,
      const shared_ptr<const std::vector<
          const Fiber::Shapeset<BasisFunctionType> *>> &testShapesets,
      const shared_ptr<const std::vector<
          const Fiber::Shapeset<BasisFunctionType> *>> &trialShapesets,
      const shared_ptr<const Fiber::OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions,
      VerbosityLevel::Level verbosityLevel, bool cacheSingularIntegrals) const;

  /** \brief Construct a local assembler suitable for this operator using a
   *  specified quadrature strategy.
   *
   *  \param[in] quadStrategy  Quadrature strategy to be used to construct the
   *assembler.
   *  \param[in] options           Assembly options.
   *
   *  This is an overloaded function, provided for convenience. It
   *  automatically constructs most of the arguments required by the other
   *  overload.
   */
  std::unique_ptr<LocalAssembler>
  makeAssembler(const QuadratureStrategy &quadStrategy,
                const AssemblyOptions &options) const;

  /** \brief Overload that takes a parameter list instead of a context object.
   */
  std::unique_ptr<LocalAssembler>
  makeAssembler(const ParameterList &parameterList) const;

  /** \brief Assemble the operator's weak form using a specified local
   *assembler.
   *
   *  This function is intended for internal use of the library. End users
   *  should not need to call it directly. They should use
   *  AbstractBoundaryOperator::assembleWeakForm() instead. */
  shared_ptr<DiscreteBoundaryOperator<ResultType_>> assembleWeakFormInternal(
      LocalAssembler &assembler,
      const Context<BasisFunctionType, ResultType> &context) const;

private:
  /** \brief Construct a local assembler suitable for this operator.
   *
   *  This function is invoked by both overloads of makeAssembler()
   *  to do the actual work. */
  virtual std::unique_ptr<LocalAssembler> makeAssemblerImpl(
      const QuadratureStrategy &quadStrategy,
      const shared_ptr<const GeometryFactory> &testGeometryFactory,
      const shared_ptr<const GeometryFactory> &trialGeometryFactory,
      const shared_ptr<const Fiber::RawGridGeometry<CoordinateType>>
          &testRawGeometry,
      const shared_ptr<const Fiber::RawGridGeometry<CoordinateType>>
          &trialRawGeometry,
      const shared_ptr<const std::vector<
          const Fiber::Shapeset<BasisFunctionType> *>> &testShapesets,
      const shared_ptr<const std::vector<
          const Fiber::Shapeset<BasisFunctionType> *>> &trialShapesets,
      const shared_ptr<const Fiber::OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions,
      VerbosityLevel::Level verbosityLevel,
      bool cacheSingularIntegrals) const = 0;

  /** \brief Assemble the operator's weak form using a specified local
   *assembler.
   *
   *  This virtual function is invoked by assembleWeakFormInternal()
   *  to do the actual work. */
  virtual shared_ptr<DiscreteBoundaryOperator<ResultType_>>
  assembleWeakFormInternalImpl2(
      LocalAssembler &assembler,
      const Context<BasisFunctionType_, ResultType_> &options) const = 0;
};

} // namespace Bempp

#endif
