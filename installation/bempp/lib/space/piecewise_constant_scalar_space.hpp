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

#ifndef bempp_piecewise_constant_scalar_space_hpp
#define bempp_piecewise_constant_scalar_space_hpp

#include "../common/common.hpp"

#include "../grid/grid_view.hpp"
#include "../grid/grid_segment.hpp"
#include "scalar_space.hpp"
#include "../common/types.hpp"
#include "../fiber/constant_scalar_shapeset.hpp"

#include <map>
#include <memory>
#include <tbb/mutex.h>

namespace Bempp {

/** \cond FORWARD_DECL */
class GridView;
/** \endcond */

/** \ingroup space
 *  \brief Space of piecewise constant scalar functions. */
template <typename BasisFunctionType>
class PiecewiseConstantScalarSpace : public ScalarSpace<BasisFunctionType> {
public:
  typedef
      typename ScalarSpace<BasisFunctionType>::CoordinateType CoordinateType;

  /** \brief Constructor.
   *
   *  Construct a space of piecewise constant scalar functions
   *  defined on the grid \p grid.
   *
   *  An exception is thrown if \p grid is a null pointer.
   */
  explicit PiecewiseConstantScalarSpace(const shared_ptr<const Grid> &grid);

  /** \brief Constructor.
   *
   *  Construct a space of piecewise constant scalar functions defined on the
   *  elements of the grid \p grid belonging to the segment \p segment.
   *
   *  An exception is thrown if \p grid is a null pointer.
   */
  PiecewiseConstantScalarSpace(const shared_ptr<const Grid> &grid,
                               const GridSegment &segment);

  virtual shared_ptr<const Space<BasisFunctionType>> discontinuousSpace(
      const shared_ptr<const Space<BasisFunctionType>> &self) const;
  virtual bool isDiscontinuous() const;

  virtual bool isBarycentric() const { return false; }

  virtual shared_ptr<const Space<BasisFunctionType>> barycentricSpace(
      const shared_ptr<const Space<BasisFunctionType>> &self) const;

  virtual int domainDimension() const;
  virtual int codomainDimension() const;

  /** \brief Return the variant of element \p element.
   *
   *  Possible return values:
   *    - 2: one-dimensional segment,
   *    - 3: triangular element,
   *    - 4: quadrilateral element. */
  virtual ElementVariant elementVariant(const Entity<0> &element) const;
  virtual void setElementVariant(const Entity<0> &element,
                                 ElementVariant variant);

  virtual const Fiber::Shapeset<BasisFunctionType> &
  shapeset(const Entity<0> &element) const;

  virtual bool spaceIsCompatible(const Space<BasisFunctionType> &other) const;

  virtual SpaceIdentifier spaceIdentifier() const {
    return PIECEWISE_CONSTANT_SCALAR;
  }

  virtual size_t globalDofCount() const;
  virtual size_t flatLocalDofCount() const;
  virtual void getGlobalDofs(const Entity<0> &element,
                             std::vector<GlobalDofIndex> &dofs) const;
  virtual void
  global2localDofs(const std::vector<GlobalDofIndex> &globalDofs,
                   std::vector<std::vector<LocalDof>> &localDofs) const;
  virtual void
  flatLocal2localDofs(const std::vector<FlatLocalDofIndex> &globalDofs,
                      std::vector<LocalDof> &localDofs) const;

  virtual void
  getGlobalDofInterpolationPoints(Matrix<CoordinateType> &points) const;
  virtual void getNormalsAtGlobalDofInterpolationPoints(
      Matrix<CoordinateType> &normals) const;

  virtual void
  getGlobalDofPositions(std::vector<Point3D<CoordinateType>> &positions) const;
  virtual void getFlatLocalDofPositions(
      std::vector<Point3D<CoordinateType>> &positions) const;

  virtual void getGlobalDofBoundingBoxes(
      std::vector<BoundingBox<CoordinateType>> &bboxes) const;
  virtual void getFlatLocalDofBoundingBoxes(
      std::vector<BoundingBox<CoordinateType>> &bboxes) const;

  virtual void
  getGlobalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const;
  virtual void
  getFlatLocalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const;

  virtual void
  dumpClusterIds(const char *fileName,
                 const std::vector<unsigned int> &clusterIdsOfGlobalDofs) const;
  virtual void
  dumpClusterIdsEx(const char *fileName,
                   const std::vector<unsigned int> &clusterIdsOfGlobalDofs,
                   DofType dofType) const;

private:
  void assignDofsImpl(const GridSegment &segment);

private:
  std::unique_ptr<GridView> m_view;
  GridSegment m_segment;
  Fiber::ConstantScalarShapeset<BasisFunctionType> m_shapeset;
  std::vector<std::vector<GlobalDofIndex>> m_local2globalDofs;
  std::vector<std::vector<LocalDof>> m_global2localDofs;
  mutable shared_ptr<Space<BasisFunctionType>> m_barycentricSpace;
  mutable tbb::mutex m_barycentricSpaceMutex;
};

/** \brief Define a PiecewiseConstantScalarSpace that has an update method for
 * grid refinement. */
template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptivePiecewiseConstantScalarSpace(const shared_ptr<const Grid> &grid);

/** \brief Overload to define a set of domains for the space and whether the
 space contains boundary entities
 (\p open = true) or not. */
template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptivePiecewiseConstantScalarSpace(const shared_ptr<const Grid> &grid,
                                     const std::vector<int> &domains,
                                     bool open);

/** \brief Overlad. */
template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptivePiecewiseConstantScalarSpace(const shared_ptr<const Grid> &grid,
                                     int domain, bool open);

} // namespace Bempp

#endif
