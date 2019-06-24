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

#include "space_helper.hpp"

#include "../common/acc.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../common/bounding_box_helpers.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "space.hpp"

namespace Bempp {

template <typename BasisFunctionType>
void SpaceHelper<BasisFunctionType>::
    getGlobalDofInterpolationPoints_defaultImplementation(
        const Space<BasisFunctionType> &space, Matrix<CoordinateType> &points) {
  std::vector<Point3D<CoordinateType>> vPoints;
  space.getGlobalDofPositions(vPoints);

  const size_t pointCount = vPoints.size();
  const size_t worldDim = space.worldDimension();
  points.resize(worldDim, pointCount);
  for (size_t p = 0; p < pointCount; ++p) {
    points(0, p) = acc(vPoints, p).x;
    if (worldDim > 1)
      points(1, p) = acc(vPoints, p).y;
    if (worldDim > 2)
      points(2, p) = acc(vPoints, p).z;
  }
}

template <typename BasisFunctionType>
void SpaceHelper<BasisFunctionType>::
    getNormalsAtGlobalDofInterpolationPoints_defaultImplementation(
        const Space<BasisFunctionType> &space,
        Matrix<CoordinateType> &normals) {
  std::vector<Point3D<CoordinateType>> vNormals;
  space.getGlobalDofNormals(vNormals);

  const size_t pointCount = vNormals.size();
  const size_t worldDim = space.worldDimension();
  normals.resize(worldDim, pointCount);
  for (size_t p = 0; p < pointCount; ++p) {
    normals(0, p) = acc(vNormals, p).x;
    if (worldDim > 1)
      normals(1, p) = acc(vNormals, p).y;
    if (worldDim > 2)
      normals(2, p) = acc(vNormals, p).z;
  }
}

template <typename BasisFunctionType>
void SpaceHelper<BasisFunctionType>::
    getGlobalDofBoundingBoxes_defaultImplementation(
        const GridView &view,
        const std::vector<std::vector<LocalDof>> &global2localDofs,
        std::vector<BoundingBox<CoordinateType>> &bboxes) {
  const IndexSet &indexSet = view.indexSet();
  const int elementCount = view.entityCount(0);

  std::vector<Matrix<CoordinateType>> elementCorners(elementCount);
  std::unique_ptr<EntityIterator<0>> it = view.entityIterator<0>();
  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    int index = indexSet.entityIndex(e);
    const Geometry &geo = e.geometry();
    geo.getCorners(acc(elementCorners, index));
    it->next();
  }

  BoundingBox<CoordinateType> model;
  const CoordinateType maxCoord = std::numeric_limits<CoordinateType>::max();
  model.lbound.x = model.lbound.y = model.lbound.z = maxCoord;
  model.ubound.x = model.ubound.y = model.ubound.z = -maxCoord;

  const int globalDofCount_ = global2localDofs.size();
  bboxes.resize(globalDofCount_, model);
  for (int i = 0; i < globalDofCount_; ++i) {
    const std::vector<LocalDof> &localDofs = acc(global2localDofs, i);
    BoundingBox<CoordinateType> &bbox = acc(bboxes, i);
    for (int j = 0; j < localDofs.size(); ++j)
      extendBoundingBox(bbox,
                        acc(elementCorners, acc(localDofs, j).entityIndex));
    assert(!localDofs.empty());
    setBoundingBoxReference<CoordinateType>(
        bbox, acc(elementCorners, localDofs[0].entityIndex)
                  .col(localDofs[0].dofIndex));
  }

#ifndef NDEBUG
  for (size_t i = 0; i < globalDofCount_; ++i) {
    assert(bboxes[i].reference.x >= bboxes[i].lbound.x);
    assert(bboxes[i].reference.y >= bboxes[i].lbound.y);
    assert(bboxes[i].reference.z >= bboxes[i].lbound.z);
    assert(bboxes[i].reference.x <= bboxes[i].ubound.x);
    assert(bboxes[i].reference.y <= bboxes[i].ubound.y);
    assert(bboxes[i].reference.z <= bboxes[i].ubound.z);
  }
#endif // NDEBUG
}

template <typename BasisFunctionType>
void SpaceHelper<BasisFunctionType>::getGlobalDofNormals_defaultImplementation(
    const GridView &view,
    const std::vector<std::vector<LocalDof>> &global2localDofs,
    std::vector<Point3D<CoordinateType>> &normals) {
  const int gridDim = view.dim();
  const int globalDofCount_ = global2localDofs.size();
  const int worldDim = view.dimWorld();
  normals.resize(globalDofCount_);

  const IndexSet &indexSet = view.indexSet();
  int elementCount = view.entityCount(0);

  Matrix<CoordinateType> elementNormals(worldDim, elementCount);
  std::unique_ptr<EntityIterator<0>> it = view.entityIterator<0>();
  Vector<CoordinateType> center(gridDim);
  // Note: we assume here that elements are flat and so the position at which
  // the normal is calculated does not matter.
  center.fill(0.5);
  Matrix<CoordinateType> normal;
  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    int index = indexSet.entityIndex(e);
    e.geometry().getNormals(center, normal);

    for (int dim = 0; dim < worldDim; ++dim)
      elementNormals(dim, index) = normal(dim, 0);
    it->next();
  }

  if (gridDim == 1)
    for (size_t g = 0; g < globalDofCount_; ++g) {
      const std::vector<LocalDof> &ldofs = acc(global2localDofs, g);
      normals[g].x = 0.;
      normals[g].y = 0.;
      for (size_t l = 0; l < ldofs.size(); ++l) {
        normals[g].x += elementNormals(0, acc(ldofs, l).entityIndex);
        normals[g].y += elementNormals(1, acc(ldofs, l).entityIndex);
      }
      double len =
          std::sqrt(normals[g].x * normals[g].x * normals[g].y * normals[g].y);
      normals[g].x /= len;
      normals[g].y /= len;
    }
  else // gridDim == 2
    for (size_t g = 0; g < globalDofCount_; ++g) {
      const std::vector<LocalDof> &ldofs = acc(global2localDofs, g);
      normals[g].x = 0.;
      normals[g].y = 0.;
      normals[g].z = 0.;
      for (size_t l = 0; l < ldofs.size(); ++l) {
        normals[g].x += elementNormals(0, acc(ldofs, l).entityIndex);
        normals[g].y += elementNormals(1, acc(ldofs, l).entityIndex);
        normals[g].z += elementNormals(2, acc(ldofs, l).entityIndex);
      }
      double len =
          std::sqrt(normals[g].x * normals[g].x + normals[g].y * normals[g].y +
                    normals[g].z * normals[g].z);
      normals[g].x /= len;
      normals[g].y /= len;
      normals[g].z /= len;
    }
}

template <typename BasisFunctionType>
void SpaceHelper<BasisFunctionType>::initializeLocal2FlatLocalDofMap(
    size_t flatLocalDofCount,
    const std::vector<std::vector<GlobalDofIndex>> &local2globalDofs,
    std::vector<LocalDof> &flatLocal2localDofs) {
  flatLocal2localDofs.clear();
  flatLocal2localDofs.reserve(flatLocalDofCount);
  for (size_t e = 0; e < local2globalDofs.size(); ++e)
    for (size_t dof = 0; dof < acc(local2globalDofs, e).size(); ++dof)
      if (acc(acc(local2globalDofs, e), dof) >= 0)
        flatLocal2localDofs.push_back(LocalDof(e, dof));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(SpaceHelper);

} // namespace Bempp
