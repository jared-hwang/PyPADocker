#ifndef bempp_grid_triangle_imp_data_container_impl_hpp
#define bempp_grid_triangle_imp_data_container_impl_hpp

#include "bempp_grid_data_container.hpp"
#include <dune/geometry/type.hh>
#include <cassert>
#include <utility>

namespace BemppGrid {

inline DataContainer::DataContainer() : m_levels(0), m_idCounter(0){};

inline DataContainer::~DataContainer(){};

inline void DataContainer::init(
    const shared_ptr<DataContainer::NodesContainer> &nodes,
    const shared_ptr<DataContainer::ElementsContainer> &elements) {

  int nelements = elements->size();
  int nnodes = nodes->size();

  m_nodes.push_back(nodes);
  m_elements.push_back(elements);
  m_edges.push_back(EdgesContainer());
  auto &edges = m_edges[m_levels];
  m_element2Edges.push_back(std::vector<std::array<unsigned int, 3>>());
  auto &element2Edges = m_element2Edges[m_levels];

  element2Edges.resize(nelements);
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>>
      nodes2EdgeIndexPair;
  nodes2EdgeIndexPair.resize(nodes->size());

  // Numbering of edges in Dune is not standard.
  // localEdgeNodes[i] has the reference element node numbers of the ith edge.
  unsigned int localEdgeNodes[3][2]{{0, 1}, {2, 0}, {1, 2}};

  for (int elementIndex = 0; elementIndex < nelements; ++elementIndex) {
    const auto &element = (*elements)[elementIndex];
    for (int i = 0; i < 3; ++i) {
      unsigned int n0 = element[localEdgeNodes[i][0]];
      unsigned int n1 = element[localEdgeNodes[i][1]];
      if (n1 > n0)
        std::swap(n0, n1); // Number edges from smaller to larger vertex
      // Check if edge already exists
      bool edgeExists = false;
      unsigned int edgeIndex;
      for (const auto &indexPair : nodes2EdgeIndexPair[n0])
        if (n1 == indexPair.first) {
          edgeExists = true;
          edgeIndex = indexPair.second;
          break;
        }
      if (!edgeExists) {
        // Create Edge
        edgeIndex = edges.size();
        nodes2EdgeIndexPair[n0].push_back(
            std::pair<unsigned int, unsigned int>(n1, edgeIndex));
        edges.push_back(std::array<unsigned int, 2>({n0, n1}));
      }
      element2Edges[elementIndex][i] = edgeIndex;
    }
  }

  // Fill the connectivity data arrays
  m_edge2Elements.push_back(std::vector<std::vector<unsigned int>>());
  m_node2Elements.push_back(std::vector<std::vector<unsigned int>>());
  m_node2Edges.push_back(std::vector<std::vector<unsigned int>>());

  auto &edge2Elements = m_edge2Elements[m_levels];
  auto &node2Elements = m_node2Elements[m_levels];
  auto &node2Edges = m_node2Edges[m_levels];

  edge2Elements.resize(edges.size());
  node2Elements.resize(nnodes);
  node2Edges.resize(nnodes);

  for (unsigned int i = 0; i < nelements; ++i) {
    for (int j = 0; j < 3; ++j) {
      edge2Elements[element2Edges[i][j]].push_back(i);
    }
  }

  for (unsigned int i = 0; i < nelements; ++i)
    for (int j = 0; j < 3; ++j)
      node2Elements[(*elements)[i][j]].push_back(i);

  for (unsigned int i = 0; i < edges.size(); ++i)
    for (int j = 0; j < 2; ++j)
      node2Edges[edges[i][j]].push_back(i);

  // Now generate the Ids

  m_nodeIds.push_back(std::vector<unsigned int>(m_nodes[m_levels]->size()));
  m_elementIds.push_back(
      std::vector<unsigned int>(m_elements[m_levels]->size()));
  m_edgeIds.push_back(std::vector<unsigned int>(m_edges[m_levels].size()));

  for (unsigned int i = 0; i < m_nodes[m_levels]->size(); ++i)
    m_nodeIds[m_levels][i] = m_idCounter++;

  for (unsigned int i = 0; i < m_elements[m_levels]->size(); ++i)
    m_elementIds[m_levels][i] = m_idCounter++;

  for (unsigned int i = 0; i < m_edges[m_levels].size(); ++i)
    m_edgeIds[m_levels][i] = m_idCounter++;

  m_levels++;

  // Compute Geometries

  computeGeometries<0>(0);
  computeGeometries<1>(0);
  computeGeometries<2>(0);
}

inline int DataContainer::levels() const { return m_levels; }

inline const DataContainer::NodesContainer &
DataContainer::nodes(int level) const {
  assert(level < m_levels);
  return *(m_nodes[level]);
}

inline const DataContainer::ElementsContainer &
DataContainer::elements(int level) const {
  assert(level < m_levels);
  return *(m_elements[level]);
}

inline const DataContainer::EdgesContainer &
DataContainer::edges(int level) const {
  assert(level < m_levels);
  return m_edges[level];
}

inline int DataContainer::numberOfNodes(int level) const {

  assert(level < m_levels);
  return m_nodes[level]->size();
}

inline int DataContainer::numberOfElements(int level) const {

  assert(level < m_levels);
  return m_elements[level]->size();
}

inline int DataContainer::numberOfEdges(int level) const {

  assert(level < m_levels);
  return m_edges[level].size();
}

template <> inline int DataContainer::numberOfEntities<0>(int level) const {

  return numberOfElements(level);
}

template <> inline int DataContainer::numberOfEntities<1>(int level) const {

  return numberOfEdges(level);
}

template <> inline int DataContainer::numberOfEntities<2>(int level) const {

  return numberOfNodes(level);
}

inline const std::array<unsigned int, 3> &
DataContainer::element2Edges(int level, unsigned int elementIndex) const {
  assert(level < m_levels);
  return m_element2Edges[level][elementIndex];
}

inline const std::vector<unsigned int> &
DataContainer::edge2Elements(int level, unsigned int edgeIndex) const {
  assert(level < m_levels);
  return m_edge2Elements[level][edgeIndex];
}

inline const std::vector<unsigned int> &
DataContainer::node2Elements(int level, unsigned int nodeIndex) const {
  assert(level < m_levels);
  return m_node2Elements[level][nodeIndex];
}

inline const std::vector<unsigned int> &
DataContainer::node2Edges(int level, unsigned int nodeIndex) const {
  assert(level < m_levels);
  return m_node2Edges[level][nodeIndex];
}

inline const DataContainer::NodesContainer
DataContainer::getEntityNodes(int codim, int level, int index) const {

  assert(level < m_levels);
  DataContainer::NodesContainer nodes;

  if (codim == 0) {
    // Element
    const auto &nodeIds = this->elements(level)[index];
    for (auto index : nodeIds)
      nodes.push_back(this->nodes(level)[index]);
    return nodes;
  }

  if (codim == 1) {
    // Edge
    const auto &nodeIds = this->edges(level)[index];
    for (auto index : nodeIds)
      nodes.push_back(this->nodes(level)[index]);
    return nodes;
  }

  if (codim == 2) {
    // Vertex
    nodes.push_back(this->nodes(level)[index]);
    return nodes;
  }

  throw std::runtime_error("DataContainer::getEntityNodes(): codim not valid.");
}

inline int
DataContainer::getElementFatherIndex(int level,
                                     unsigned int elementIndex) const {

  assert(level > 0);
  assert(level < m_levels);

  return m_fatherElements[level][elementIndex];
}

inline const std::vector<unsigned int> &
DataContainer::getElementSons(int level, unsigned int elementIndex) const {

  assert(m_levels > 1);

  return m_sonElements[level][elementIndex];
}

template <int cd>
DataContainer::DuneGeometry<cd>
DataContainer::geometry(const Entity<cd> &entity) {

  const auto &levelNodes =
      *(std::get<cd>(m_geometries)[TriangleGrid::entityLevel<cd>(entity)]);
  return DuneGeometry<cd>(levelNodes[TriangleGrid::entityIndex<cd>(entity)]);
}

template <> inline unsigned int DataContainer::id<0>(const Entity<0> &entity) {

  const auto &realEntity = TriangleGrid::getEntityImp<0>(entity);
  return m_elementIds[TriangleGrid::entityLevel<0>(
      entity)][TriangleGrid::entityIndex<0>(entity)];
}

template <> inline unsigned int DataContainer::id<1>(const Entity<1> &entity) {

  const auto &realEntity = TriangleGrid::getEntityImp<1>(entity);
  return m_edgeIds[TriangleGrid::entityLevel<1>(
      entity)][TriangleGrid::entityIndex<1>(entity)];
}

template <> inline unsigned int DataContainer::id<2>(const Entity<2> &entity) {

  const auto &realEntity = TriangleGrid::getEntityImp<2>(entity);
  return m_nodeIds[TriangleGrid::entityLevel<2>(
      entity)][TriangleGrid::entityIndex<2>(entity)];
}

template <int cd> void DataContainer::computeGeometries(int level) {

  std::get<cd>(m_geometries)
      .push_back(shared_ptr<std::vector<Geometry<cd>>>(
          new std::vector<Geometry<cd>>()));
  auto &levelGeometries = *(std::get<cd>(m_geometries)[level]);

  int entityCount = numberOfEntities<cd>(level);
  levelGeometries.reserve(entityCount);

  for (int index = 0; index < entityCount; ++index) {

    Dune::GeometryType geometryType;
    geometryType.makeSimplex(2 - cd);
    levelGeometries.push_back(
        Geometry<cd>(geometryType, getEntityNodes(cd, level, index)));
  }
}
}

#endif
