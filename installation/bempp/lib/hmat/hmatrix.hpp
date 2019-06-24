// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_HPP
#define HMAT_HMATRIX_HPP

#include "block_cluster_tree.hpp"
#include "common.hpp"
#include "compressed_matrix.hpp"
#include "data_accessor.hpp"
#include "eigen_fwd.hpp"
#include "hmatrix_compressor.hpp"

#include <mpi.h>
#include <unordered_map>

namespace hmat {

template <typename ValueType> class HMatrixData;
template <typename ValueType, int N> class HMatrix;

template <typename ValueType> using DefaultHMatrixType = HMatrix<ValueType, 2>;

template <typename ValueType, int N> class HMatrix {
public:
  typedef tbb::concurrent_unordered_map<
      shared_ptr<BlockClusterTreeNode<N>>, shared_ptr<HMatrixData<ValueType>>,
      shared_ptr_hash<BlockClusterTreeNode<N>>>
      ParallelDataContainer;

  HMatrix(const shared_ptr<BlockClusterTree<N>> &blockClusterTree,
          MPI_Comm comm = MPI_COMM_WORLD);
  HMatrix(const shared_ptr<BlockClusterTree<N>> &blockClusterTree,
          const HMatrixCompressor<ValueType, N> &hMatrixCompressor,
          MPI_Comm comm = MPI_COMM_WORLD);

  std::size_t rows() const;
  std::size_t columns() const;

  double frobeniusNorm() const;

  void initialize(const HMatrixCompressor<ValueType, N> &hMatrixCompressor);
  bool isInitialized() const;
  void reset();

  void apply(const Eigen::Ref<Matrix<ValueType>> &X,
             Eigen::Ref<Matrix<ValueType>> Y, TransposeMode trans,
             ValueType alpha, ValueType beta) const;

  Matrix<ValueType>
  permuteMatToHMatDofs(const Eigen::Ref<Matrix<ValueType>> &mat,
                       RowColSelector rowOrColumn) const;
  Matrix<ValueType>
  permuteMatToOriginalDofs(const Eigen::Ref<Matrix<ValueType>> &mat,
                           RowColSelector rowOrColumn) const;

  shared_ptr<const BlockClusterTree<N>> blockClusterTree() const;
  shared_ptr<const hmat::HMatrixData<ValueType>>
  data(shared_ptr<const BlockClusterTreeNode<N>> node) const;

  int numberOfDenseBlocks() const;
  int numberOfLowRankBlocks() const;
  int numberOfBlocks() const;

  double memSizeKb() const;

private:
  void apply_impl_parallel(const shared_ptr<BlockClusterTreeNode<N>> &node,
                           const Eigen::Ref<Matrix<ValueType>> &X,
                           Eigen::Ref<Matrix<ValueType>> Y, TransposeMode trans,
                           int levelCount) const;

  void apply_impl_serial(const shared_ptr<BlockClusterTreeNode<N>> &node,
                         const Eigen::Ref<Matrix<ValueType>> &X,
                         Eigen::Ref<Matrix<ValueType>> Y,
                         TransposeMode trans) const;

  double
  frobeniusNorm_impl(const shared_ptr<BlockClusterTreeNode<N>> &node) const;

  bool coarsen_impl(const shared_ptr<BlockClusterTreeNode<N>> &node,
                    double coarsen_accuracy);

  shared_ptr<BlockClusterTree<N>> m_blockClusterTree;
  ParallelDataContainer m_hMatrixData;

  int m_numberOfDenseBlocks;
  int m_numberOfLowRankBlocks;
  int m_memSizeKb;
  int m_nproc;
  int m_rank;
  MPI_Comm m_comm;

  std::vector<shared_ptr<BlockClusterTreeNode<N>>> m_myLeafs;
};
}

#include "hmatrix_impl.hpp"

#endif // HMATRIX_HPP
