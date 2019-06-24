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

#ifndef bempp_concrete_domain_index_hpp
#define bempp_concrete_domain_index_hpp

#include "domain_index.hpp"
#include "concrete_index_set.hpp"

namespace Bempp {

template <typename DuneGrid> class ConcreteDomainIndex : public DomainIndex {
public:
  ConcreteDomainIndex(const DuneGrid &grid,
                      const std::vector<int> &elementIndexToPhysicalEntity)
      : DomainIndex(std::unique_ptr<IndexSet>(
                        new ConcreteIndexSet<typename DuneGrid::LevelGridView>(
                            &grid.levelIndexSet(0))),
                    elementIndexToPhysicalEntity) {}
};

} // namespace

#endif
