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

#ifndef bempp_domain_index_hpp
#define bempp_domain_index_hpp

#include "../common/common.hpp"
#include <memory>
#include <vector>
#include <utility>

#include "index_set.hpp"

namespace Bempp {

class DomainIndex {
public:
  DomainIndex(std::unique_ptr<IndexSet> &&level0IndexSet,
              const std::vector<int> &elementIndexToPhysicalEntity)
      : m_level0IndexSet(std::move(level0IndexSet)),
        m_domainIndices(elementIndexToPhysicalEntity) {}

  int domain(const Entity<0> &entity) const;
  std::vector<int> domainIndices() const;

private:
  std::unique_ptr<IndexSet> m_level0IndexSet;
  std::vector<int> m_domainIndices;
};

} // namespace

#endif
