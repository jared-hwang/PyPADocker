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

#ifndef bempp_concrete_entity_pointer_hpp
#define bempp_concrete_entity_pointer_hpp

#include "../common/common.hpp"

#include "entity_pointer.hpp"
#include "concrete_entity_decl.hpp"

namespace Bempp {

/**
 \ingroup grid_internal
 \brief Wrapper of a Dune entity pointer of type \p DuneEntityPointer.

 EntityPointers have been removed in newer Dune versions. This class is only
 left for compatibility.
 */
template <typename DuneEntity_>
class ConcreteEntityPointer
    : public EntityPointer<DuneEntity_::codimension> {
private:
  typedef DuneEntity_ DuneEntity;
  ConcreteEntity<DuneEntity::codimension, DuneEntity> m_entity;

public:
  /** \brief Constructor */
  ConcreteEntityPointer(const DuneEntity &dune_entity,
                        const DomainIndex &domain_index)
      : m_entity(dune_entity, domain_index) {
  }

  virtual const Entity<DuneEntity::codimension> &entity() const {
    return m_entity;
  }
};

} // namespace Bempp

#endif
