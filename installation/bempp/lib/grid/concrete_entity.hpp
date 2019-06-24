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

#ifndef bempp_concrete_entity_hpp
#define bempp_concrete_entity_hpp

#include "../common/common.hpp"
#include "dune.hpp"

#include "concrete_entity_decl.hpp"
#include "concrete_entity_pointer.hpp"
#include "concrete_range_entity_iterator.hpp"
#include "concrete_subentity_iterator.hpp"

namespace Bempp {

template <typename DuneEntity>
template <int codimSub>
inline typename boost::disable_if_c<
    (codimSub <= DuneEntity::dimension),
    std::unique_ptr<EntityIterator<codimSub>>>::type
ConcreteEntity<0, DuneEntity>::subEntityCodimNIterator() const {
  throw std::logic_error(
      "Entity::subEntityIterator(): invalid subentity codimension");
}

template <typename DuneEntity>
template <int codimSub>
inline
    typename boost::enable_if_c<(codimSub <= DuneEntity::dimension),
                                std::unique_ptr<EntityIterator<codimSub>>>::type
    ConcreteEntity<0, DuneEntity>::subEntityCodimNIterator() const {
  typedef ConcreteSubentityIterator<DuneEntity, codimSub> ConcIterator;
  return std::unique_ptr<EntityIterator<codimSub>>(
      new ConcIterator(m_dune_entity, m_domain_index));
}

template <typename DuneEntity>
std::unique_ptr<EntityPointer<0>>
ConcreteEntity<0, DuneEntity>::father() const {
  if (m_dune_entity.hasFather())
    return std::unique_ptr<EntityPointer<0>>(
        new ConcreteEntityPointer<DuneEntity>(m_dune_entity.father(), m_domain_index));
  else
    return nullptr;
}

/*
template <typename DuneEntity>
std::unique_ptr<EntityIterator<0>>
ConcreteEntity<0, DuneEntity>::sonIterator(int maxlevel) const {
  typedef typename DuneEntity::HierarchicIterator DuneIterator;
  typedef typename DuneEntity::EntityPointer DuneEntityPointer;
  typedef ConcreteRangeEntityIterator<DuneIterator, DuneEntityPointer>
      ConcIterator;
  return std::unique_ptr<EntityIterator<0>>(
      new ConcIterator(m_dune_entity->hbegin(maxlevel),
                       m_dune_entity->hend(maxlevel), m_domain_index));
}
*/

} // namespace Bempp

#endif
