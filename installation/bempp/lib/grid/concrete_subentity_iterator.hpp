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

#ifndef bempp_concrete_subentity_iterator_hpp
#define bempp_concrete_subentity_iterator_hpp

#include "../common/common.hpp"

#include "entity_iterator.hpp"
#include "concrete_entity_decl.hpp"

namespace Bempp {

/** \ingroup grid_internal
 *  \brief Iterator over the subentities of codimension \p codimSub
 *  of a given Dune entity of type \p DuneEntity. */
template <typename DuneEntity, int codimSub>
class ConcreteSubentityIterator : public EntityIterator<codimSub> {
  static_assert(DuneEntity::codimension == 0,
                     "ConcreteSubentityIterator: only codim-0 entities "
                     "support iteration over subentities");
  static_assert((int)DuneEntity::codimension <
                         (int)ConcreteSubentityIterator::codimension,
                     "ConcreteSubentityIterator: subentity codimension "
                     "must exceed entity codimension");

public:
  typedef typename DuneEntity::template Codim<codimSub>::Entity DuneSubentity;

private:
  const DuneEntity m_dune_entity;
  std::unique_ptr<ConcreteEntity<ConcreteSubentityIterator::codimension, DuneSubentity>>
      m_subentity;
  const DomainIndex &m_domain_index;
  int m_cur_n;

  void updateFinished() {
    const int codim = codimSub;
    this->m_finished = (m_cur_n == m_dune_entity.subEntities(codim));
  }

  void updateSubentity() {
    if (!this->finished()) {
      m_subentity.reset(
        new ConcreteEntity<ConcreteSubentityIterator::codimension, DuneSubentity>(
          m_dune_entity.template subEntity<codimSub>(m_cur_n), m_domain_index));
    }
  }

public:
  /** \brief Constructor.

      \param dune_entity Entity whose subentities will be iterated over */
  explicit ConcreteSubentityIterator(const DuneEntity &dune_entity,
                                     const DomainIndex &domain_index)
      : m_dune_entity(dune_entity),
        m_subentity(new ConcreteEntity<ConcreteSubentityIterator::codimension, DuneSubentity>(dune_entity.template subEntity<codimSub>(0), domain_index)),
        m_domain_index(domain_index), m_cur_n(0) {
    updateFinished();
  }

  virtual void next() {
    ++m_cur_n;
    updateFinished();
    updateSubentity();
  }

  virtual const Entity<ConcreteSubentityIterator::codimension> &entity() const {
    return *m_subentity;
  }

  virtual std::unique_ptr<EntityPointer<ConcreteSubentityIterator::codimension>>
  frozen() const {
    const int codim = ConcreteSubentityIterator::codimension;
    return std::unique_ptr<EntityPointer<codim>>(
        new ConcreteEntityPointer<DuneSubentity>(
            m_subentity->duneEntity(), m_domain_index));
  }
};

} // namespace Bempp

#endif
