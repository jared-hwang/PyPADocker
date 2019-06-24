"""Define the EntityIterator class."""


#pylint: disable=too-few-public-methods

class EntityIterator(object):
    """Implements iterators over entities."""

    def __init__(self, codimension, impl):
        """Constructor. Should not be called by the user."""
        self._codimension = codimension
        self._impl = impl

    def next(self):
        """Return the next element."""
        from bempp.api.grid import entity as _entity
        return _entity.Entity(self._codimension,
                              self._impl.__next__())

    def __next__(self):
        """Special method. Return the next element."""
        return self.next()

    def __iter__(self):
        """Special method for the Python iterator concept."""
        return self
