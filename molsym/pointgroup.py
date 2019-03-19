#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
"""
Molecular symmetry point group and irreducible representation objects.

Examples
--------
    >>> d2h = D2h()
    >>> print(repr(d2h))
    D2h()

    >>> print(d2h.b2g * d2h.b1u)
    b3u

    >>> print(d2h.b2g**2)
    ag

    >>> print(repr(d2h.b1g))
    IrreducibleRepresentation(pg=D2h(), irrep=(1, 1, -1, -1, 1, 1, -1, -1), degenerate=False)

    >>> print(D2h('b2g') * D2h('b1u'))
    b3u

    >>> d6h = PointGroup('D6h')
    >>> print(repr(d6h))
    PointGroup(pg="dnh", n=6)

    >>> product = d6h.e1g**2
    >>> print(*(str(x) for x in product))
    a1g a2g e2g

    >>> product = d6h.e2u + d6h.e1g * d6h.e2u * d6h.e2u + d6h.a1g
    >>> print(*(str(x) for x in product))
    a1g b1g b2g e1g e1g e1g e2u
"""

from __future__ import division
from bidict import FrozenOrderedBidict
from collections import Iterable
from functools import reduce, total_ordering
from operator import mul
import re


def flatten(items):
    """
    Yield items from any nested iterable.
    Modified from Beazley, D. and B. Jones. Recipe 4.14, Python Cookbook 3rd Ed., O'Reilly Media Inc., Sebastopol, CA (2013)
    by StackOverflow user pylang (cf. https://stackoverflow.com/a/40857703)

    Parameters
    ----------
        items : Iterable
            any (nested) iterable

    Yields
    ------
        all items from the iterable consecutively
    """
    for x in items:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes, IrreducibleRepresentation)):
            yield from flatten(x)
        else:
            yield x


class PointGroup(object):
    """
    Representation of a molecular symmetry point group.

    Parameters
    ----------
        pg : string
            Name of the point group, case-insensitive.
            If a complete name is given, eg. ´D2h´, then the second parameter n is ignored.
            If the name is given as ´Dnh´, the point group order n has to be specified.
        n : int
            Order of the highest-order symmetry operation.

    Raises
    ------
        AttributeError
            If an invalid point group is specified.

    Examples
    --------
        >>> d2h = PointGroup('D2h')
        >>> print(d2h.b2g)
        b2g
        >>> print(d2h('b2g'))
        b2g
    """
    def __init__(self, pg='c', n=1):
        pg = pg.lower()
        if pg.isalpha():
            self.pg = pg
            self.n = n
        else:
            # if pg contains numbers, ignore n and parse pg as specific point group
            regexp_digit = r'\d+'
            self.n = int(re.search(regexp_digit, pg).group())
            self.pg = re.sub(regexp_digit, 'n', pg)
        self._set_irreps()

    def __str__(self):
        """Name of the point group."""
        pg_str = self.pg.title()
        return pg_str.replace('n', str(self.n))

    def __repr__(self):
        """Complete string representation of the point group."""
        return '{}(pg="{}", n={})'.format(self.__class__.__name__, self.pg, self.n)

    def __hash__(self):
        return hash((self.pg, self.n))

    def __eq__(self, other):
        return self.pg == other.pg and self.n == other.n

    def __getattr__(self, symbol):
        """
        Make the irreducible representations of the point group directly
        accessible by their Mulliken symbol as attributes of the PointGroup object.

        Parameters
        ----------
            symbol : str
                The Mulliken symbol of the irreducible representation to be returned.

        Returns
        -------
            IrreducibleRepresentation

        Raises
        ------
            AttributeError
                If the Mulliken symbol does not exist in the point group.

        Example
        -------
            >>> d2h = PointGroup('D2h')
            >>> print(d2h.b1g)
            b1g
        """
        symbol = symbol.lower()
        if symbol in self.elements:
            return self.irrep(symbol)
        else:
            raise AttributeError('character "{}" does not exist in point group {}'.format(symbol, str(self)))

    def __call__(self, symbol):
        """
        Make the PointGroup object callable, return the irreducible representation
        corresponding to the Mulliken symbol string given as parameter.

        Parameters
        ----------
            symbol : str
                The Mulliken symbol of the irreducible representation to be returned.

        Returns
        -------
            IrreducibleRepresentation

        Raises
        ------
            AttributeError
                If the Mulliken symbol does not exist in the point group.

        Example
        -------
            >>> d2h = PointGroup('D2h')
            >>> print(d2h('b1g')
            b1g
        """
        symbol = symbol.lower()
        return self.irrep(symbol)

    # def __iter__(self):
    #     yield from self.elements.values()

    def symbol(self, irrep):
        """
        Return the Mulliken symbol of an irreducible representation.

        Parameters
        ----------
            irrep : IrreducibleRepresentation
                The irreducible representation of which the Mulliken symbol is wanted.

        Returns
        -------
            str
                The Mulliken symbol of ´irrep´.
        """
        return self.elements.inv[irrep]

    def irrep(self, symbol):
        """
        Return the irreducible representation of a Mulliken symbol.

        Parameters
        ----------
            symbol : str
                The Mulliken symbol of the irreducible representation to be returned.

        Returns
        -------
            IrreducibleRepresentation
        """
        return self.elements[symbol]

    @property
    def ts(self):
        """Returns the totally symmetric irreducible representation of the point group."""
        return next(iter(self.elements.inv))

    def _set_irreps(self):
        """Sets the irreducible representations of the point group."""
        # the totally symmetric irrep must always be the first element of elements
        if self.pg == 'cn':
            if self.n == 1:
                # C1
                a = IrreducibleRepresentation(self, (1,))
                self.elements = FrozenOrderedBidict({ 'a': a })
                # besser: dimensions?
                self.symop_multiplicity = (1,)
                self.order = 1

        if self.pg == 'cnv':
            if self.n == 2:
                # C2v
                a1 = IrreducibleRepresentation(self, (1, 1, 1, 1))
                a2 = IrreducibleRepresentation(self, (1, 1,-1,-1))
                b1 = IrreducibleRepresentation(self, (1,-1, 1,-1))
                b2 = IrreducibleRepresentation(self, (1,-1,-1, 1))
                self.elements = FrozenOrderedBidict({ 'a1': a1, 'a2': a2, 'b1': b1, 'b2': b2 })
                self.symop_multiplicity = (1, 1, 1, 1)
                self.order = 4

        if self.pg == 'dn':
            if self.n == 3:
                # D3
                a1 = IrreducibleRepresentation(self, (1, 1, 1))
                a2 = IrreducibleRepresentation(self, (1, 1,-1))
                e =  IrreducibleRepresentation(self, (2,-1, 0), degenerate=True)
                self.elements = FrozenOrderedBidict({ 'a1': a1, 'a2': a2, 'e': e })
                self.symop_multiplicity = (1, 2, 3)
                self.order = 6

        if self.pg == 'dnh':
            if self.n == 2:
                # D2h
                ag =  IrreducibleRepresentation(self, (1, 1, 1, 1, 1, 1, 1, 1))
                b1g = IrreducibleRepresentation(self, (1, 1,-1,-1, 1, 1,-1,-1))
                b2g = IrreducibleRepresentation(self, (1,-1,-1, 1, 1,-1, 1,-1))
                b3g = IrreducibleRepresentation(self, (1,-1, 1,-1, 1,-1,-1, 1))
                au =  IrreducibleRepresentation(self, (1, 1, 1, 1,-1,-1,-1,-1))
                b1u = IrreducibleRepresentation(self, (1, 1,-1,-1,-1,-1, 1, 1))
                b2u = IrreducibleRepresentation(self, (1,-1,-1, 1,-1, 1,-1, 1))
                b3u = IrreducibleRepresentation(self, (1,-1, 1,-1,-1, 1, 1,-1))
                self.elements = FrozenOrderedBidict({ 'ag': ag, 'b1g': b1g, 'b2g': b2g, 'b3g': b3g,
                                                      'au': au, 'b1u': b1u, 'b2u': b2u, 'b3u': b3u })
                self.symop_multiplicity = (1, 1, 1, 1, 1, 1, 1, 1)
                self.order = 8

            if self.n == 6:
                # D6h
                a1g = IrreducibleRepresentation(self, (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
                a2g = IrreducibleRepresentation(self, (1, 1, 1, 1,-1,-1, 1, 1, 1, 1,-1,-1))
                b1g = IrreducibleRepresentation(self, (1,-1, 1,-1, 1,-1, 1, 1,-1,-1,-1, 1))
                b2g = IrreducibleRepresentation(self, (1,-1, 1,-1,-1, 1, 1, 1,-1,-1, 1,-1))
                e1g = IrreducibleRepresentation(self, (2, 1,-1,-2, 0, 0, 2,-1, 1,-2, 0, 0), degenerate=True)
                e2g = IrreducibleRepresentation(self, (2,-1,-1, 2, 0, 0, 2,-1,-1, 2, 0, 0), degenerate=True)
                a1u = IrreducibleRepresentation(self, (1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1))
                a2u = IrreducibleRepresentation(self, (1, 1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1))
                b1u = IrreducibleRepresentation(self, (1,-1, 1,-1, 1,-1,-1,-1, 1, 1, 1,-1))
                b2u = IrreducibleRepresentation(self, (1,-1, 1,-1,-1, 1,-1,-1, 1, 1,-1, 1))
                e1u = IrreducibleRepresentation(self, (2, 1,-1,-2, 0, 0,-2, 1,-1, 2, 0, 0), degenerate=True)
                e2u = IrreducibleRepresentation(self, (2,-1,-1, 2, 0, 0,-2, 1, 1,-2, 0, 0), degenerate=True)
                self.elements = FrozenOrderedBidict({ 'a1g': a1g, 'a2g': a2g, 'b1g': b1g, 'b2g': b2g, 'e1g': e1g, 'e2g': e2g,
                                                      'a1u': a1u, 'a2u': a2u, 'b1u': b1u, 'b2u': b2u, 'e1u': e1u, 'e2u': e2u })
                self.symop_multiplicity = (1, 2, 2, 1, 3, 3, 1, 2, 2, 1, 3, 3)
                self.order = 24


@total_ordering
class IrreducibleRepresentation(object):
    """
    Representation of an irreducible representation in a molecular symmetry point group.

    If two degenerate irreducible representations are multiplied, a list of irrep's is returned.
    This class implements all dunder functions of an Iterable, so that the result of any multiplication
    can be treated as an Iterable without further checking for the number of yielded irrep's.

    Parameters
    ----------
        pg : PointGroup
            Point group this irreducible representation is an element of.
        irrep : tuple of int
            Characters of the irreducible representation.
        degenerate: bool
            True if the irreducible representation is degenerate, else false.

    Examples
    --------
        >>> d6h = PointGroup('D6h')
        >>> b2g = IrreducibleRepresentation(d6h, (1,-1, 1,-1,-1, 1, 1, 1,-1,-1, 1,-1))
        >>> e1g = IrreducibleRepresentation(d6h, (2, 1,-1,-2, 0, 0, 2,-1, 1,-2, 0, 0), degenerate=True)
    """
    def __init__(self, pg, irrep, degenerate=False):
        self.pg = pg
        self.irrep = irrep
        self.degenerate = degenerate

    def _irrep_product(self, other):
        """Product of two irreducible representations from the same point group."""
        return tuple(a*b for a, b in zip(self.irrep, other.irrep))

    def __mul__(self, other):
        """
        Multiplication of an irreducible representation with another irreducible
        representation or a list of irreducible representation from the same point group.

        Parameters
        ----------
            other : (list of) IrreducibleRepresentation
                All multiplied irreducible representations have to be elements of the same point group.

        Returns
        -------
            IrreducibleRepresentation
                If the product is only one irreducible representation.
            list of IrreducibleRepresentation
                If the product yields multiple irreducible representations,
                eg. after multiplication of degenerate irreducible representations
                or multiplication with a list of irreducible representations.

        Raises
        ------
            NotImplemented
                If the irreducible representations are elements of different point groups.
            TypeError
                If a parameter with type other than (list of) IrreducibleRepresentation is passed.

        Example
        -------
        """
        if isinstance(other, IrreducibleRepresentation):
            if self.pg == other.pg:
                product = self._irrep_product(other)
                if self.degenerate and other.degenerate:
                    # The product of two degenerate irrep's is a reducible representation
                    # int() weg, dafür // -> from future import division
                    # eleganter als skalarprodukt schreibbar?
                    ai = [int(1/self.pg.order * sum([i*j*k for i,j,k in zip(self.pg.symop_multiplicity, irrep.irrep, product)])) for irrep in self.pg.elements.values()]
                    product_irreps = [irrep for a, irrep in zip(ai, self.pg.elements.values()) if a != 0]
                    return product_irreps
                else:
                    return IrreducibleRepresentation(pg=self.pg, irrep=product, degenerate=(self.degenerate or other.degenerate))
            else:
                raise NotImplementedError('no automatic lowering of symmetry: cannot multiply irreps of non-identical point groups {} and {}.'.format(self.pg, other.pg))
        elif isinstance(other, list):
            # besser: self * list_product(other) = self * reduce(mul, other)
            # was der konvention entspricht, dass in a x (b x c) die klammer zuerst ausgerechnet wird
            all_products = list(flatten([self*irrep for irrep in other]))
            all_products.sort()
            return all_products
        else:
            raise TypeError('IrreducibleRepresentation can only be multiplied with (list of) IrreducibleRepresentation')

    def __rmul__(self, other):
        if isinstance(other, list):
            return self * other

    def __add__(self, other):
        """Addition of two irreducible representations gives a list of IrreducibleRepresentation."""
        if isinstance(other, list):
            irrep_list = [self] + other
        else:
            irrep_list = [self, other]
        irrep_list.sort()
        return irrep_list

    def __radd__(self, other):
        if isinstance(other, list):
            irrep_list = other + [self]
            irrep_list.sort()
            return irrep_list

    def __pow__(self, other):
        if isinstance(other, int):
            if other < 2:
                return self
            else:
                return reduce(mul, [self] * other)
        else:
            raise ValueError('IrreducibleRepresentation() can only be raised to scalar of type <int>')

    def __str__(self):
        return self.pg.symbol(self)

    def __repr__(self):
        return '{}(pg={!r}, irrep={}, degenerate={})'.format(self.__class__.__name__, self.pg, self.irrep, self.degenerate)

    def __hash__(self):
        return hash((self.pg, self.irrep))

    def __eq__(self, other):
        if isinstance(other, IrreducibleRepresentation):
            return self.pg == other.pg and self.irrep == other.irrep

    def __lt__(self, other):
        """The order of the elements of a point group is only dictated by the canonical convention.
        It is therefore looked up from the point group's ordered dictionary of its elements."""
        if isinstance(other, IrreducibleRepresentation):
            for irrep in self.pg.elements.values():
                if irrep == self:
                    return True
                elif irrep == other:
                    return False

    def __len__(self):
        return 1

    def __contains__(self, item):
        return self == item

    def __iter__(self):
        yield self

class PointGroupABC(PointGroup):
    """
    Abstract base class for the named point group classes.

    The named point group classes inherit from this class and provide an abbreviation
    to the PointGroup() class, eg. C2v() as a shorthand for PointGroup('C2v'). The desired
    point group is determined by the __name__ attribute of the inheriting class.

    This class also implements in its constructor the functionality of the named point
    group classes to directly return an irreducible representation when called with
    a Mulliken symbol that represents an element of the point group.

    Parameters
    ----------
        mulliken_symbol : str, optional
            Mulliken symbol of the IrreducibleRepresentation that is to be returned.

    Returns
    -------
        PointGroup
            If mulliken_symbol is None.
        IrreducibleRepresentation
            If a Mulliken symbol of the specified point group is passed.

    Examples
    --------
        >>> D2h() == PointGroup('D2h')
        True
        >>> print(D2h('b3g'))
        b3g
        >>> print(repr(D2h('b3g')))
        IrreducibleRepresentation(pg=D2h(), irrep=(1, -1, 1, -1, 1, -1, -1, 1), degenerate=False)
    """
    def __new__(cls, mulliken_symbol=None):
        pgobj = super().__new__(cls)
        super().__init__(pgobj, pg=cls.__name__)
        if mulliken_symbol:
            # get the IrreducibleRepresentation of ´mulliken_symbol´
            return getattr(pgobj, mulliken_symbol)
        else:
            return pgobj

    def __init__(self):
        # initializiation is already done in __new__(), but called by interpreter anyway
        pass

    def __repr__(self):
        return '{}()'.format(self.__class__.__name__)


# Shorthands to the PointGroup object via named point group classes
class C1(PointGroupABC): pass
class C2v(PointGroupABC): pass
class D3(PointGroupABC): pass
class D2h(PointGroupABC): pass
class D6h(PointGroupABC): pass
