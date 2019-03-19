# molsym

Molecular symmetry point group and irreducible representation objects supporting algebraic operations for Python.

This package is under development. Names and functionality are likely to change.

Currently supported point groups:
C1, C2v, D3, D2h, D6h.

Additional point groups can be easily added in `pointgroup.py`.

Resources used for creating this package:
* J.M. Hollas, "Symmetry in Molecules," Chapman&Hall, London (1972)
* [Gernot Katzer, "Character Tables for Point Groups used in Chemistry"](http://www.gernot-katzers-spice-pages.com/character_tables/)


Examples
--------
```
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
````

TODO / Backlog
--------------
    * General:
        * add __all__ = [...]
        * unit tests

    * `PointGroup`:
        * implement comparison methods by group order and self.n
              x.__lt__(y), x<=y calls x.__le__(y), x==y calls x.__eq__(y), x!=y calls x.__ne__(y),
              x>y calls x.__gt__(y), and x>=y calls x.__ge__(y).
        * add attribute to return totally symmetric irrep
        * in __init__:
              * raise AttributeError if pg invalid
              * check lower bound for n depending on point group
        * generate irreps and character table from the generating symmetry operations

    * `IrreducibleRepresentation`:
        * remove `degenerate` parameter, not neccessary:
            Degeneracy = Dimension = character of symmetry operation `E`
            simply reduce every Representation before multiplication, then
            degeneracy does not matter anymore
        * implement correct Mulliken symbols:
            1-dimensional characters `A` und `B`
            2-dimensional charakters `E`
            3-, 4-, 5-dimensional charakters `T`, `G`, `H`
        * lowering of symmetry when multiplied with an irrep from a different point group, if possible.

    * `SymmetryOperation` class
        * analytic implementation of symmetry operation application
        * can then be used to generate all point groups
