import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from .context import molsym

d2h = molsym.D2h()
print(repr(d2h))
# D2h()

product = d2h.b2g * d2h.b1u
print(product)
# b3u

print(len(product))
# 1

print(d2h.ag in product)
# False

print(d2h.b3u in product)
# True

for irrep in product:
    print('iter:', irrep)

print(d2h.b2g**2)
# ag

print(repr(d2h.b1g))
# IrreducibleRepresentation(pg=D2h(), irrep=(1, 1, -1, -1, 1, 1, -1, -1), degenerate=False)

print(molsym.D2h('b2g') * molsym.D2h('b1u'))
# b3u

d6h = molsym.PointGroup('D6h')
print(repr(d6h))
# PointGroup(pg="dnh", n=6)

product = d6h.e1g**2
print(*(str(x) for x in product))
# a1g a2g e2g

product = d6h.e2u + d6h.e1g * d6h.e2u * d6h.e2u + d6h.a1g
print(*(str(x) for x in product))
# a1g b1g b2g e1g e1g e1g e2u

print(len(product))
# 1

print(d6h.a1u in product)
# False

print(d6h.b2g in product)
# False

for irrep in product:
    print('iter:', irrep)
