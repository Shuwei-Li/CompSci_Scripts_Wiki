from ase.build import surface
from ase.io import read
from ase.io import write
cell = read('POSCAR')
s2 = surface(cell, (1,0,0), 2)
s2 = s2*(3, 3, 1)
s2.center(vacuum=10, axis=2)
s2.write('surface.traj')
