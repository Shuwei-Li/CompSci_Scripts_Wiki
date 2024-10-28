from ase.build import fcc111, root_surface
atoms = fcc111('Pt', size=(2, 2, 4), a=3.17, vacuum=10, orthogonal=True)
atoms = root_surface(atoms, 3)
atoms.write('POSCAR.root')
