### This is a python script that can be used to add adsorbates on surface automatically. ###
import os
from ase.io import read, write
import numpy as np

def count_oxygen_neighbors(poscar_path, metal_atoms, threshold_distance, layer_thickness):
    # Read the POSCAR file
    structure = read(poscar_path)

    # Determine the z-coordinate range for the topmost layer
    z_coordinates = [atom.position[2] for atom in structure if atom.symbol in metal_atoms]
    max_z = max(z_coordinates)
    min_z_for_top_layer = max_z - layer_thickness

    # Identify the topmost metal atoms
    topmost_metal_atoms = [atom.index for atom in structure if atom.symbol in metal_atoms and atom.position[2] >= min_z_for_top_layer]

    # Initialize count dictionary for topmost metal atoms
    oxygen_neighbors_count = {index: 0 for index in topmost_metal_atoms}

    # Calculate distances and count neighbors for topmost metal atoms
    for index in topmost_metal_atoms:
        distances = structure.get_distances(index, range(len(structure)), mic=True)
        for j, distance in enumerate(distances):
            if structure[j].symbol == 'O' and distance <= threshold_distance and distance > 0:
                oxygen_neighbors_count[index] += 1

    # Find the metal atoms with the minimum number of oxygen neighbors
    min_neighbors = min(oxygen_neighbors_count.values())
    min_indices = [index for index, count in oxygen_neighbors_count.items() if count == min_neighbors]

    # Filter out duplicate metal types, keeping only one of each type
    unique_metal_types = set()
    unique_min_indices = []
    for index in min_indices:
        metal_type = structure[index].symbol
        if metal_type not in unique_metal_types:
            unique_metal_types.add(metal_type)
            unique_min_indices.append(index)

    return oxygen_neighbors_count, unique_min_indices

def add_O2(poscar_path, index, directory):
    # Read the POSCAR file
    slab = read(poscar_path)
    total_atom_number = slab.get_global_number_of_atoms()
    metal_ads_site_pos = slab.positions[index]
    slab.append('O')  # add the first oxygen atom
    slab.positions[total_atom_number] = metal_ads_site_pos + [0, 0, 1.8]
    slab.append('O')  # add the second oxygen atom
    slab.positions[total_atom_number + 1] = slab.positions[total_atom_number] + [0, 0, 1.0]

    # Save in the specified directory
    slab.write(f'{directory}/POSCAR_O2.vasp')

def add_OOH(poscar_path, index, directory):
    # Read the POSCAR file
    slab = read(poscar_path)
    total_atom_number = slab.get_global_number_of_atoms()
    metal_ads_site_pos = slab.positions[index]
    slab.append('O')  # add the first oxygen atom
    slab.positions[total_atom_number] = metal_ads_site_pos + [0, 0, 1.8]
    slab.append('O')  # add the second oxygen atom
    slab.positions[total_atom_number + 1] = slab.positions[total_atom_number] + [0, 0, 1.0]
    slab.append('H')  # add the hydrogen atom
    slab.positions[total_atom_number + 2] = slab.positions[total_atom_number + 1] + [0, 0, 0.9]

    # Save in the specified directory
    slab.write(f'{directory}/POSCAR_OOH.vasp')

def create_directory_and_save_files(poscar_path, min_indices):
    upper_level_path = os.path.join(os.path.dirname(poscar_path), "..")  # One level up from the POSCAR file's directory

    for index in min_indices:
        metal_type = read(poscar_path)[index].symbol
        directory = os.path.join(upper_level_path, f"{metal_type}_{index}_adsorption")

        # Create directory if it doesn't exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        # Add O2 and OOH and save in the respective directories
        add_O2(poscar_path, index, directory)
        add_OOH(poscar_path, index, directory)
        print(f"Files saved in {directory}")

# Example usage
poscar_path = './POSCAR'  # Replace with your POSCAR file path
metal_atoms = ['Fe', 'Co', 'Ni', 'Zn', 'Sn', 'Sb', 'Bi']  # Replace with your list of metal atom symbols
threshold_distance = 2.5  # Replace with your threshold distance in Angstrom
layer_thickness = 0.5  # Thickness of the topmost layer in Angstrom

neighbors_count, min_indices = count_oxygen_neighbors(poscar_path, metal_atoms, threshold_distance, layer_thickness)
print("Oxygen neighbors count:", neighbors_count)
print("Unique indices of metal atoms with minimum oxygen neighbors:", min_indices)

create_directory_and_save_files(poscar_path, min_indices)
