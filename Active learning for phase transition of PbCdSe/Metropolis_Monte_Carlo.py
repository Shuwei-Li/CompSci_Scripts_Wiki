import numpy as np
from ase.io import read,write
from ase.neighborlist import NeighborList
import ase.db
import numpy as np
from tpot import TPOTRegressor
from sklearn.model_selection import train_test_split
import math
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.gaussian_process as gp
from sklearn.pipeline import Pipeline
import pickle
from random import choice
import copy
from ase import Atoms
from ase import units

db=ase.db.connect('cdpbs.db')
#db1=ase.db.connect('cdpbs_mc.db')
#db2=ase.db.connect('cdpbs_mc_min.db')
ini_config = db.get(id=2)

# print(ini_config.volume)

# feature extract
list = []
distance = []
P = read('pure.traj')

# P_c = read('pure.traj')

for Pi in P:
    if Pi.symbol == 'Cd':
        list.append(Pi.index)
# print('list',list)
a = list[0]

for l in list:
    distance.append(P.get_distance(a,l))
# print(sorted(distance))
distance_first = sorted(distance)[1]
radius_first = (distance_first + 0.1)/2
distance_second = sorted(distance)[4]
radius_second = (distance_second + 0.1)/2

neighbor_first_list = []
neighbor_second_list = []
for g_atom in list:
    nl = NeighborList([radius_first]*216, self_interaction = False, bothways = True)
    nl.update(P)
    indices, offsets = nl.get_neighbors(g_atom)
    nlist = []
    nlist.append(g_atom)
    for i in indices:
        if P[i].symbol == 'Cd':
            nlist = nlist + [i]
    neighbor_first_list.append(nlist)
    #the format of "neighbor_first_list" is : the index of guest atom + the index of its first neighbor atom

    nl2 = NeighborList([radius_second]*216, self_interaction = False, bothways = True)
    nl2.update(P)
    indices, offsets = nl2.get_neighbors(g_atom)
    n2list = []
    n2list.append(g_atom)
    for i in indices:
        if P[i].symbol == 'Cd':
            n2list = n2list + [i]
    for j in n2list[1:]:
        if j in nlist[1:]:
            n2list.remove(j)
    neighbor_second_list.append(n2list)

# print('first',neighbor_first_list)
# print('second',neighbor_second_list)

sid = ini_config.sid
begin_index = sid.find('0x')
config_code = sid[begin_index:-3]
config_str = '{:0216b}'.format(int(config_code,16))
# print(config_str)
config_int = [int(i) for i in config_str]
# print(len(config_int))
# print(config_int)
# P_symbol = P.get_chemical_symbols()


# from here, carry out the loop compute
# K = 1.380649e-4/1.602176620898
T = 100

N_sample = 200000
Energy = []
Energy_min = []

config_int_MC_steps = []
config_int_MC_steps.append(config_int)

config_int_min = []
config_int_min.append(config_int)

# X_pre = [[0]*41]*(N_sample+10)
for a in range(N_sample):
    # print('sid-0', config_int)
    P_symbol_1 = []
    P_symbol_2 = []
    # for i in range(216):
    #     if P[i].symbol == 'Cd':
    #         P_symbol_1.append('Cd')
    #         P_symbol_2.append('Cd')
    for i in range(216):
        P_symbol_1.append(P[i].symbol)
        P_symbol_2.append(P[i].symbol)
    # P_symbol_1 = copy.deepcopy(P.get_chemical_symbols())

    for i in range(len(config_int)):
        if config_int[i] == 1:
            P_symbol_1[i] = 'Pb'
    # print('P_symbol_1',P_symbol_1)
    # print(len(config_int))
    n_pb_1 = P_symbol_1.count('Pb')
    pb_neigh_first_num_1 = [0] * 13
    cd_neigh_first_num_1 = [0] * 13
    for row in neighbor_first_list:
        row = [int(i) for i in row]
        row_atomic_symbol_1 = [P_symbol_1[i] for i in row]
        n_pb_row_first_1 = row_atomic_symbol_1[1:].count('Pb')
        if row_atomic_symbol_1[0] == 'Pb':
            pb_neigh_first_num_1[n_pb_row_first_1] += 1
        else:
            cd_neigh_first_num_1[n_pb_row_first_1] += 1

    pb_neigh_second_num_1 = [0] * 7
    cd_neigh_second_num_1 = [0] * 7

    for row in neighbor_second_list:
        row = [int(i) for i in row]
        row_atomic_symbol_11 = [P_symbol_1[i] for i in row]
        n_pb_row_second_1 = row_atomic_symbol_11[1:].count('Pb')
        if row_atomic_symbol_11[0] == 'Pb':
            pb_neigh_second_num_1[n_pb_row_second_1] += 1
        else:
            cd_neigh_second_num_1[n_pb_row_second_1] += 1
    atom_feature_1 = [n_pb_1] + pb_neigh_first_num_1 + cd_neigh_first_num_1 + pb_neigh_second_num_1 + cd_neigh_second_num_1
    # print(atom_feature)
    Atom_feature_1 = []
    Atom_feature_1.append(atom_feature_1)
    fr = open('best_model', 'rb+')
    model = pickle.load(fr)
    fr.close()
    X_pre_1 = Atom_feature_1
    energy = model.predict(X_pre_1)
    # print(energy)
    # print(Atom_feature_1)
    Energy.append(energy[0])


    Cd_list = []
    Pb_list = []

    for i in range(len(config_int)):
        if config_int[i] == 1:
            Pb_list.append(i)
        # else:
        #     if P[i].symbol == 'Cd':
        #         Cd_list.append(i)



    # print('sid-0', config_int)

    new_config_int = [i for i in config_int]

    # new_config_int = config_int

    # print('sid-1', config_int)

    # ex_cd = choice(Cd_list)
    while len(Cd_list)==0:
        ex_pb = choice(Pb_list)
        for row in neighbor_first_list:
            if row[0] == ex_pb:
                first_neighbor_list_atom = [row[i] for i in range(1,13)]
                for i in first_neighbor_list_atom:
                    if config_int[i] == 0:
                        Cd_list.append(i)

    ex_cd = choice(Cd_list)
    new_config_int[ex_cd] = 1
    new_config_int[ex_pb] = 0

    # print('sid-2', config_int)
    # print('sid-n', new_config_int)

    # P_symbol_2 = copy.deepcopy(P.get_chemical_symbols())


    for i in range(len(new_config_int)):
        if new_config_int[i] == 1:
            P_symbol_2[i] = 'Pb'
    # print('P_symbol_2', P_symbol_2)
    n_pb_2 = P_symbol_2.count('Pb')
    # print(n_pb_2)
    pb_neigh_first_num_2 = [0] * 13
    cd_neigh_first_num_2 = [0] * 13

    for row in neighbor_first_list:
        row = [int(i) for i in row]
        row_atomic_symbol_2 = [P_symbol_2[i] for i in row]
        n_pb_row_first_2 = row_atomic_symbol_2[1:].count('Pb')
        if row_atomic_symbol_2[0] == 'Pb':
            pb_neigh_first_num_2[n_pb_row_first_2] += 1
        else:
            cd_neigh_first_num_2[n_pb_row_first_2] += 1

    pb_neigh_second_num_2 = [0] * 7
    cd_neigh_second_num_2 = [0] * 7
    Atom_feature_2 = []

    for row in neighbor_second_list:
        row = [int(i) for i in row]
        row_atomic_symbol_22 = [P_symbol_2[i] for i in row]
        n_pb_row_second_2 = row_atomic_symbol_22[1:].count('Pb')
        if row_atomic_symbol_22[0] == 'Pb':
            pb_neigh_second_num_2[n_pb_row_second_2] += 1
        else:
            cd_neigh_second_num_2[n_pb_row_second_2] += 1

    atom_feature_2 = [n_pb_2] + pb_neigh_first_num_2 + cd_neigh_first_num_2 + pb_neigh_second_num_2 + cd_neigh_second_num_2
    Atom_feature_2.append(atom_feature_2)
    # print(Atom_feature_2)


    fr = open('best_model', 'rb+')
    model = pickle.load(fr)
    fr.close()
    X_pre_2 = Atom_feature_2
    energy_2 = model.predict(X_pre_2)
    # print('test-energy', energy)

    # print('energya+1',energy)
    #Energy.append(energy[a+1])
    energy_new = energy_2
    # print('energy',Energy)
    # print('energy_new',energy_new)
    # print('energy[0]',energy[0])

    delt_energy = energy[0] - energy_new
    theta = np.random.rand()
    mi = delt_energy/(units.kB*T)
    mi = np.float(mi)
    p_a = np.exp(mi)
    panduan = theta - p_a
    # print('puan', panduan)
    # print('delt_energy',delt_energy)
    if delt_energy<0 and panduan>0:
        pass
    else:
        config_int = new_config_int

    config_int_MC_steps.append(config_int)

    if delt_energy > 0:
        config_int_min.append(config_int)






    # print(config_int)



np.savetxt('energy.csv',Energy,delimiter=',')

for i in range(1):
    Energy_min.append(Energy[i])
# print(type(Energy_min[0]))

for i in range(N_sample-1):
    if Energy[i+1] < Energy_min[i]:
        Energy_min.append(Energy[i+1])
    else:
        Energy_min.append(Energy_min[i])

np.savetxt('energy_min.csv',Energy_min,delimiter=',')


Plot = []
# Energy_min_plot = []
for i in range(N_sample):
    if i % 100 == 0:
        Plot.append([i,Energy[i],Energy_min[i]])
        # Energy_min_plot.append(Energy_min[i])
np.savetxt('Plot.csv',Plot,delimiter=',')
# np.savetxt('energy_min_plot.csv',Energy_min_plot,delimiter=',')


np.savetxt('config_int_MC_steps.csv',config_int_MC_steps,delimiter=',')
np.savetxt('config_int_min.csv',config_int_min,delimiter=',')


#alpha = 7.73/100

## Storage the configurations along with MC steps

# for i in range(N_sample):
#     P_c = P.copy()
#     for a in range(216):
#
#         if config_int_MC_steps[i][a] == 1:
#             P_c[a].symbol = 'Pb'
#     scale_cell = P_c.get_cell()*(1 + alpha * sum(config_int_MC_steps[i])/108)
#     P_c.set_cell(scale_cell,scale_atoms=True)
#     db1.write(P_c,MC_id=i,MC_energy=Energy[i])

## Storage the configurations for the minimal energy

#P_c2 = P.copy()
#for a in range(216):
#    if config_int_min[-1][a] == 1:
#        P_c2[a].symbol = 'Pb'
## print(P_c2)
#scale_cell2 = P_c2.get_cell()*(1 + alpha * sum(config_int_min[-1])/108)
#P_c2.set_cell(scale_cell2,scale_atoms=True)
#db2.write(P_c2,MC_energy=Energy_min[-1])





# print(energy)
# print(energy_2)


# print(Energy)
# print(np.size(Atom_feature,1))
# print(Energy[0])
# print(Energy[1])
