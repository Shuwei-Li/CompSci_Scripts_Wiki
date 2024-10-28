import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import csv

######### 1. Extract the adsorption energy from 'scf.pwo', and distance between gas molecule and surface from 'scf.pwi' ##########
regex_energy = re.compile(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry")
regex_coordinates = re.compile(r"\s+(N|H|S)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)")
E_sur_NH3 = []
Eads = []
S_heights = []
NH3_heights = []
E_sur = -2343.6746789321      ## Ry, the energy of bilayer MoS2
E_NH3 = -22.8129700664        ## Ry,  the energy of NH3 molecule

for i in [1,2,3,4,5,6,7,8,15,16,17,18,19,20,21,22,23,24]:        ## the numbers mean which directory you will read
    directory = f"{i:02d}"
    if os.path.exists(directory):
        pwi_filepath = os.path.join(directory, "scf.pwi")
        if not os.path.exists(pwi_filepath):
            print(f"No scf.pwi file found in {directory} directory")
            continue
        with open(pwi_filepath, 'r') as f:
            pwi_content = f.read()
            coordinates = regex_coordinates.findall(pwi_content)
            if not coordinates:
                print(f"No N, H or S atoms found in {directory}/scf.pwi file")
                continue
            N_coordinates = []
            H_coordinates = []
            S_coordinates = []
            for atom, x, y, z in coordinates:
                x, y, z = float(x), float(y), float(z)
                if atom == "N":
                    N_coordinates.append([x, y, z])
                elif atom == "H":
                    H_coordinates.append([x, y, z])
                elif atom == "S":
                    S_coordinates.append([x, y, z])
            # Calculate center of mass for NH3
            mass_N = 14.01
            mass_H = 1.01
            mass_total = 3 * mass_H + mass_N
            NH3_center_of_mass = (
                (mass_N * np.mean(N_coordinates, axis=0)) +
                (3 * mass_H * np.mean(H_coordinates, axis=0))
            ) / mass_total
            # Get the height of the center of mass
            NH3_heights.append(NH3_center_of_mass[2])
            # Get the height of the highest 'S' layer
            max_S_height = max(S_coordinates, key=lambda coord: coord[2])[2]
            tolerance = 0.5 ## define a tolerance value
            S_height_range = [max_S_height - tolerance, max_S_height + tolerance] ## the range for the top layer of 'S'
            # Get the height of the 'S' atoms within the defined range
            S_heights_within_range = [coord[2] for coord in S_coordinates if S_height_range[0] <= coord[2] <= S_height_range[1]]
            # Calculate the average height of the 'S' atoms within the defined range
            avg_S_height = sum(S_heights_within_range) / len(S_heights_within_range)
            S_heights.append(avg_S_height)
        pwo_filepath = os.path.join(directory, "scf.pwo")
        if not os.path.exists(pwo_filepath):
            print(f"No scf.pwo file found in {directory} directory")
            continue
        with open(pwo_filepath, 'r') as f:
            pwo_content = f.read()
            energy_match = regex_energy.search(pwo_content)
            if not energy_match:
                print(f"No total energy found in {directory}/scf.pwo file")
                continue
            energy = float(energy_match.group(1))
            ads_energy = (energy - E_sur - E_NH3) * 13.605693009
            E_sur_NH3.append(energy)
            Eads.append(ads_energy)
dis=[]
for i in range(len(S_heights)):
    distance=NH3_heights[i]-S_heights[i]
    dis.append(distance)

########### 2. Using above adsorption energy and distance to fit 'Morse potential' #########
E_sur_NH3_equi = -2366.4974202699   ## unit, Ry 
Eads_equi = (E_sur_NH3_equi - E_sur -E_NH3) * 13.605693009   ## unit, eV, the adsorption energy of NH3 on bilayer MoS2
dis_equi = dis[0]-0.2  ## unit, angstrom
print('Eads_equi: ', Eads_equi)
print('dis_equi: ',dis_equi)
Eads.append(Eads_equi)   ## save the adsorption energy to the list
dis.append(dis_equi)   ## save the equilibrium distance to the list
rows = zip(dis, Eads) ## use zip() to combine the two lists into tuples
with open('morse_fit.dat', 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')   ## use the csv module to write the tuples to the file
    writer.writerow(['distance/anstrom', 'adsorption_energy/eV'])
    for row in rows:
        writer.writerow(row)

De = -Eads_equi*1.602176634E-19 ## the adosorption energy, unit: J, convert eV to J
z_e = dis_equi  ## angstrom, the corresponding distance between NH3 and MoS2
def func(z,y1):
    return De*(np.exp(-2*y1*(z-z_e))-2*np.exp(-y1*(z-z_e)))
Eads = [i*1.602176634E-19 for i in Eads]
popt, pcov = curve_fit(func, dis, Eads)
print(popt)
y_fit = popt[0]
print("y_fit:", y_fit)
z = np.arange(1, 5, 0.1)
E_fit = []
for z_x in z:
    E_ads_fitting = De*(np.exp(-2*y_fit*(z_x-z_e))-2*np.exp(-y_fit*(z_x-z_e)))
    E_fit.append(E_ads_fitting)

y_fit_m = y_fit*1E10  ## the fitting parameter, convert the angstrom to m, we will use this one for subsequent calculation for adsorption density
z_m = [i*1E-10 for i in z]  ## the distance, unit: m
z_e_m = z_e*1E-10 ## the equilibrium distance, unit: m
E_fit_m = []
for z_x in z_m:
    E_ads_fitting_m = De*(np.exp(-2*y_fit_m*(z_x-z_e_m))-2*np.exp(-y_fit_m*(z_x-z_e_m)))
    E_fit_m.append(E_ads_fitting_m)

plt.plot(dis, Eads, 'o', mfc='blue')
plt.plot(z, E_fit, 'k-')
plt.plot(z, E_fit_m, 'r--')
plt.xlabel('Distance', fontsize=12)
plt.ylabel('Adsorption energy', fontsize=12)
plt.show()

############# 3. Calculate the adsorption density #############
### Parameters or constant for calculation ###
h_bar = 1.054571800E-34   ### reduced Planck constant, unit: J*s
k = 1.380649E-23  ### Boltzmann constant, unit: J/K
mg = (mass_N + 3*mass_H)*1.660539040E-27  ### the mass of gas molecule, unit: Kg
T = 300  ### the temperature, unit: K
Pg_ppm = [0.001, 0.01, 0.1, 1, 10, 100, 200, 500, 1000, 5] ### partial pressure, unit: ppm
Pg = [i*1E-6*1E5 for i in Pg_ppm]  ### convert the ppm to Pa
wg = y_fit_m*(2*De/mg)**0.5  ### the characteristic frequency of a gas molecule with mass mg
delta = De/h_bar/wg
max_n = int(np.round(2*delta - 0.5))   ### calculate the maximum value of n in En
print('max_n: ', max_n)
lam = (2*math.pi*h_bar**2/mg/k/T)**0.5   ### the thermal wavelength, m
print('the thermal wavelength (m) of NH3 is : ', lam)

### Define the function E(n) which is used to calculate the eigenvalues En
def E(n):
    En = -h_bar*wg*(-2*(2*n+1)+8*delta)**2/64/delta
    return En

### Define the function q_sum which is used to calculate the sum of exp(-En/k/T)
def q_sum(n):
    q_sum = 0
    for i in range(1, n+1):
        q_sum_i = np.exp(-E(i)/k/T)
        q_sum = q_sum + q_sum_i
    return q_sum

ng = []
for i in Pg:
    ng_i = i*lam/k/T*np.exp(De/k/T)*q_sum(max_n)*1E-4  ### calculate the adsorption density for each Pg, convert m-2 to cm-2
    ng.append(ng_i)

print('the temperature is: ', T)
print('the gas concentration (ppm) is: ', Pg_ppm)
print('the gas adsorption density (cm-2) is: ')
for i  in ng:
    print('{:.6E}'.format(i))
