import numpy as np
import math

#### The input parameters for calculate the carrier concentration ####
me_eff = 0.6211492033381033     ### the effective mass of electron, unit: me, the mass of an electron
mp_eff = 1.461381561903015     ### the effective mass of hole, unit: me, the mass of an electron
me = 9.1093837015E-31   ### the mass of an eletron, unit: Kg
h = 6.62607015E-34  ### Planck constant, unit: J*s
k = 1.380649E-23   ### Boltzmann constant, unit: J/K
T = [300]   ### temperature, unit: K
band_gap = 2.35969 - 0.93995 ### bandgap, unit: eV
band_gap = band_gap*1.602176634E-19 ## bandgap, unit: J
Ec = 2.35969 ### CBM, unit: eV
Ec = Ec*1.602176634E-19 ### unit: J
Ev = 0.93995 ### VBM, unit: eV
Ev = Ev*1.602176634E-19 ### unit: J
Ef = 2.2116296 ### Fermi level, unit: eV Note: this value was obtained via fitting the calculated carrier concentration with that of experimental resutls
Ef = Ef*1.602176634E-19 ### unit: J

Ef_change_down = (0.44346-1.05132)/2 - (0.43029-0.98596)/2
Ef_NO2 = Ef + Ef_change_down*1.602176634E-19

Ef_change_up = (0.53895-0.90059)/2 - (0.43029-0.98596)/2
Ef_NH3 = Ef + Ef_change_up*1.602176634E-19

for Ti in T:
    Nc_i = 2*(2*math.pi*me_eff*me*k*Ti/(h**2))**1.5*1E-6
    n0_i = Nc_i*np.exp(-(Ec-Ef)/k/Ti)      ### the electron concentration, unit: cm-3
    print('The temperature (K) is: ', Ti)
    print('The thermal-equilibrium electron concentration (cm-3) is :', n0_i)
    n0_2D_i = n0_i*25*1E-8    ### convert the cm-3 to cm-2
    print('The thermal-equilibrium electron concentration (cm-2) is :')
    print('{:.5E}'.format(n0_2D_i))
    Nv_i = 2*(2*math.pi*mp_eff*me*k*Ti/(h**2))**1.5*1E-6
    p0_i = Nv_i*np.exp(-(Ef-Ev)/k/Ti)   ### the hole cocentration, unit: cm-3
    print('The thermal-equilibrium hole concentration (cm-3) is :', p0_i)

    ### calculate the carrier concentration for NO2@MoS2 based on the change of fermi level ###
    n0_NO2_i = Nc_i*np.exp(-(Ec-Ef_NO2)/k/Ti)
    print('The electron concentration (cm-3) of NO2@MoS2 is: ', n0_NO2_i)
    p0_NO2_i = Nv_i*np.exp(-(Ef_NO2-Ev)/k/Ti)
    print('The hole concentration (cm-3) of NO2@MoS2 is: ', p0_NO2_i)

    ### calculate the carrier concentration for NH3@MoS2 based on the change of fermi level ###
    n0_NH3_i = Nc_i*np.exp(-(Ec-Ef_NH3)/k/Ti)
    print('The electron concentration (cm-3) of NH3@MoS2 is: ', n0_NH3_i)
    p0_NH3_i = Nv_i*np.exp(-(Ef_NH3-Ev)/k/Ti)
    print('The hole concentration (cm-3) of NH3@MoS2 is: ', p0_NH3_i)
