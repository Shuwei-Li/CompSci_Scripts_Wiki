#!/usr/local/bin/pythona
######################################################################
#This is the script to calculate the activity of H2+O2--H2O over Pt-111 surface

#### elementary reaction
#1. 1/2H2(g) + * -- H*
#2. O2(g) + * -- O2*
#3. O2* + * -- 2O*
#4. O2* + H* -- OOH* + *          rate-determining step
#5. OOH* + H* -- 2OH*
#6. OH* + H* -- H2O* + *
#7. H2O* -- H2O(g) + *
#The following parameters are needed 
# dG1, dG2, dG3, dG4, dG5, dG6, dG7; dEa4; T; ph2, po2, ph2o
######################################################################

import numpy as np


T = 300 # K
ph2 = 0.04 # bar
po2 = 0.21 # bar
ph2o = 0.035 # bar

h=4.136E-15 # planck constant, eV*s
kb=8.617E-5 # boltzmann constant, eV/K
mu=kb*T/h

## Strain = +3.0%
#### import your dG
dG1 = -0.28
dG2 = -0.16
dG3 = -1.58
dG4 = -0.07
dG5 = -2.33
dG6 = -0.4
dG7 = -0.48


### import your Ea3
Ea4 = 0.45

# calculate equilibrium constant
K1 = np.exp(-dG1/kb/T)
K2 = np.exp(-dG2/kb/T)
K3 = np.exp(-dG3/kb/T)
K4 = np.exp(-dG4/kb/T)
K5 = np.exp(-dG5/kb/T)
K6 = np.exp(-dG6/kb/T)
K7 = np.exp(-dG7/kb/T)

# calculate the coverage degree of sites, cov_site

cov_site = 1/(1 + K1*ph2**0.5 + K2*po2 + (K2*K3*po2)**0.5 + ph2o/K7 + ph2o/(K1*K6*K7*ph2**0.5) + ph2o**2/(K1**3*K5*K6**2*K7**2*ph2**1.5))


# calculate the coverage degree of H, O2, OH, OOH

cov_H = K1*ph2**0.5*cov_site
cov_O2 = K2*po2*cov_site
cov_O = (K2*K3*po2)**0.5*cov_site
cov_OH = ph2o*cov_site/(K1*K6*K7*ph2**0.5)
cov_OOH = ph2o**2*cov_site/(K1**3*K5*K6**2*K7**2*ph2**1.5)
cov_H2O = ph2o*cov_site/K7






# calculate the forward rate constant for reaction3-----rate-determining step

k4 = mu*np.exp(-Ea4/kb/T)



# calculate the y

Keq = K1**4*K2*K4*K5*K6**2*K7**2
y = ph2o**2/Keq/ph2**2/po2


# calculate the rate

rate = k4*K1*K2*ph2**0.5*po2*cov_site**2*(1-y)

print('k4    ', k4)
print('K1    ', K1)
print('K2    ', K2)
print('K3    ', K3)
print('K4    ', K4)
print('K5    ', K5)
print('K6    ', K6)
print('K7    ', K7)

print('Keq    ', Keq)
 
print('cov_site    ', cov_site)
print('cov_H    ', cov_H)
print('cov_O2    ', cov_O2)
print('cov_O    ', cov_O)
print('cov_OH    ', cov_OH)
print('cov_OOH    ', cov_OOH)
print('cov_H2O    ', cov_H2O)

print('y    ', y)

print('rate    ', rate)


