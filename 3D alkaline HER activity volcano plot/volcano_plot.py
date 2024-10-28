import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from pylab import *


Gh = [-0.30, -0.20, -0.20, -0.16, -0.13]  # binding energy of H*, eV. The data can be obtained from the reference 'Nature Energy, 2020, 5(11): 891-899.'
Goh = [0.69, 0.21, -0.43, -0.49, -0.08]   # binding energy of OH*, eV. The data can be obtained from the reference 'Nature Energy, 2020, 5(11): 891-899.'
Element = ['Pt(111)', 'Pt(553)', 'Pt(553) Ru*', 'Pt(111) Ru*', 'PtRu(111)']


xmin=-1.0
xmax=1.5001
ymin=-1.0
ymax=1.5001
x = np.arange(xmin, xmax, 0.01) # GH
y = np.arange(ymin, ymax, 0.01)  # GOH
act = np.zeros((len(y),len(x)))
kb=8.617E-5 # boltzmann constant, eV/K
T = 300 # temperature, K
A=3E+10 # pre-exponential factor s-1
RT = kb*T

def cal_rate1(gh, goh):
    dG_ts = 0.72*gh + 0.51*goh + 0.38
    rate1 = A*exp(-dG_ts/RT)*1/(1+exp(-gh/RT))
    rate1_log = math.log(rate1)
    return rate1_log

def cal_rate2(gh, goh):
    rate2 = A*exp(goh/RT)
    rate2_log = math.log(rate2)
    return rate2_log

def cal_rate3(gh, goh):
    rate3 = A*exp(gh/RT)
    rate3_log = math.log(rate3)
    return rate3_log


X, Y = np.meshgrid(x, y)
for j in range(len(x)):
    for i in range(len(y)):
        rate1 = cal_rate1(X[i][j], Y[i][j])
        rate2 = cal_rate2(X[i][j], Y[i][j])
        rate3 = cal_rate3(X[i][j], Y[i][j])
        act[i][j] = min(rate1, rate2, rate3)

act_xyZ = [X, Y, act]

Xt = act_xyZ[0]
Yt = act_xyZ[1]
E  = act_xyZ[2]
plt.figure()
plt.rc('font', family='Arial')
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
CS0= plt.contour(Xt, Yt, E, 30,
        colors = 'k',   linewidths= 0.1)
#CS4 = plt.contourf(Xt, Yt, E, 20,
        #origin = 'lower',
        #extend = 'both',
        #cmap = 'jet_r')

CS4 = plt.contourf(Xt, Yt, E, 20,
        origin = 'lower',
        extend = 'both')
x0=[];y0=[]
cbar =plt.colorbar(CS4,shrink=0.8)

for i in range(len(Gh)):
    g_h = Gh[i]
    g_oh = Goh[i]

    plt.plot([g_h],[g_oh],'k^',mfc='w',markersize=8)
    
    if Element[i] in ['Pt(111) Ru*']:
        plt.text(g_h+0.04,g_oh-0.01,Element[i], fontsize=9, color = "k", style = "normal", weight = "light", verticalalignment='center',)
    else:
        plt.text(g_h+0.03,g_oh+0.03,Element[i], fontsize=9, color = "k", style = "normal", weight = "light", verticalalignment='center',)


plt.xlabel(r'$\Delta \rm G_{\rm {H*}}$ (eV)', fontsize=10)
plt.ylabel(r'$\Delta \rm G_{\rm {OH*}}$ (eV)', fontsize=10)
plt.savefig('her_act_nature_energy.eps')
plt.show()
