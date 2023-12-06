#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 16:00:28 2023

@author: poletti
"""

import matplotlib.pyplot as plt
import numpy as np
import shutil
import pandas as pd

#%%
# Experimental forces estimated from the digitalization of 
# "Recent progress in flapping wing aerodynamics and aeroelasticity" from Shyy 2010

t_Ct_Shyy_ref = [0.000746269, 0.014925373, 0.032835821, 0.048507463, 0.06119403 , 0.070895522, 0.082089552, 0.091044776, 0.102238806, 0.111940299, 0.123880597, 0.132089552, 0.13880597 , 0.15       , 0.158955224, 0.170895522, 0.187313433, 0.205223881, 0.220895522, 0.23880597 , 0.262686567, 0.280597015, 0.292537313, 0.308955224, 0.319402985, 0.330597015, 0.338059701, 0.35       , 0.360447761, 0.370149254, 0.37761194 , 0.389552239, 0.397761194, 0.408955224, 0.420149254, 0.429104478, 0.439552239, 0.45       , 0.462686567, 0.473134328, 0.490298507, 0.513432836, 0.532835821, 0.552985075, 0.567910448, 0.580597015, 0.591791045, 0.602238806, 0.609701493, 0.617910448, 0.62761194 , 0.63880597 , 0.650746269, 0.661940299, 0.671641791, 0.68358209 , 0.697014925, 0.711940299, 0.731343284, 0.752238806, 0.774626866, 0.787313433, 0.800746269, 0.810447761, 0.821641791, 0.832835821, 0.83880597 , 0.848507463, 0.858955224, 0.870895522, 0.882835821, 0.889552239, 0.9        , 0.908955224, 0.920895522, 0.932089552, 0.944776119, 0.958208955, 0.976119403, 0.997761194 ]
Ct_Shyy_ref   = [-0.196502058, -0.181069959, -0.114197531, -0.016460905, 0.08127572  , 0.158436214 , 0.281893004 , 0.395061728 , 0.523662551 , 0.652263374 , 0.801440329 , 0.899176955 , 0.991769547 , 1.117798354 , 1.228395062 , 1.351851852 , 1.495884774 , 1.619341564 , 1.686213992 , 1.727366255 , 1.701646091 , 1.624485597 , 1.552469136 , 1.439300412 , 1.356995885 , 1.24382716  , 1.176954733 , 1.063786008 , 0.945473251 , 0.83744856  , 0.75        , 0.621399177 , 0.528806584 , 0.405349794 , 0.302469136 , 0.209876543 , 0.117283951 , 0.034979424 , -0.052469136, -0.109053498, -0.170781893, -0.186213992, -0.134773663, -0.026748971, 0.091563786 , 0.209876543 , 0.343621399 , 0.446502058 , 0.554526749 , 0.672839506 , 0.796296296 , 0.935185185 , 1.094650206 , 1.248971193 , 1.367283951 , 1.501028807 , 1.624485597 , 1.742798354 , 1.830246914 , 1.845679012 , 1.78909465  , 1.722222222 , 1.619341564 , 1.537037037 , 1.434156379 , 1.326131687 , 1.233539095 , 1.125514403 , 1.002057613 , 0.868312757 , 0.71399177  , 0.62654321  , 0.508230453 , 0.400205761 , 0.281893004 , 0.158436214 , 0.040123457 , -0.04218107 , -0.145061728, -0.191358025 ]

f = 0.96

#%%
# Copy the files from the postProcess folder

name_in = 'coefficient_0.dat'
path_out = 'RESULTS/'

# Source and destination file paths
source_file = './postProcessing/airFoil_coef/0/' + name_in
destination_file_airfoil = path_out + 'coefficient_a.dat'

shutil.copy(source_file, destination_file_airfoil)

source_file = './postProcessing/flatPlate_coef/0/' + name_in
destination_file_plate = path_out + 'coefficient_p.dat'

shutil.copy(source_file, destination_file_plate)

#%% 
# Read the file

# Might not be needed depending on the OF version used
def remove_parentheses(s):
    return float(s.replace('(', '').replace(')', ''))


df_airfoil = pd.read_csv(destination_file_airfoil, delimiter='\s+', skiprows=13, header=None,skipinitialspace=True,
                 names=['t', 'cx', 'cy', 'cz', 'croll', 'cpitch', 'cyaw', 'cxf', 'cxr', 'cyf', 'cyr','czf','czr'],
                 converters={float(i): remove_parentheses for i in range(13)})

df_plate = pd.read_csv(destination_file_plate, delimiter='\s+', skiprows=13, header=None,skipinitialspace=True,
                 names=['t', 'cx', 'cy', 'cz', 'croll', 'cpitch', 'cyaw', 'cxf', 'cxr', 'cyf', 'cyr','czf','czr'],
                 converters={float(i): remove_parentheses for i in range(13)})

#%%
# Few post-processing 

# Convert into numpy
t = np.array(df_airfoil['t'])
total_thurst = -(np.array(df_airfoil['cx']) + np.array(df_plate['cx']))


# Repeat the experiments (t/T in [0,1]) for the duration of the simulation (t/T in [0,x]) 
nShyy = int(np.floor(t[-1]))
t_Ct_Shyy = []; Ct_Shyy = []
for i in range(nShyy):
    t_tmp = [j+i for j in t_Ct_Shyy_ref]
    t_Ct_Shyy += t_tmp
    Ct_Shyy   += Ct_Shyy_ref

#%%
# Plot the files
NundeSample    = 5
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=22)
plt.rc('text.latex', preamble=r'\usepackage{amsmath,bm}')
plt.close('all')
plt.figure()
plt.plot(t_Ct_Shyy, Ct_Shyy, 'k', label="Experiments")  
plt.plot((t*f)[::NundeSample], total_thurst[::NundeSample], 'b', label="OpenFOAM")  
plt.ylabel("Ct [-]")
plt.xlabel("t/T [-]")
plt.tight_layout()
plt.minorticks_on(); 
plt.grid(which='both',axis='both',alpha=0.5); 
plt.legend()
plt.xticks(np.arange(0, 6, step=1))
plt.xlim(0,np.max(t*f));
plt.yticks(np.arange(-0.5, 2, step=0.5))
plt.ylim(-0.3,2);

#%%
# TE displacement
# Flexible Flapping Airfoil Propulsion at Low Reynolds Numbers, S. Heathcote 2007
ha = 0.0175
t_te_heat_ref = [0, 0.041851107, 0.081609658, 0.121368209, 0.16112676 , 0.205070423, 0.242736418, 0.28249497, 0.326438632, 0.364104628, 0.403863179, 0.447806841, 0.487565392, 0.527323944, 0.56917505 , 0.611026157, 0.652877264, 0.690543259, 0.730301811, 0.767967807, 0.809818913, 0.85167002 , 0.891428571, 0.933279678, 0.968853118, 1.016981891, 1.052555332, 1.096498994, 1.13416499 , 1.171830986, 1.215774648, 1.255533199, 1.297384305, 1.337142857, 1.378993963, 1.418752515, 1.458511066, 1.502454728, 1.542213279, 1.581971831, 1.621730382, 1.661488934, 1.70334004 , 1.743098592, 1.782857143, 1.82470825 , 1.866559356, 1.904225352, 1.946076459, 1.98583501 , 2.023501006, 2.067444668, 2.105110664]
te_heat_ref   = [0.016421039, 0.020924067, 0.0232056, 0.025186932, 0.026327699, 0.02476665, 0.02140439, 0.017321645, 0.012638496, 0.006694501, 0.000330222, -0.005913976, -0.012038093, -0.017681887, -0.021644551, -0.025126892, -0.026988143, -0.027288345, -0.025367053, -0.021344349, -0.016661201, -0.011197528, -0.005193491, 0.001230827, 0.007114783, 0.012878658, 0.017681887, 0.021704591, 0.023806004, 0.025427094, 0.026087538, 0.023625883, 0.020023461, 0.015880676, 0.010717205, 0.004653128, -0.00171115, -0.008015389, -0.013719223, -0.018642533, -0.022905399, -0.025367053, -0.027048184, -0.026868063, -0.023926085, -0.0197833, -0.01492003, -0.009396317, -0.003092079, 0.003032038, 0.008795913, 0.014259586, 0.018762614]

df_le = pd.read_csv(path_out + 'LE_disp.csv', delimiter=',', skiprows=1, header=None,skipinitialspace=True,
                 names=['t', 'le'])
df_te = pd.read_csv(path_out + 'TE_disp.csv', delimiter=',', skiprows=1, header=None,skipinitialspace=True,
                 names=['t', 'te'])

t_tmp = np.array(df_te['t'])
nt= int(np.floor(t_tmp[-1]))
t_te_heat = []; te_heat = []
for i in range(0,nt,2):
    t_tmp = [j*f+i for j in t_te_heat_ref]
    t_te_heat += t_tmp
    te_heat   += te_heat_ref
    
t_te_heat = np.array(t_te_heat)  
te_heat   = np.array(te_heat)
t_te = np.array(df_te['t'])
te = np.array(df_te['te'])+ha

t_le = np.array(df_le['t'])
le = np.array(df_le['le'])+ha

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=22)
plt.rc('text.latex', preamble=r'\usepackage{amsmath,bm}')
plt.close('all')

plt.figure()
plt.plot(t_le*f, le/ha, 'r', label="TE-OpenFOAM") 
plt.plot(t_te*f, te/ha, 'b', label="TE-OpenFOAM") 
plt.plot(t_te_heat, te_heat/ha, 'k', label="TE-Experiments") 
plt.ylabel("Displacement [-]")
plt.xlabel("t/T [-]")
plt.tight_layout()
plt.minorticks_on(); 
plt.grid(which='both',axis='both',alpha=0.5); 
plt.legend()
plt.xticks(np.arange(0, 6, step=1))
plt.xlim(0,np.max(t*f));
plt.yticks(np.arange(-2, 2, step=0.5))
plt.ylim(-1.7,1.7);
