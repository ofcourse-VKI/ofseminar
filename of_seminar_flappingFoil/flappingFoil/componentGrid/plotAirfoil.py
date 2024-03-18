#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 23:10:24 2023

@author: poletti
"""

import numpy as np 
import matplotlib.pyplot as plt

# On-off flag
of_write = 0 # 1, write edges in txt files read by blockMeshDict

x0 = 0                   # [mm] airfoil LE position
c  = 30                  # [mm] chord        
T  = 10                  # [mm] thickness    
l  = 60                  # [mm] plate length 
N  = 1000                # [-]  number of points of the airfoil 
x_airfoil = 4.784764     # [mm] arbitrary junction point separating nose part and tail part of the airfoil
alpha_te  = 1.149        # [rad] angle between the horizontal and the normal taken at the airoil TE
d = c / np.sin(alpha_te) # [mm]  distance from the component grid (CG) boundary to the airfoil TE, taken normally at the airfoil TE


# NACA0012 defined in:
# http://airfoiltools.com/airfoil/naca4digit?MNaca4DigitForm%5Bcamber%5D
a0 =  0.2969; a1 = -0.126
a2 = -0.3516; a3 =  0.2843
a4 = -0.1036; T_2 = (T/c)/0.2
x  = np.linspace(x0,c,N) # x-coord of the airfoil
y = (a0*(x/c)**0.5 + a1*(x/c) + a2*(x/c)**2 + a3*(x/c)**3 + a4*(x/c)**4) 
yp = T_2*y               # y-coord of the airfoil
ym = -yp                 # y-coord of the airfoil

# Computing the derivative of the NACA profile
Nnormal = 1000
x2 = np.linspace(x0,c,Nnormal) # x-coord of the airfoil
y2 = T_2*(a0*(x2/c)**0.5 + a1*(x2/c) + a2*(x2/c)**2 + a3*(x2/c)**3 + a4*(x2/c)**4) # NACA0012 profile
yprime = T_2*(0.5*a0*((x2/c)**(-0.5))*(1/c) + a1/c + 2*a2*(x2/c)*(1/c) + 3*a3*((x2/c)**2)*(1/c) + 4*a4*((x2/c)**3)*(1/c)) # derivative of the profile
n = -1/(yprime) # normals of the profile computed for the Nnormal points

Npoints = 5 # only used to plot the normal lines
# init
xlinesp = []; ylinesp = []
xlinesm = []; ylinesm = []
x_list_nose = []; y_list_nose = []
x_list_te   = []; y_list_te   = []

xNACA_list_nose = []; yNACA_list_nose = []
xNACA_list_te   = []; yNACA_list_te   = []
# Loop on all the points of the profile and compute their normal lines
for i in range(len(x2)): 
    alpha = np.arctan(n[i]*c/1000) # angle made between the normal and the horrizontal
    # y>0
    x0_i = x2[i] + d*np.cos(np.abs(alpha))*np.sign(alpha-1e-6)   # x-coord of the component grid boundary
    y0_i = y2[i]*c + d*np.sin(np.abs(alpha))                     # y-coord of the component grid boundary
    xn = np.linspace(x0_i,x2[i] ,Npoints)                        # xcoord of the normal lines
    yn = np.linspace(y0_i,y2[i]*c  ,Npoints)                     # ycoord of the normal lines 
    xlinesp.append(xn)
    ylinesp.append(yn)
    # y<0
    ylinesm.append(-yn)
    
    # Used to write edges of the component boundary around the airfoil
    if x2[i]<x_airfoil:
        x_list_nose.append(x0_i)
        y_list_nose.append(y0_i)
        xNACA_list_nose.append(x2[i])
        yNACA_list_nose.append(y2[i]*c)
    else:
        x_list_te.append(x0_i)
        y_list_te.append(y0_i)
        xNACA_list_te.append(x2[i])
        yNACA_list_te.append(y2[i]*c)
        
#%% 
# Plot the airfoil & component grid boundary
plt.close('all')
plt.figure()
plt.plot(x,yp*c,'b')                               # airfoil
plt.plot(x,ym*c,'b')                               # airfoil
plt.plot(np.linspace(c,c+l,100),np.zeros(100),'b') # plate
for j in range(len(xlinesp)):
    plt.plot(xlinesp[j],ylinesp[j],'k')
    plt.plot(xlinesp[j],ylinesm[j],'k')
plt.scatter(x2,y2*c,c='k')
plt.scatter(x2,-y2*c,c='k')


#%% 
# Write txt file for the edges of the blockMeshDict

if of_write:
    path = 'system/edges/'
    
    #%% 
    # Component grid boundaries around the airoil
    
    ## UP PART
    ## NOSE PART
    f= open(path + "nose_z0_ypos.txt","w+")
    f.write("polyLine 31 29\n")
    f.write("( \n")
    for i in range(len(x_list_nose)-1):
         f.write("( %.8f %.8f $z0)\n" % (x_list_nose[i],y_list_nose[i]))
    f.write(") \r")   
    f.close();
    
    f= open(path + "nose_z1_ypos.txt","w+")
    f.write("polyLine 32 30\n")
    f.write("( \n")
    for i in range(len(x_list_nose)-1):
         f.write("( %.8f %.8f $z1)\n" % (x_list_nose[i],y_list_nose[i]))
    f.write(") \r")   
    f.close();
    
    ## TAIL PART
    f= open(path + "tail_z0_ypos.txt","w+")
    f.write("polyLine 29 5\n")
    f.write("( \n")
    for i in range(len(x_list_te)-1):
         f.write("( %.8f %.8f $z0)\n" % (x_list_te[i],y_list_te[i]))
    f.write(") \r")   
    f.close();
    
    f= open(path + "tail_z1_ypos.txt","w+")
    f.write("polyLine 30 14\n")
    f.write("( \n")
    for i in range(len(x_list_te)-1):
         f.write("( %.8f %.8f $z1)\n" % (x_list_te[i],y_list_te[i]))
    f.write(") \r")   
    f.close();
    
    ## DOWN PART
    ## NOSE PART
    f= open(path + "nose_z0_yneg.txt","w+")
    f.write("polyLine 31 18\n")
    f.write("( \n")
    for i in range(len(x_list_nose)-1):
         f.write("( %.8f %.8f $z0)\n" % (x_list_nose[i],-y_list_nose[i]))
    f.write(") \r")   
    f.close();
    
    f= open(path + "nose_z1_yneg.txt","w+")
    f.write("polyLine 32 28\n")
    f.write("( \n")
    for i in range(len(x_list_nose)-1):
         f.write("( %.8f %.8f $z1)\n" % (x_list_nose[i],-y_list_nose[i]))
    f.write(") \r")   
    f.close();
    
    ## TAIL PART
    f= open(path + "tail_z0_yneg.txt","w+")
    f.write("polyLine 18 19\n")
    f.write("( \n")
    for i in range(len(x_list_te)-1):
         f.write("( %.8f %.8f $z0)\n" % (x_list_te[i],-y_list_te[i]))
    f.write(") \r")   
    f.close();
    
    f= open(path + "tail_z1_yneg.txt","w+")
    f.write("polyLine 28 27\n")
    f.write("( \n")
    for i in range(len(x_list_te)-1):
         f.write("( %.8f %.8f $z1)\n" % (x_list_te[i],-y_list_te[i]))
    f.write(") \r")   
    f.close();
    
    #%% 
    # NACA profile
    
    ## UP PART
    ## NOSE PART
    f= open(path + "noseNACA_z0_ypos.txt","w+")
    f.write("spline 7 6\n")
    f.write("( \n")
    for i in range(len(xNACA_list_nose)-1):
         f.write("( %.8f %.8f $z0)\n" % (xNACA_list_nose[i],yNACA_list_nose[i]))
    f.write(") \r")   
    f.close();
    
    f= open(path + "noseNACA_z1_ypos.txt","w+")
    f.write("spline 16 15\n")
    f.write("( \n")
    for i in range(len(xNACA_list_nose)-1):
         f.write("( %.8f %.8f $z1)\n" % (xNACA_list_nose[i],yNACA_list_nose[i]))
    f.write(") \r")   
    f.close();
    
    ## TAIL PART
    f= open(path + "tailNACA_z0_ypos.txt","w+")
    f.write("spline 6 0\n")
    f.write("( \n")
    for i in range(len(xNACA_list_te)-1):
         f.write("( %.8f %.8f $z0)\n" % (xNACA_list_te[i],yNACA_list_te[i]))
    f.write(") \r")   
    f.close();
    
    f= open(path + "tailNACA_z1_ypos.txt","w+")
    f.write("spline 15 9\n")
    f.write("( \n")
    for i in range(len(xNACA_list_te)-1):
         f.write("( %.8f %.8f $z1)\n" % (xNACA_list_te[i],yNACA_list_te[i]))
    f.write(") \r")   
    f.close();
    
    ## DOWN PART
    ## NOSE PART
    f= open(path + "noseNACA_z0_yneg.txt","w+")
    f.write("spline 7 8\n")
    f.write("( \n")
    for i in range(len(xNACA_list_nose)-1):
         f.write("( %.8f %.8f $z0)\n" % (xNACA_list_nose[i],-yNACA_list_nose[i]))
    f.write(") \r")   
    f.close();
    
    f= open(path + "noseNACA_z1_yneg.txt","w+")
    f.write("spline 16 17\n")
    f.write("( \n")
    for i in range(len(xNACA_list_nose)-1):
         f.write("( %.8f %.8f $z1)\n" % (xNACA_list_nose[i],-yNACA_list_nose[i]))
    f.write(") \r")   
    f.close();
    
    ## TAIL PART
    f= open(path + "tailNACA_z0_yneg.txt","w+")
    f.write("spline 8 0\n")
    f.write("( \n")
    for i in range(len(x_list_te)-1):
         f.write("( %.8f %.8f $z0)\n" % (xNACA_list_te[i],-yNACA_list_te[i]))
    f.write(") \r")   
    f.close();
    
    f= open(path + "tailNACA_z1_yneg.txt","w+")
    f.write("spline 17 9\n")
    f.write("( \n")
    for i in range(len(xNACA_list_te)-1):
         f.write("( %.8f %.8f $z1)\n" % (xNACA_list_te[i],-yNACA_list_te[i]))
    f.write(") \r")   
    f.close();