#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 16:37:06 2018

@author: herbiek9
"""
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
from scipy import interpolate as intp



time_step = 0.001
N = 100
#constants:
rho_al = 2.7 #g/cm^3
k_al = 2.05 #W/cm*K
cp_al = .9 #J/g*K
dx = 2.95 #cm
dxN = dx/N
dx2 = dxN**2

#Find the time_steps need for PDE Solver
alpha_al = k_al/(cp_al*rho_al)
tFom = dx2/alpha_al
Fom = 1/515
time_step = tFom*Fom


data = np.load('Data/blowtorch_center_w_foam.npy')

init_con_raw = data[:,:,4]

#setup to translate experimental grid to simulation grid

x = py.arange(0,10,1)
y = x
xx, yy = py.meshgrid(x, y)
f = intp.interp2d(x, y, init_con_raw, kind='cubic')



#Now apply interp2d to the bigger simulation grid
xnew = np.linspace(0, 10, 100)
ynew = np.linspace(0, 10, 100)
new_vals = f(xnew, ynew)


run_time = 0.001#seconds
Nx,Ny,Nt = 100,100, int(run_time/time_step)
Temperature = py.zeros((Nx,Ny,Nt))
Temperature[:,:,0] = new_vals


for i in range(1,Nx-1):
    for j in range(1,Ny-1):
        for n in range(Nt-1):
            Temperature[i,j,n+1]= (1-4*Fom) * Temperature[i,j,n] + \
            Fom*(Temperature[i+1,j,n] + Temperature[i,j+1,n] + Temperature[i-1,j,n] \
                 + Temperature[i,j-1,n])
            
'''           
for ti in range(1,Nt-1):
    for xi in range (1,Nx-1):
        for yi in range (1,Ny-1):
            sol[ti+1,xi,yi] = ratio * (sol[ti,xi+1,yi] -2*sol[ti,xi,yi]+sol[ti,xi-1,yi]
            +sol[ti,xi,yi+1]-2*sol[ti,xi,yi]+sol[ti,xi,yi-1])+2*sol[ti,xi,yi]-sol[ti-1,xi,yi]
'''          
            