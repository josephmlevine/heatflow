#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 12:19:18 2018

@author: steven
"""

import pylab as py
from scipy import interpolate as intp

# bring in data
dataN = 'matches'
data = py.load(dataN +'.npy')
st_frame = 0 #input the initial starting frame (0)
initial_data = data[:,:,st_frame]

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

#setup to translate experimental grid to simulation grid
x = py.arange(0,10,1)
y = x
xx, yy = py.meshgrid(x, y)
f = intp.interp2d(x, y, initial_data, kind='cubic')
#Now apply interp2d to the bigger simulation grid
xnew = py.arange(0, 100, 1)
ynew = xnew
new_vals = f(xnew, ynew)

#main loop to calculate new temperatures:
run_time = 0.01#seconds
Nx,Ny,Nt = 100,100, int(run_time/time_step)
Temperature = py.zeros((Nx,Ny,Nt))
Temperature[:,:,0] = new_vals
minT = Temperature.min()
maxT = Temperature.max()

'''
for i in range(1,Nx-1):
    for j in range(1,Ny-1):
        for n in range(Nt-1):
            Temperature[i,j,n+1]= (1-4*Fom) * Temperature[i,j,n] + Fom*(Temperature[i+1,j,n] + Temperature[i,j+1,n] + Temperature[i-1,j,n] + Temperature[i,j-1,n])
'''

time = py.arange(0,120,time_step)
print_time = py.arange(0,120,2)

itt_number = 0
outputName = 'simulation'

exper_grid = py.zeros([N,N])

for t in range(len(time)):
	filename = outputName + "%04d.txt" %itt_number

#now we are looking at how the temperature changes for each cell.

	#this will check to see if the time for the simulation
	#corresponds to the time that the experiment will read
	#these are the only files we want
	for p in print_time:

		if( time[t] > p - time_step/2. and time[t] < p + time_step/2. ):

			for i in range(N):

				for j in range(N):

					exper_grid[i,j] = Temperature[i*N + N/2,j*N + N/2]
                    
			figname = dataN+ "_" + '%04d' %(time[t]) + ".png"
			py.contourf(exper_grid,py.linspace(minT,maxT,N))
			py.colorbar()
			py.savefig(figname)
			py.clf()

			print ("Picture " + figname + " saved")

	#takes into account the adiabatic condition for the boundary 
	for i in range(NN):

		for j in range(NN):

			if i == 99 and j == 99:

				Temp_new[i,j] = (1 - 4*Fom)*Temperature[i,j] + Fom*(Temperature[i-1,j] + Temperature[i,j-1] + 2*Temperature[i,j])

			elif i == 99:

				Temp_new[i,j] = (1 - 4*Fom)*Temperature[i,j] + Fom*(Temperature[i,j] + Temperature[i-1,j] + Temperature[i,j+1] + Temperature[i,j-1])
				
			elif j == 99:

				Temp_new[i,j] = (1 - 4*Fom)*Temperature[i,j] + Fom*(Temperature[i+1,j] + Temperature[i-1,j] + Temperature[i,j] + Temperature[i,j-1])
				
			elif i == 0:

				Temp_new[i,j] = (1 - 4*Fom)*Temperature[i,j] + Fom*(Temperature[i+1,j] + Temperature[i,j] + Temperature[i,j+1] + Temperature[i,j-1])

			elif j == 0:

				Temp_new[i,j] = (1 - 4*Fom)*Temperature[i,j] + Fom*(Temperature[i+1,j] + Temperature[i-1,j] + Temperature[i,j+1] + Temperature[i,j])
	
			elif j == 0 and i == 0:

				Temp_new[i,j] = (1 - 4*Fom)*Temperature[i,j] + Fom*(Temperature[i+1,j] + Temperature[i,j+1] + 2*Temperature[i,j])
		
			elif i == 0 and j == 99:

				Temp_new[i,j] = (1 - 4*Fom)*Temperature[i,j] + Fom*(Temperature[i+1,j] + Temperature[i,j-1] + 2*Temperature[i,j])

			elif i == 99 and j == 0:


				Temp_new[i,j] = (1 - 4*Fom)*Temperature[i,j] + Fom*(Temperature[i-1,j] + Temperature[i,j+1] + 2*Temperature[i,j])
			else:

				Temp_new[i,j] = (1 - 4*Fom)*Temperature[i,j] + Fom*(Temperature[i+1,j] + Temperature[i-1,j] + Temperature[i,j+1] + Temperature[i,j-1])

	#copy the grid, rinse and repeat
	Temperature = copy(Temp_new)
	itt_number += 1