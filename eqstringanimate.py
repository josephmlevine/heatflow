""" From "COMPUTATIONAL PHYSICS" & "COMPUTER PROBLEMS in PHYSICS"
    by RH Landau, MJ Paez, and CC Bordeianu (deceased)
    Copyright R Landau, Oregon State Unv, MJ Paez, Univ Antioquia, 
    C Bordeianu, Univ Bucharest, 2017. 
    Please respect copyright & acknowledge our work."""
	 
# EqStringAnimate.py:    Animated leapfrog solution of wave equation

#from visual import *
import numpy as np
import matplotlib.pyplot as plt

# Parameters
rho   = 0.01                                              
ten   = 40.                                               
c     = np.sqrt(ten/rho)                                  


tf = 0.1
tsteps = 1000
dt = tf/tsteps

xf = 1
xsteps = 100
dx = xf/xsteps

c1    = dx/dt                                                 
ratio =  c*c/(c1*c1) # CFL criterium, set to < .5 for stability

pluckpoint=81
# Initialization
xi = np.zeros((xsteps), float)        
                       
for i in range(0, pluckpoint):     
    xi[i] = 0.00125*i;                      
for i in range (pluckpoint, xsteps-1):  
    xi[i] = 0.1 - 0.005*(i - 80)       #find the initial conditions  
           
       
dis = np.zeros((tsteps,xsteps))  #displacement array


for xin in range (0,xsteps):         #sets initial cond.
    dis[0,xin],dis[1,xin] =xi[xin],xi[xin] #the displacement before pluck = displacement at pluck
 
for tin in range(1,tsteps-1):
    for xin in range (1,xsteps-1):
        dis[tin+1,xin] = 2*dis[tin,xin] - dis[tin-1,xin] + ratio * (dis[tin,xin+1] + dis[tin,xin-1] - 2*dis[tin,xin])

plt.imshow(dis)



























'''
# Later time steps
for i in range(1, 100): 
    xi[i,1] = xi[i,0] + 0.5*ratio*(xi[i+1,0]+xi[i-1,0]-2*xi[i,0])
while 1:                               
#    rate(50)                                             # Plotting delay
    for i in range(1, 100):              
        xi[i,2] = 2.*xi[i,1] - xi[i,0] + ratio * (xi[i+1,1]+xi[i-1,1]-2*xi[i, 1])
#    for i in range(1, 100):
#         vibst.x[i] = 2.*i - 100.0                       # Scale for plot
#         vibst.y[i] = 300.*xi[i, 2]                             
#    vibst.pos                                                
    for i in range(0, 101):
        xi[i, 0] = xi[i, 1]                                
        xi[i, 1] = xi[i, 2]          
                   

print("Done!")
'''