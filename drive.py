

# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 15:11:22 2017

@author: jlevine7
"""

import numpy as np
import matplotlib.pyplot as plt


# Parameters
rho   = 0.01                                              
ten   = 40                                               
c     = np.sqrt(ten/rho)                                  



#step sizes and such
tf = 0.08
tsteps = int(tf*10000)
dt = tf/tsteps

xf = 1
xsteps = int(xf*101)
dx = xf/xsteps
x1d = np.linspace(-.5*xf,.5*xf,xsteps)

yf = 1
ysteps = int(yf*101)
dy = yf/ysteps
y1d = np.linspace(-.5*yf,.5*yf,ysteps)

c1    = dx/dt                                                 
ratio =  c*c/(c1*c1) #could be dx or dy as long as grid spacing is const. CFL criterium, set to 1 for stability




#initial conditions
x0 , y0 = .1 , .2 #point of initial disturbance
A = 0 #amplitude 
sigma = .05 #FWHM. bigger is wider

x2d, y2d = np.meshgrid(x1d,y1d)
r  = np.sqrt((x1d-x0)**2+(y2d-y0)**2) #this is the distance of every point from the center of displacement
initial = A * np.exp(-r**2/(2*sigma**2)) 

#drive peramaters   
driveamplitude = .1
omegadrive = 28 * (2*np.pi)
drivephase = 0
xdrive , ydrive = 50 , 50 #must be in terms of grid position 
#nxdrive , nydrive = xsteps*xdrive , ysteps*ydrive
driveon = y



#PDE solver
sol = np.zeros((tsteps,xsteps,ysteps))  #solution array

for xi in range (0,xsteps):
    for yi in range (0,ysteps):         #sets solution array to initial cond.
        sol[0,xi,yi],sol[1,xi,yi] =initial[xi,yi],initial[xi,yi] #the displacement before pluck = displacement at pluck
        

for ti in range(1,tsteps-1):
    for xi in range (1,xsteps-1):
        for yi in range (1,ysteps-1):
            if xi == xdrive and yi == ydrive and driveon = y:
                sol[ti+1,xi,yi] = driveamplitude * np.sin((omegadrive * ti*dt) + drivephase ) 
            else:
                sol[ti+1,xi,yi] = ratio * (sol[ti,xi+1,yi] -2*sol[ti,xi,yi]+sol[ti,xi-1,yi]
                +sol[ti,xi,yi+1]-2*sol[ti,xi,yi]+sol[ti,xi,yi-1])+2*sol[ti,xi,yi]-sol[ti-1,xi,yi]

trmsf = tsteps
trmsi = 0

solsquare = sol**2
solmean = solsquare.mean(axis=0)
solsqrt = solmean**(1/2)
plt.imshow(solsqrt)





def rms(ti,tf,xa,xb,ya,yb):
    som = sol
    rms = np.sqrt(np.mean((som)**2)) #calculates and plots rms difference between the homebrew function and scipy's
    rmsarr[m,l] = rms

"""
#needs fixing           
def rms(xa,xb,ya,yb):
    sol.size
    np.mean(y**2)
    np.sqrt(np.mean(y**2))
    sol[:,(xa,xb),(ya,yb)]
    return     
"""
            
#plotting
timeresolution = 10 #must be int. 5-10 is good. 1 for testing
timeslices = (tsteps * timeresolution) // 100 
plottype = "pix" 
timestepplot = tf//timeslices
minvalue = sol.min() 
maxvalue = sol.max() 
times = np.linspace(0,tsteps-1,timeslices)
plt.close('all')
figindex = 0



if plottype == 'pix':
    for ti in range(timeslices):
        t = times[ti]
        plt.close('all')
        plt.imshow(sol[t,:,:],vmin=minvalue, vmax=maxvalue)
        plt.colorbar()
        name1 = 'pix/ztemp'
        number = str(t)
        name2 = '.png'
        timeequalstring = 't= '
        timestepnumbtemp = str(t*dt)
        secondsstring = 'seconds'
        plt.title(timeequalstring+timestepnumbtemp+secondsstring)
        filename = name1+number+name2
        plt.savefig(filename, dpi=100)
else:
    for ti in range(timeslices):
        t = times[ti]
        plt.figure(figindex)            
        plt.pcolor(sol[t,:,:],vmin=minvalue, vmax=maxvalue)
        timeequalstring = 't= '
        timestepnumbtemp = str(t*dt)
        secondsstring = 'seconds'
        plottitle = timeequalstring+timestepnumbtemp+secondsstring
        plt.title(plottitle)
        plt.colorbar()
        figindex += 1












































