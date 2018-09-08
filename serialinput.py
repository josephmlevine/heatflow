#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 12:15:46 2018

This script takes in data from the McDougall/Ayars 2d heat flow 
apparatus and generates a 10x10x(time steps) array of data that is exported
as a .nyp file for later plotting/analysis 
"""

import serial as pyser
import numpy as np
import matplotlib.pyplot as plt

#establish a connection with arduino. Port will probably need to be changed.
#try for linux try port = '/dev/ttyUSB1' for USB1 - USB4
port = '/dev/cu.usbserial-AH01872J'
ser = pyser.Serial(port, 57600, timeout=2)

#user inputs
time = .1 #in minutes
outfile = "samplefilename" #name your output file 


timesteps = int((time*60)/2 + 3)
print ("expect", timesteps, "frames")



array = np.zeros((10,10,timesteps))
   


        
k = 0
while k <= (timesteps - 1):
    try:
        data = ser.readline()
        temptemp = [float(T) for T in data.split()] #temporary temperature. yes really.
        
    
        for i in range(10):
            for j in range(10):
                array[i,j,k] = temptemp[i*10+j]
                
        k += 1
        print("got frame", k)
                
   
    #error stuff.        
    except ValueError:
       # error at float conversion
       print ("Conversion error: frame dropped, continuing collection.")
       continue
    except IndexError:
        # not all values received for a measurement set
        print ("Dont panic. Partial frame received: frame dropped, continuing collection.\
               If this continues, quit program, unplug and replug adruino")
        print('')
        print('')
        continue
    except pyser.serialutil.SerialException:
        # Someone just unplugged the device, or other loss of communication.
        print ("Device disconnected, exiting.")
        break
    except: 
        # Catch-all
        print ("Unknown error! Exiting program... Panic...")
        break

np.save(outfile, array)
            


