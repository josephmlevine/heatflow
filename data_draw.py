#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 11:57:33 2018

@author: steven
"""

import serial as s
import pylab as py

port = ('/dev/ttyUSB0')
num  = 2
ser  = s.Serial(port, 57600, timeout= 2)
tem  = py.zeros((10,10,num))
junk = ser.readline()
for n in range(num):
        
        line = ser.readline()
    
        # split the line into a list of floating-point values
        temperatures = [float(T) for T in line.split()]
