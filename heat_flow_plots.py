#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 21:08:57 2018

Imports array files (.npy) and plots each frame.

You may have to change the file directories, I was working
from my own folder.

@author: steven
"""
import pylab as py

data_1 = py.load('matches.npy') 
data_2 = py.load('icecube_lefton.npy')

plotting = True
if plotting == True:
    for i in range(len(data_1[0,0,:])):
        py.figure()
        sufix = str(i)
        path = 'data1_images/'
        name = 'data1'
        fname = path + name + sufix
        py.imshow(data_1[:,:,i])
        py.savefig(fname,format='png')
        py.imsave(fname, data_1[:,:,i], format='jpg')
        py.close('all')
    for i in range(len(data_2[0,0,:])):
        py.figure()
        sufix = str(i)
        path = 'data2_images/'
        name = 'data2'
        fname = path + name + sufix
        py.imshow(data_2[:,:,i])
        py.savefig(fname,format='png')
        py.close('all')        