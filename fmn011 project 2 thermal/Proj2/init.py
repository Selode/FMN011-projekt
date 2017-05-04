#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 11:42:41 2017

@author: selode
"""

import math
import scipy as sp
import numpy as np
import scipy.interpolate as cs
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
#from scipy.interpolate import CubicSpline

def tasks():
    T = np.array([ 70, 70, 55, 22, 13, 10, 10])
    z = np.array([ 0, 0.5, 1, 1.5, 2, 2.5, 3])
   #z2 = np.array([ 0, 50, 100, 150, 200, 250, 300])
    
    p = cs.CubicSpline(z, T, bc_type='clamped')
    
    xs = np.arange(0,3,0.01)
        
    root = fsolve(p, 1, args=2)
    
    k = 0.01
    
    plt.plot(z, T, 'o', label = "data points")
    plt.plot(xs, p(xs), label="S")
    
    plt.axvline(x=root, linestyle='dashed', color='b', label="Thermal depth")
   
    
    plt.plot(xs, -p(xs,1), label="Heat flux")
    
    plt.plot(xs, p(xs,1), label="S'")       #measures the derivate
    plt.plot(xs, p(xs,2), label="S'\'")     #measures the second derivate

    plt.axhline(0, linestyle='dotted', color='black')
    plt.axvline(0, linestyle='dotted', color='black')    
    plt.legend(loc='upper right', ncol=2)
    
    plt.plot(root, 0, 'd', color='r')
    plt.show()
    
    

tasks()
