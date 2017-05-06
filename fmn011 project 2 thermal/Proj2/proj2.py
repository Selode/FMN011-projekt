#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 11:42:41 2017

@author: selode & boodle
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
    
    temp = np.array([50, 50, 50, 50, 50, 50, 50])
        
    root = fsolve(p, 1, args=2)
    
    k = 0.01
    
    plt.plot(z, T, 'o', label = "data points")
    
    plt.plot(xs, p(xs), label="S")          #Uppgift 1
    
    plt.axvline(x=root, linestyle='dashed', color='b', label="Thermocline depth")
    
    plt.plot(xs, -p(xs,1), label="Heat flux")
    
    plt.plot(xs, p(xs,1), label="S'")       #measures the derivate
    plt.plot(xs, p(xs,2), label="S'\'")     #measures the second derivate

    plt.axhline(50, linestyle='dotted', color='black', label="Temperature at 50 Celsius")
    plt.axvline(1.7, linestyle='dotted', color='black', label="Depth at 1.7 m")    
    
    plt.plot(root, 0, 'd', color='r')       #Thermocline depth point
    plt.plot(p.solve(50), 50, 'o', label="Function at T = 50 Celsius")   
    plt.plot(1.7, p(1.7), 'o', label = "Function at z = 1.7 m")
    
    
    plt.legend(loc='upper right', ncol=2)
    plt.show()
    
    print(p.solve(50))
    print(p(1.7))
    
    
tasks()
