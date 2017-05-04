import scipy as sp
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as cs



def tasks():
    T = np.array([70, 70, 55, 22, 13, 10, 10])
    D = np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    p = cs.interp1d(D, T, kind = 'cubic')
    x = np.arange(0, 7)
    y = np.exp(-x/3.0)
    xnew = np.arange(0, 10, 0.1)
    ynew = cs.spline(D, T, (0, 3), 3, 'smoothest')
    plt.plot(np.array([0, 3]), np.array([0, 100]), 'o', xnew, ynew, '-')
    plt.show()
    
    
    
    
tasks()