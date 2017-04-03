import numpy as np
import math

L = (8, 8, 8, 8, 8, 8);
P = 
P = constructP(L, 15, 1)
print(constructH(L, P))


def geval(l, b, d):
    L = np.zeros(6)
    L.fill(l)
    P = constructP(L, b, d)

def constructP(L, b, d):
    P = np.zeros(3)
    for i in range(1,4):
        P[i] = (1 / (2 * b) * (b**2 + L[2 * i-1]**2 - L[i]**2))
    return P

def constructH(L,P):
    H = np.zeros(3)
    for i in range(1,4):
        H[i] = math.sqrt(L[2 * i - 1]**2 - P[i]**2)
    return H
        
    