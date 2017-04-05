import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import math





def task1(L, b, d):
    print("task 1:")
    print("minimum leg length: ")
    print(b/2.0)
    print("Height:")
    geval(L, b, d)

def task2(L, b, d, a):
    XT = (2.4537, -4.9075, 2.7424)
    print(XT)
    print stf(XT, L, b, d, a)
    

def task4(L, b, d, a):
    XT = (2.4537, -4.9075, 2.7424)
    print fsolve(stf, XT, (L,b,d,a))
    
def task5(L, b, d, a):
    startXT = (2.4537, -4.9075, 2.7424)
    #Legs at min length:
    L = (8, 8, 8, 8, 8, 8)
    XT = (fsolve(stf, startXT, (L,b,d,a)))
    P = constructP(L, b, d)
    H = constructH(L, P)
    XP = XPequations(P, b, d)
    YP = YPequations(P, b, d)
    YT = calcYT(XP, XT, YP)
    ZT = calcZT(H, XP, XT)
    #Legs at max length:
    L = (15, 15, 15, 15, 15, 15)
    calculatedXT = (fsolve(stf, startXT, (L,b,d,a)))
    
    #Maximally tilted platform:
    L = (15, 15, 8, 8, 8, 8)
    calculatedXT = (fsolve(stf, startXT, (L,b,d,a)))
    
    #Maximally twisted platform:
    L = (8, 15, 8, 15, 8, 15)
    calculatedXT = (fsolve(stf, startXT, (L,b,d,a)))
    
     
    
def geval(l, b, d):
    L = l
    P = constructP(L, b, d)
    H = constructH(L, P)

    print(H)
    print("XP: ")
    print(XPequations(P, b, d))
    print("YP: ")
    print(YPequations(P, b, d))

def constructP(L, b, d):
    P = np.zeros(3)
    for i in range(3):
        P[i] = (1.0 / (2 * b) * (b**2 + L[2 * (i+1) - 1]**2 - L[i]**2))
    return P

def constructH(L,P):
    H = np.zeros(3)
    for i in range(3):
        H[i] = math.sqrt(L[2 * (i+1) - 1]**2 - P[i]**2)
    return H

def calcYT(XP, XT, YP):
    YT = np.zeros(3)
    YT[0] = np.sqrt(3) * XT[0] - (np.sqrt(3) * XP[0] - YP[0])
    YT[1] = YP[1]
    YT[2] = -np.sqrt(3) * XT[2] + (np.sqrt(3) * XP[2] + YP[2])
    return YT
    
def calcZT(H, XP, XT):
    ZT = np.zeros(3)
    ZT[0] = np.sqrt(H[0]**2 - 4 * (XT[0] - XP[0])**2)
    ZT[1] = np.sqrt(H[1]**2 - (XT[1] - XP[1])**2)
    ZT[2] = np.sqrt(H[2]**2 - 4 * (XT[2] - XP[2])**2)
    return ZT

def XPequations(P, b, d):
    XP = np.zeros(3)
    XP[0] = math.sqrt(3) / 6.0 * (2 * b + d - 3 * P[0])
    XP[1] = - math.sqrt(3) / 6.0 * (b + 2 * d)
    XP[2] = - math.sqrt(3) / 6.0 * (b - d - 3 * P[2])
    return XP

def YPequations(P, b, d):
    YP = np.zeros(3)
    YP[0] = 0.5 * (d + P[0])
    YP[1] = 1 / 2.0 * (b - 2.0 * P[1])
    YP[2] = - 1 / 2.0 * (b + d - P[2])
    return YP

def stf(XT, L, b, d, a):
    finalXT = np.zeros(3)
    P = constructP(L, b, d)
    H = constructH(L, P)
    XP = XPequations(P, b, d)
    YP = YPequations(P, b, d)
    finalXT[0] = a**2 + 2 * XT[0] * XT[1] - 2 * XT[0] * (XP[0] + np.sqrt(3) * (YP[0] - YP[1])) - 2 * XP[1] * XT[1] - ((np.sqrt(3) * XP[0] - YP[0] + YP[1])**2 + (H[0]**2 + H[1]**2) - 4 * XP[0]**2 - XP[1]**2) + 2 * np.sqrt((H[0]**2 - 4 * (XT[0] - XP[0])**2) * (H[1]**2 - (XT[1] - XP[1])**2))
    finalXT[1] = a**2 - 4 * XT[0] * XT[2] - 2 * XT[0] *  (XP[0] - 3 * XP[2] + np.sqrt(3) * (YP[0] - YP[2])) - 2 * XT[2] * ((-3) * XP[0] + XP[2] + np.sqrt(3) * (YP[0] - YP[2])) - ((np.sqrt(3) * (XP[0] + XP[2]) - YP[0] + YP[2])**2 + (H[0]**2 + H[2]**2) - 4 * XP[0]**2 - 4 * XP[2]**2) + 2 * np.sqrt((H[0]**2 - 4 * (XT[0] - XP[0])**2) * (H[2]**2 - 4 * (XT[2] - XP[2])**2)) 
    finalXT[2] = a**2 + 2 * XT[1] * XT[2] - 2 * XT[2] * (XP[2] + np.sqrt(3) * (YP[1] - YP[2])) - 2 * XP[1] * XT[1] - ((np.sqrt(3) * XP[2] - YP[1] + YP[2])**2 + (H[1]**2 + H[2]**2) - XP[1]**2 - 4 * XP[2]**2) + 2 * np.sqrt((H[1]**2 - (XT[1] - XP[1])**2) * (H[2]**2 - 4 * (XT[2] - XP[2])**2))
    return finalXT


L = (8, 8, 8, 8, 8, 8)
b = 15
d = 1

task1(L, b, d)
L = (11.5, 11.5, 11.5, 11.5, 11.5, 11.5)
a = 10
task2(L, b, d, a)
task4(L, b, d, a)



