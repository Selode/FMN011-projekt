import math
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection






def task1(L, b, d):
    print("task 1:")
    print("minimum leg length: ")
    print(b / 2.0)
    print("Height:")
    geval(L, b, d)


def task2(L, b, d, a):
    XT = (2.4537, -4.9075, 2.7424)
    print(XT)
    P = constructP(L, b, d)
    H = constructH(L, P)
    XP = XPequations(P, b, d)
    YP = YPequations(P, b, d)


    print (stf(XT, a, b, d, H, L, P, XP, YP))
    return



def task4(L, b, d, a):
    P = constructP(L, b, d)
    H = constructH(L, P)
    XP = XPequations(P, b, d)
    YP = YPequations(P, b, d)
    XT = (2.4537, -4.9075, 2.7424)

    print (fsolve(stf, XT, (a, b, d, H, L, P, XP, YP)))


def task5(L, b, d, a):
    startXT = (2.4537, -4.9075, 2.7424)
    # Legs at min length:
    #L = (8, 8, 8, 8, 8, 8)
    #Legs at max length:
    #L = (15, 15, 15, 15, 15, 15)
    # Maximally tilted platform:
    L = (15, 15, 8, 8, 8, 8)
    # Maximally twisted platform:
    #L = (8, 15, 8, 15, 8, 15)
    P = constructP(L, b, d)
    H = constructH(L, P)
    XP = XPequations(P, b, d)
    YP = YPequations(P, b, d)
    XT = (fsolve(stf, startXT, (a, b, d, H, L, P, XP, YP)))
    YT = calcYT(XP, XT, YP)
    ZT = calcZT(H, XP, XT)
    XB = calcBaseX(b, d)
    YB = calcBaseY(b, d)
    ZB = np.zeros(6)
    (L1, L2, L3) = calcLegs(XT, YT, ZT, XB, YB, ZB)
    plotGraph(XT, YT, ZT, XB, YB, L1, L2, L3)


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
        P[i] = (1.0 / (2 * b) * (b**2 + L[2 * i]**2 - L[2 * (i + 1) - 1]**2))
    return P


def constructH(L, P):
    H = np.zeros(3)
    for i in range(3):
        H[i] = math.sqrt(L[2 * (i + 1) - 1]**2 - P[i]**2)
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

def calcBaseX(b,d):
    BX = np.zeros(6)
    BX[0] = np.sqrt(3) * (2 * b + d)/6.0
    BX[1] = -np.sqrt(3) * (b - d)/6.0
    BX[2] = -np.sqrt(3) * (b + 2 * d)/6.0
    BX[3] = -np.sqrt(3) * (b + 2 * d)/6.0
    BX[4] = -np.sqrt(3) * (b - d)/6.0
    BX[5] = np.sqrt(3) * (2 * b + d)/6.0
    return BX

def calcBaseY(b, d):
    BY= np.zeros(6)
    BY[0] = d / 2.0
    BY[1] = (b + d) / 2.0
    BY[2] = b / 2.0
    BY[3] = -b / 2.0
    BY[4] = -(b+d) / 2.0
    BY[5] = -d / 2.0
    return BY

def calcLegs(xt,yt,zt,bx,by,bz):
    lx = np.zeros((6,2))
    ly = np.zeros((6,2))
    lz = np.zeros((6,2))
    for i in range(2):
        lx[i] = np.array([bx[i],xt[0]])
        ly[i] = np.array([by[i],yt[0]])
        lz[i] = np.array([bz[i],zt[0]])
    for i in range(2,4):
        lx[i] = np.array([bx[i],xt[1]])
        ly[i] = np.array([by[i],yt[1]])
        lz[i] = np.array([bz[i],zt[1]])
    for i in range(4,6):
        lx[i] = np.array([bx[i],xt[2]])
        ly[i] = np.array([by[i],yt[2]])
        lz[i] = np.array([bz[i],zt[2]])
    return lx,ly,lz

def stf(XT, a, b, d, H, L, P, XP, YP):
    finalXT = np.zeros(3)
    finalXT[0] = a**2 + 2 * XT[0] * XT[1] - 2 * XT[0] * (XP[0] + np.sqrt(3) * (YP[0] - YP[1])) - 2 * XP[1] * XT[1] - ((np.sqrt(3) * XP[0] - YP[0] + YP[1])**2 + (H[0]**2 + H[1]**2) - 4 * XP[0]**2 - XP[1]**2) + 2 * np.sqrt((H[0]**2 - 4 * (XT[0] - XP[0])**2) * (H[1]**2 - (XT[1] - XP[1])**2))
    finalXT[1] = a**2 - 4 * XT[0] * XT[2] - 2 * XT[0] * (XP[0] - 3 * XP[2] + np.sqrt(3) * (YP[0] - YP[2])) - 2 * XT[2] * ((-3) * XP[0] + XP[2] + np.sqrt(3) * (YP[0] - YP[2])) - ((np.sqrt(3) * (XP[0] + XP[2]) - YP[0] + YP[2])**2 + (H[0]**2 + H[2]**2) - 4 * XP[0]**2 - 4 * XP[2]**2) + 2 * np.sqrt((H[0]**2 - 4 * (XT[0] - XP[0])**2) * (H[2]**2 - 4 * (XT[2] - XP[2])**2))
    finalXT[2] = a**2 + 2 * XT[1] * XT[2] - 2 * XT[2] * (XP[2] + np.sqrt(3) * (YP[1] - YP[2])) - 2 * XP[1] * XT[1] - ((np.sqrt(3) * XP[2] - YP[1] + YP[2])**2 + (H[1]**2 + H[2]**2) - XP[1]**2 - 4 * XP[2]**2) + 2 * np.sqrt((H[1]**2 - (XT[1] - XP[1])**2) * (H[2]**2 - 4 * (XT[2] - XP[2])**2))
    
    return finalXT

def plotGraph(XT, YT, ZT, XB, YB, LX, LY, LZ):
    fig = plt.figure()
  
    ZB = np.zeros(6)
    
    ax = Axes3D(fig)
    top = ax.plot_trisurf(XT, YT, ZT, color= 'blue')
    bottom = ax.plot_trisurf(XB, YB, ZB, color= 'red')
    ax.add_collection3d(top)
    ax.add_collection3d(bottom)
    for i in range(6):
        ax.add_collection3d(ax.plot_wireframe(LX[i], LY[i], LZ[i], linewidths = 5, color = 'black'))

    plt.show()
    

L = (8, 8, 8, 8, 8, 8)
b = 15
d = 1

task1(L, b, d)
L = (11.5, 11.5, 11.5, 11.5, 11.5, 11.5)
a = 10
task2(L, b, d, a)
task4(L, b, d, a)
task5(L, b, d, a)
