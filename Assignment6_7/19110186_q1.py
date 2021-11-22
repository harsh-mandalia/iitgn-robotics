import cv2
import numpy as np
from numpy.lib.function_base import append
from scipy.integrate import solve_ivp, ode
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import sympy as sym 
from time import time
import pickle

y0=0.06
y_0=0
y__0=0
y1=0.01
y_1=0
y__1=0

t0=0
t1=10

b=np.array([[y0],
            [y_0],
            [y__0],
            [y1],
            [y_1],
            [y__1]])
M=np.array([[1, t0, t0**2, t0**3, t0**4, t0**5],
            [0, 1, 2*t0, 3*t0**2, 4*t0**3, 5*t0**4],
            [0, 0, 2, 6*t0, 12*t0**2, 20*t0**3],
            [1, t1, t1**2, t1**3, t1**4, t1**5],
            [0, 1, 2*t1, 3*t1**2, 4*t1**3, 5*t1**4],
            [0, 0, 2, 6*t1, 12*t1**2, 20*t1**3]])

a=np.linalg.inv(M)@b

def state(a,t):
    a=a.T[0]
    x = 0.4
    y = a[0] + a[1]*t + a[2]*t**2 + a[3]*t**3 + a[4]*t**4 + a[5]*t**5
    z = 0.1

    x_ = 0
    y_ = a[1] + 2*a[2]*t + 3*a[3]*t**2 + 4*a[4]*t**3 + 5*a[5]*t**4
    z_ = 0

    x__ = 0
    y__ = 2*a[2] + 6*a[3]*t + 12*a[4]*t**2 + 20*a[5]*t**3
    z__ = 0

    pos=[x,y,z]
    pos_=[x_,y_,z_]
    pos__=[x__,y__,z__]
    return pos,pos_,pos__

t=0
x,y,z = [],[],[]
time=[]
fig = plt.figure()
while(True):
    temp,temp2,temp3 = state(a,t)
    x.append(temp3[0])
    y.append(temp3[1])
    z.append(temp3[2])
    time.append(t)
    plt.plot(time,x)
    plt.plot(time,y)
    plt.plot(time,z)
    plt.pause(0.001)
    t+=0.1
    if(t>10):
        break
plt.legend(['x_dot_dot','y_dot_dot','z_dot_dot'])
plt.ylabel('acceleration') 
plt.xlabel('time')
plt.title('acceleration vs time')

plt.show()