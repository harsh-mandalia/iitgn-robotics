import cv2
import numpy as np
from numpy.lib.function_base import append
from scipy.integrate import solve_ivp, ode
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import sympy as sym 
from time import time
import pickle

def SCARA_invkin(p):
    x,y,z=p
    D=(x**2+y**2-l2**2-l3**2)/(2*l2*l3)
    q2=np.arctan2(np.sqrt(1-D**2),D)
    q1=np.arctan2(y,x)-np.arctan2(l3*np.sin(q2),l2+l3*np.cos(q2))
    
    d3=z-l1
    
    jp=np.array([np.rad2deg(q1),np.rad2deg(q2),d3])
    return q1,q2,d3

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

def joit_params(q1,q2,d3,vx,vy,vz,ax,ay,az):
    Jv=np.array([[-l3*np.sin(q1+q2)-l2*np.sin(q1), -l3*np.sin(q1+q2), 0],
                 [l3*np.cos(q1+q2)+l2*np.cos(q1), l3*np.cos(q1+q2), 0],
                 [0, 0, 1]])
    v=np.array([[vx],
                [vy],
                [vz]])
    qs=np.linalg.inv(Jv)@v
    q_dot=qs.T[0]
    
    q1_dot,q2_dot,d3_dot = q_dot
    a=np.array([[ax + l3*np.cos(q1+q2)*(q1_dot+q2_dot)**2 + l2*np.cos(q1)*q1_dot**2],
                [ay + l3*np.sin(q1+q2)*(q1_dot+q2_dot)**2 + l2*np.sin(q1)*q1_dot**2],
                [az]])
    qds=np.linalg.inv(Jv)@a
    q_dot_dot=qds.T[0]
    q1_dot_dot,q2_dot_dot,d3_dot_dot = q_dot_dot
    return q1_dot,q2_dot,d3_dot, q1_dot_dot,q2_dot_dot,d3_dot_dot



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

l1=0.25
l2=0.25
l3=0.25
m1=0.8
m2=0.8
m3=0.8
I1=0.005
I2=0.005
I3=0.8

t=0
q,q_,q__ = [],[],[]
time=[]
fig = plt.figure()
while(True):
    pos,pos_,pos__ = state(a,t)
    q1,q2,d3=SCARA_invkin(pos)
    q1_dot,q2_dot,d3_dot, q1_dot_dot,q2_dot_dot,d3_dot_dot=joit_params(q1,q2,d3, pos_[0],pos_[1],pos_[2], pos__[0],pos__[1],pos__[2])
    # q.append(q1)
    # q_.append(q2)
    # q__.append(d3)
    # q.append(q1_dot)
    # q_.append(q2_dot)
    # q__.append(d3_dot)
    q.append(q1_dot_dot)
    q_.append(q2_dot_dot)
    q__.append(d3_dot_dot)
    time.append(t)
    plt.plot(time,q)
    plt.plot(time,q_)
    plt.plot(time,q__)
    plt.pause(0.001)
    t+=0.1
    if(t>10):
        break
plt.legend(['q1_dot_dot','q2_dot_dot','d3_dot_dot'])
plt.ylabel('joint acceleration') 
plt.xlabel('time')
plt.title('joint acceleration vs time')

plt.show()