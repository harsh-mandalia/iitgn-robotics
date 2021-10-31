from operator import eq
import numpy as np
import matplotlib.pyplot as plt
from time import time
import pickle

def SCARA_fkin(q1,q2,d3,k):
    n=3
    dh=np.array([[q1, l1, l2, 0],
                 [q2, 0, l2, 0],
                 [0, d3, 0, 0]])

    A=[0]*n
    for i in range(n):
        theta, d, a, alpha = dh[i]
        A[i]=np.array([[np.cos(theta), -np.sin(theta)*np.cos(alpha), np.sin(theta)*np.sin(alpha), a*np.cos(theta)],
                       [np.sin(theta), np.cos(theta)*np.cos(alpha), -np.cos(theta)*np.sin(alpha), a*np.sin(theta)],
                       [0, np.sin(alpha), np.cos(alpha), d],
                       [0, 0, 0, 1]])
    ans=np.identity(4)

    for i in range(k):
        ans=ans@A[i]
    pos=ans[0:3,3]
    return pos

def PUMA_fkin(q1,q2,q3,k):
    n=3
    dh=np.array([[q1, l1, 0, np.pi/2],
                 [q2, 0, l2, 0],
                 [q3, 0, l3, 0]])

    A=[0]*n
    for i in range(n):
        theta, d, a, alpha = dh[i]
        A[i]=np.array([[np.cos(theta), -np.sin(theta)*np.cos(alpha), np.sin(theta)*np.sin(alpha), a*np.cos(theta)],
                       [np.sin(theta), np.cos(theta)*np.cos(alpha), -np.cos(theta)*np.sin(alpha), a*np.sin(theta)],
                       [0, np.sin(alpha), np.cos(alpha), d],
                       [0, 0, 0, 1]])
    ans=np.identity(4)

    for i in range(k):
        ans=ans@A[i]
    pos=ans[0:3,3]
    return pos

def Stanford_fkin(q1,q2,d3,k):
    n=3
    dh=np.array([[q1, l1, 0, np.pi/2],
                 [np.pi/2+q2, 0, 0, np.pi/2],
                 [0, l2+d3, 0, 0]])
    A=[0]*n
    for i in range(n):
        theta, d, a, alpha = dh[i]
        A[i]=np.array([[np.cos(theta), -np.sin(theta)*np.cos(alpha), np.sin(theta)*np.sin(alpha), a*np.cos(theta)],
                       [np.sin(theta), np.cos(theta)*np.cos(alpha), -np.cos(theta)*np.sin(alpha), a*np.sin(theta)],
                       [0, np.sin(alpha), np.cos(alpha), d],
                       [0, 0, 0, 1]])
    ans=np.identity(4)

    for i in range(k):
        ans=ans@A[i]
    pos=ans[0:3,3]
    return pos

def animate_angle(q1,q2,d3,dt,func):
    plt.clf()
    ax = plt.axes(projection='3d')
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_zlim(0, 2)
    O0=(0, 0, 0)
    O1=(0, 0, l1)
    # O1=SCARA_fkin(q1,q2,d3,1)
    # O2=(l2*np.cos(q1), l2*np.sin(q1), l1)
    O2=func(q1,q2,d3,1)
    # O3=(l2*np.cos(q1) + l3*np.cos(q1+q2), l2*np.sin(q1) + l3*np.sin(q1+q2), l1)
    O3=func(q1,q2,d3,2)
    # O4=(l2*np.cos(q1) + l3*np.cos(q1+q2), l2*np.sin(q1) + l3*np.sin(q1+q2), l1+d3)
    O4=func(q1,q2,d3,3)
    ax.plot3D([O0[0], O1[0]], [O0[1], O1[1]], [O0[2], O1[2]],'-o')
    ax.plot3D([O1[0], O2[0]], [O1[1], O2[1]], [O1[2], O2[2]],'-o')
    ax.plot3D([O2[0], O3[0]], [O2[1], O3[1]], [O2[2], O3[2]],'-o')
    ax.plot3D([O3[0], O4[0]], [O3[1], O4[1]], [O3[2], O4[2]],'-o')

    ax.plot3D(p1[0], p1[1], p1[2],'*')
    ax.plot3D(p2[0], p2[1], p2[2],'*')

    return plt.waitforbuttonpress(dt)

bot_name="stanford"
func=SCARA_fkin
path=""
if(bot_name=="stanford"):
    p1=(1,1.5,0)
    p2=(1.5,1,1)
    func=Stanford_fkin
    path=path="data\stanford.pkl"
elif(bot_name=="scara"):
    p1=(1,1.5,0)
    p2=(1.5,1,1)
    func=SCARA_fkin
    path=path="data\scara.pkl"
elif(bot_name=="puma"):
    p1=(0.5,0.5,0)
    p2=(1,0.7,0.5)
    func=PUMA_fkin
    path=path="data\puma.pkl"


l1,l2,l3=1,1,1

plt.ion()
plt.show()
fig = plt.figure()

myvar=0
with open(path, 'rb') as file:
    myvar = pickle.load(file)
q1s,q2s,d3s=myvar

print(len(q1s), "iterations")
for i in range(len(q1s)):
    if animate_angle(q1s[i],q2s[i],d3s[i],0.0001,func):
        break

plt.plot(q1s)
plt.plot(q2s)
plt.plot(d3s)

plt.ioff()
plt.show()