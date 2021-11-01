from operator import eq
import cv2
import numpy as np
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

def SCARA_fkin(q1,q2,d3,k):
    n=3
    dh=np.array([[q1, l1, l2, 0],
                 [q2, 0, l3, 0],
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

def D_calculator():
    q1 = sym.Symbol('q1')
    q1_dot = sym.Symbol('q1_dot')
    q1_dot_dot = sym.Symbol('q1_dot_dot')
    q2 = sym.Symbol('q2')
    q2_dot = sym.Symbol('q2_dot')
    q2_dot_dot = sym.Symbol('q2_dot_dot')
    q3 = sym.Symbol('q3')
    q3_dot = sym.Symbol('q3_dot')
    q3_dot_dot = sym.Symbol('q3_dot_dot')

    jv1=np.array([[-l2/2*sym.sin(q1), 0, 0],
                  [l2/2*sym.cos(q1), 0, 0],
                  [0, 0, 0]])
    jv2=np.array([[-l2*sym.sin(q1)-l3/2*sym.sin(q1+q2), -l3/2*sym.sin(q1+q2), 0],
                  [l2*sym.cos(q1)+l3/2*sym.cos(q1+q2), l3/2*sym.cos(q1+q2), 0],
                  [0, 0, 0]])
    jv3=np.array([[-l2*sym.sin(q1)-l3*sym.sin(q1+q2), -l3*sym.sin(q1+q2), 0],
                  [l2*sym.cos(q1)+l3*sym.cos(q1+q2), l3*sym.cos(q1+q2), 0],
                  [0, 0, 1/2]])
    D1=m1*jv1.T@jv1+m2*jv2.T@jv2+m3*jv3.T@jv3
    D2_1=np.array([[I1, 0, 0],
                   [0, 0, 0],
                   [0, 0, 0]])
    D2_2=np.array([[I2, I2, 0],
                   [I2, I2, 0],
                   [0, 0, 0]])
    D2=D2_1+D2_2
    D=D1+D2
    return D

def dynamical_equation(D):
    n=3
    q1 = sym.Symbol('q1')
    q1_dot = sym.Symbol('q1_dot')
    q1_dot_dot = sym.Symbol('q1_dot_dot')
    q2 = sym.Symbol('q2')
    q2_dot = sym.Symbol('q2_dot')
    q2_dot_dot = sym.Symbol('q2_dot_dot')
    q3 = sym.Symbol('q3')
    q3_dot = sym.Symbol('q3_dot')
    q3_dot_dot = sym.Symbol('q3_dot_dot')

    # D=np.array([[m1*l2**2/3+m2*l2**2, m2*l1*l2/2*sym.cos(q2-q1)],
    #             [m2*l1*l2/2*sym.cos(q2-q1), m2*l2**2/3]])
    # V=m1*g*l1/2*sym.sin(q1)+m2*g*(l1*sym.sin(q1)+l2/2*sym.sin(q2))
    V = g*(m1*l1+m2*l1+m3*(l1+q3/2))

    phi=np.array([[sym.diff(V, q1)],
                  [sym.diff(V, q2)],
                  [sym.diff(V, q3)]])
    q=np.array([[q1],
                [q2],
                [q3]])
    q_dot=np.array([[q1_dot],
                    [q2_dot],
                    [q3_dot]])
    q_dot_dot=np.array([[q1_dot_dot],
                        [q2_dot_dot],
                        [q3_dot_dot]])
    c=[0]*n
    for k in range(n):
        for i in range(n):
            for j in range(n):
                sum=sym.diff(D[k][j], "q"+str(i+1))+sym.diff(D[k][i], "q"+str(j+1))+sym.diff(D[i][j], "q"+str(k+1))
                c[k]+=0.5*(sum)*sym.Symbol("q"+str(i+1))*sym.Symbol("q"+str(j+1))
    eqn=sym.Array(D@q_dot_dot+phi+np.transpose([c]))
    return eqn

def func(t,y):
    q1=y[0]
    q1_dot=y[1]
    q2=y[2]
    q2_dot=y[3]
    q3=y[4]
    q3_dot=y[5]
    q_dot=[q1_dot, q2_dot, q3_dot]
    q=[q1,q2,q3]

    temp=eqn.subs([('q1_dot',q_dot[0]),('q2_dot',q_dot[1]),('q3_dot',q_dot[2]), ('q1',q[0]),('q2',q[1]),('q3',q[2])])
    # temp=fsolve(equation_solve, (0, 0, 0), ([q1,q2,q3], [q1_dot,q2_dot,q3_dot]),xtol=1)
    # print(temp)
    a1=temp[0][0].coeff('q1_dot_dot')
    b1=temp[0][0].coeff('q2_dot_dot')
    c1=temp[0][0].coeff('q3_dot_dot')
    d1=temp[0][0].coeff('1')
    a2=temp[1][0].coeff('q1_dot_dot')
    b2=temp[1][0].coeff('q2_dot_dot')
    c2=temp[1][0].coeff('q3_dot_dot')
    d2=temp[1][0].coeff('1')
    a3=temp[2][0].coeff('q1_dot_dot')
    b3=temp[2][0].coeff('q2_dot_dot')
    c3=temp[2][0].coeff('q3_dot_dot')
    d3=temp[2][0].coeff('1')
    
    M=np.array([[a1, b1, c1],
                [a2, b2, c2],
                [a3, b3, c3]],dtype="float")
    T=np.array([[t1-d1],
                [t2-d2],
                [f3-d3]])
    temp=np.linalg.inv(M)@T
    # print(M)
    # print(M@temp-T)
    # print(temp)
    # print(t1,t2,f3)
    # temp=[t1,t2,f3]
    q1_dot_dot=temp[0][0]
    q2_dot_dot=temp[1][0]
    q3_dot_dot=temp[2][0]

    return [q1_dot,q1_dot_dot, q2_dot,q2_dot_dot, q3_dot,q3_dot_dot]

def animate_angle(q1,q2,d3,dt):
    plt.clf()
    ax = plt.axes(projection='3d')
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_zlim(0, 2)
    O0=(0, 0, 0)
    O1=(0, 0, l1)
    # O1=SCARA_fkin(q1,q2,d3,1)
    # O2=(l2*np.cos(q1), l2*np.sin(q1), l1)
    O2=SCARA_fkin(q1,q2,d3,1)
    # O3=(l2*np.cos(q1) + l3*np.cos(q1+q2), l2*np.sin(q1) + l3*np.sin(q1+q2), l1)
    O3=SCARA_fkin(q1,q2,d3,2)
    # O4=(l2*np.cos(q1) + l3*np.cos(q1+q2), l2*np.sin(q1) + l3*np.sin(q1+q2), l1+d3)
    O4=SCARA_fkin(q1,q2,d3,3)
    ax.plot3D([O0[0], O1[0]], [O0[1], O1[1]], [O0[2], O1[2]],'-o')
    ax.plot3D([O1[0], O2[0]], [O1[1], O2[1]], [O1[2], O2[2]],'-o')
    ax.plot3D([O2[0], O3[0]], [O2[1], O3[1]], [O2[2], O3[2]],'-o')
    ax.plot3D([O3[0], O4[0]], [O3[1], O4[1]], [O3[2], O4[2]],'-o')

    ax.plot3D(p1[0], p1[1], p1[2],'*')
    ax.plot3D(p2[0], p2[1], p2[2],'*')

    return plt.waitforbuttonpress(dt)

def animate_torque(q1,q2,d3,t1,t2,f3,dt):
    # q1+=0.1*t1+t2/100
    # q2+=0.1*t2+f3/10
    # d3+=0.1*f3+(t1+t2)/100
    newstate=ode_eqn.integrate(ode_eqn.t+0.1)
    q1=newstate[0]
    q1_dot=newstate[1]
    q2=newstate[2]
    q2_dot=newstate[3]
    d3=newstate[4]
    d3_dot=newstate[5]
    # print(q1)
    return q1,q2,d3, q1_dot,q2_dot,d3_dot, animate_angle(q1,q2,d3,dt)
    # return q1,q2,d3,False

g=9.81
l1,l2,l3=1,1,1
m1,m2,m3=1,1,1
I1=1/12*m1*l2**2
I2=1/12*m2*l3**2

p1=(1,1.5,0)
p2=(1.5,1,1)
# p2=(1.1,1.6,0.1)

q1_1,q2_1,d3_1=SCARA_invkin(p1)
q1_2,q2_2,d3_2=SCARA_invkin(p2)
D=sym.simplify(D_calculator())
eqn1=dynamical_equation(D)
eqn=sym.simplify(eqn1)
print(eqn)
ode_eqn=ode(func).set_integrator('vode', nsteps=5, method='bdf')
state = [q1_1,0, q2_1,0, d3_1,0]
ode_eqn.set_initial_value(state,0)

plt.ion()
plt.show()
fig = plt.figure()

# q1=np.linspace(q1_1, q1_2, 100)
# q2=np.linspace(q2_1, q2_2, 100)
# d3=np.linspace(d3_1, d3_2, 100)

# for i in range(100):
#     if animate_angle(q1[i],q2[i],d3[i],0.001):
#         break

# animate_angle(q1_1,q2_1,d3_1,1)

ki=0.01
kp=100
kd=10
q1,q2,d3=q1_1,q2_1,d3_1
q1_dot,q2_dot,d3_dot=0,0,0
dt=0
q1_i,q2_i,d3_i=0,0,0
q1s=[]
q2s=[]
d3s=[]
# for i in range(1000):
while(True):
    time1=time()
    q1_i+=(q1_2-q1)*dt
    q2_i+=(q2_2-q2)*dt
    d3_i+=(d3_2-d3)*dt
    
    t1 = ki*q1_i + kp*(q1_2-q1) - kd*q1_dot
    t2 = ki*q2_i + kp*(q2_2-q2) - kd*q2_dot
    f3 = ki/100*d3_i + kp/10*(d3_2-d3) - kd/10*d3_dot
    
    flag=True
    q1,q2,d3, q1_dot,q2_dot,d3_dot, flag=animate_torque(q1,q2,d3,t1,t2,f3,0.0001)
    print(q1,t1, "->", ki*q1_i,kp*(q1_2-q1), -kd*q1_dot)
    # print(q1, t1)
    q1s.append(q1)
    q2s.append(q2)
    d3s.append(d3)
    time2=time()
    dt=time2-time1
    if(flag):
        break

with open('data\scara7.pkl', 'wb') as file:
    pickle.dump([q1s,q2s,d3s], file)

# myvar=0
# with open('data\scara7.pkl', 'rb') as file:
#     myvar = pickle.load(file)
# q1s,q2s,d3s=myvar
plt.plot(q1s)
plt.plot(q2s)
plt.plot(d3s)

# print(len(q1s))
# for i in range(len(q1s)):
#     if animate_angle(q1s[i],q2s[i],d3s[i],0.0001):
#         break

plt.ioff()
plt.show()