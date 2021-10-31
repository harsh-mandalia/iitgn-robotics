import cv2
import numpy as np
from scipy.integrate import solve_ivp, ode
import matplotlib.pyplot as plt
from time import time
import sympy as sym
import pickle

def PUMA_invkin(p):
    x,y,z=p
    q1=np.arctan2(y,x)
    r=np.sqrt(x**2+y**2)
    s=z-l1
    
    D=(r**2+s**2-l2**2-l3**2)/(2*l2*l3)
    q3=-np.arctan2(np.sqrt(1-D**2),D)  # negative because we want the end effector to have downwards orientation
    q2=np.arctan2(s,r)-np.arctan2(l3*np.sin(q3),l2+l3*np.cos(q3))
    
    jp=np.array([np.rad2deg(q1),np.rad2deg(q2),np.rad2deg(q3)])
    return q1,q2,q3

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

    jv1=np.array([[0, 0, 0],
                  [0, 0, 0],
                  [0, 0, 0]])
    jv2=np.array([[-l2/2*sym.sin(q1)*sym.cos(q2), -l2/2*sym.cos(q1)*sym.sin(q2), 0],
                  [l2/2*sym.cos(q1)*sym.cos(q2), -l2/2*sym.sin(q1)*sym.sin(q2), 0],
                  [0, l2/2*sym.cos(q2), 0]])
    jv3=np.array([[-l2*sym.sin(q1)*sym.cos(q2)-l3/2*sym.sin(q1)*sym.cos(q3 + q2), -l2*sym.cos(q1)*sym.sin(q2)-l3/2*sym.sin(q3 + q2)*sym.cos(q1), -l3/2*sym.sin(q3 + q2)*sym.cos(q1)],
                  [l2*sym.cos(q1)*sym.cos(q2)+l3/2*sym.cos(q1)*sym.cos(q3 + q2), -l2*sym.sin(q1)*sym.sin(q2)-l3/2*sym.sin(q3 + q2)*sym.sin(q1), -l3/2*sym.sin(q3 + q2)*sym.sin(q1)],
                  [0, l2*sym.cos(q2)+l3/2*sym.cos(q2+q3), l3/2*sym.cos(q2+q3)]])

    D1=m1*jv1.T@jv1+m2*jv2.T@jv2+m3*jv3.T@jv3
    return D1

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
    V=m1*g*l1/2 + m2*g*(l1+l2/2*sym.sin(q1)) + m3*g*(l1+l2*sym.sin(q2)+l3/2*sym.sin(q2+q3))

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
                sum=sym.diff(D[k][j], "q"+str(i))+sym.diff(D[k][i], "q"+str(j))+sym.diff(D[i][j], "q"+str(k))
                c[k]+=0.5*(sum)*sym.Symbol("q"+str(i))*sym.Symbol("q"+str(j))
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

def animate_angle(q1,q2,q3,dt):
    plt.clf()
    ax = plt.axes(projection='3d')
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_zlim(0, 2)
    O0=(0, 0, 0)
    O1=(0, 0, l1)
    O2=PUMA_fkin(q1,q2,q3,1)
    O3=PUMA_fkin(q1,q2,q3,2)
    O4=PUMA_fkin(q1,q2,q3,3)
    ax.plot3D([O0[0], O1[0]], [O0[1], O1[1]], [O0[2], O1[2]],'-o')
    ax.plot3D([O1[0], O2[0]], [O1[1], O2[1]], [O1[2], O2[2]],'-o')
    ax.plot3D([O2[0], O3[0]], [O2[1], O3[1]], [O2[2], O3[2]],'-o')
    ax.plot3D([O3[0], O4[0]], [O3[1], O4[1]], [O3[2], O4[2]],'-o')

    ax.plot3D(p1[0], p1[1], p1[2],'*')
    ax.plot3D(p2[0], p2[1], p2[2],'*')

    return plt.waitforbuttonpress(dt)

def animate_torque(q1,q2,q3,t1,t2,t3,dt):
    newstate=ode_eqn.integrate(ode_eqn.t+1)
    q1=newstate[0]
    q1_dot=newstate[1]
    q2=newstate[2]
    q2_dot=newstate[3]
    q3=newstate[4]
    q3_dot=newstate[5]
    # print(q1)
    return q1,q2,q3, q1_dot,q2_dot,q3_dot, animate_angle(q1,q2,q3,dt)
    # return q1,q2,q3,False

g=9.81
l1,l2,l3=1,1,1
m1,m2,m3=1,1,1

p1=(0.5,0.5,0)
p2=(1,0.7,0.5)
# p1=(1,1.5,0)
# p2=(1.5,1,1)

q1_1,q2_1,q3_1=PUMA_invkin(p1)
q1_2,q2_2,q3_2=PUMA_invkin(p2)

D=sym.simplify(D_calculator())
eqn1=dynamical_equation(D)
eqn=sym.simplify(eqn1)
print(eqn)
ode_eqn=ode(func).set_integrator('vode', nsteps=5, method='bdf')
state = [q1_1,0, q2_1,0, q3_1,0]
ode_eqn.set_initial_value(state,0)

plt.ion()
plt.show()
fig = plt.figure()

# q1=np.linspace(q1_1, q1_2, 100)
# q2=np.linspace(q2_1, q2_2, 100)
# q3=np.linspace(q3_1, q3_2, 100)

# for i in range(100):
#     if animate_angle(q1[i],q2[i],q3[i],0.001):
#         break

# animate_angle(0,0,np.pi/4,0.0001)

ki=0.01
kp=1000
kd=10
q1,q2,q3=q1_1,q2_1,q3_1
q1_dot,q2_dot,q3_dot=0,0,0
dt=0
q1_i,q2_i,q3_i=0,0,0
q1s=[]
q2s=[]
q3s=[]
# for i in range(1000):
while(True):
    time1=time()
    q1_i+=(q1_2-q1)*dt
    q2_i+=(q2_2-q2)*dt
    q3_i+=(q3_2-q3)*dt
    
    t1 = ki*q1_i + kp*(q1_2-q1) - kd*q1_dot
    t2 = ki*q2_i + kp*(q2_2-q2) - kd*q2_dot
    f3 = ki/100*q3_i + kp/10*(q3_2-q3) - kd/10*q3_dot
    
    print(ki*q1_i,kp*(q1_2-q1), -kd*q1_dot)

    flag=True
    q1,q2,q3, q1_dot,q2_dot,q3_dot, flag=animate_torque(q1,q2,q3,t1,t2,f3,0.0001)
    print(q1, t1)
    q1s.append(q1)
    q2s.append(q2)
    q3s.append(q3)
    time2=time()
    dt=time2-time1
    if(flag):
        break

with open('data\puma1.pkl', 'wb') as file:
    pickle.dump([q1s,q2s,q3s], file)

# myvar=0
# with open('data\scara5.pkl', 'rb') as file:
#     myvar = pickle.load(file)
# q1s,q2s,q3s=myvar
plt.plot(q1s)
plt.plot(q2s)
plt.plot(q3s)

# print(len(q1s))
# for i in range(len(q1s)):
#     if animate_angle(q1s[i],q2s[i],q3s[i],0.0001):
#         break

plt.ioff()
plt.show()