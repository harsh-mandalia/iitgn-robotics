# ANSWER 1
import numpy as np
import math
p = [] #end effector vector
for i in range(0,3):
    ele = float(input("Enter the End-Effector Coordinates:"))
    p.append(ele)
a_2 = float(input("Enter Link Length a_2:")) # Link Length 2
d_1 = float(input("Enter Distance d_1:")) #Distance between O_1 and O_0
r = math.sqrt(p[0]**2 + p[1]**2)
s = p[2]-d_1
theta_1 = math.atan2(p[1],p[0])
theta_2 = math.atan2(s,r)
d_3 = math.sqrt(r**2+s**2)-a_2
joint_var = [theta_1*180/math.pi, theta_2*180/math.pi, d_3] 
print("Joint Variables:",joint_var)

# ANSWER 2

# Rotation matrix relating end-effector frame to base frame
rot_matrix = [] 
print("Enter the entries row-wise:")
  
# For user input
for i in range(3):          # A for loop for row entries
    a =[]
    for j in range(3):      # A for loop for column entries
         a.append(int(input()))
    rot_matrix.append(a)
d = [] #end effector position vector
for i in range(3):
    s = "d" + "_" + str(i+1) + ":"
    ele = float(input(s)) 
    d.append(ele)
a_1 = float(input("Link Length 1:"))
a_2 = float(input("Link Length 2:"))
d_4 = float(input("Distance between O_3 and O_4:"))
c_2 = 0.5*(d[0]**2+d[1]**2-a_1**2-a_2**2)/(a_1*a_2)
theta_2 = math.atan2(c_2,math.sqrt(1-c_2**2))
theta_1 = math.atan2(d[1],d[0])-math.atan2(math.sin(theta_2),a_1+a_2*math.cos(theta_2))
theta_4 = theta_1 + theta_2 - math.atan2(rot_matrix[0][1],rot_matrix[0][0])
d_3 =  d[2] - d_4
joint_var = [theta_1, theta_2, theta_4, d_3]
print("Joint Variables:", joint_var)

# ANSWER 3

# The following segment of the code is for the Jacobian

D_H = []
A = []
o = []
Z = []
J = []
R_P = []
n = int(input("Number of links:"))
D_H = [0]*n
R_P = []
print("Link |   theta_i |   d_i |   a_i |   alpha_")
for i in range(n):
    a = []
    for j in range(4):
        a.append(int(input("Enter the parameters row-wise): ")))
    D_H.append(a)
for i in range(0,n):
    ele = input("Enter R (Revolute) or P (Prismatic) for joint " + str(i+1))
    R_P.append(ele)


for i in range(n):
    theta,d,a,alpha = D_H[i]
    A[i] = np.array([[np.cos(theta), -np.sin(theta)*np.cos(alpha), np.sin(theta)*np.sin(alpha), a*np.cos(theta)],
                [np.sin(theta), np.cos(theta)*np.cos(alpha), -np.cos(theta)*np.sin(alpha), a*np.sin(theta)],
                [0, np.sin(alpha), np.cos(alpha), d],
                [0, 0, 0, 1]])
# initialization
T = np.identity(4)
Z = [0]*(n +1)
Z[0] = np.array([[0],[0],[1]])
o = [0]*(n+1)
o[0] = np.array([0,0,0])

for i in range(n):
    T = T@A[i]
    Z[i+1]=T[0:3, 0:3]@Z[0]
    o[i+1]=T[0:3,3].T

End_pos = T[0:3,3] # end effector position 

for i in range(n):
    if R_P[i] == "P":
        J[i] = np.concatenate((Z[i], np.array([[0,0,0]]).T), axis=0)
    else:
        cross=np.cross(Z[i].T,(o[n]-o[i]).T)
        J[i]=np.concatenate((cross.T,Z[i]),axis=0)

final_J = J[0]
for i in range(1, n):
    final_J = np.concatenate((final_J, J[i]), axis=1)
J = final_J

# The above Jacobian was inspired from Harsh Mandalia's code 

# Now the following code is for Answer 3 of the Assignment 4

Xdot = []
qdot = []
for i in range(n):
    Xdot.append(int(input("Enter the velocties in the form (xdot, ydot, zdot): ")))
inv_J = np.linalg.pinv(J)
qdot = inv_J@Xdot
print(qdot)


# ANSWER 7

# Using the Jacobian again, we solve Question 7 and we can get the end-effector position from the Jacobian itself

End_vel = J@qdot # end effector velocity 

print(End_pos,End_vel)

# ANSWER 6
def Spherical_Wrist(U):
    if U[0][2]!=0 or U[1][2]!=0 :
        i = int(input("Enter 1 for sin_theta > 0 or 0 for sin_theta <0 :"))
        if i == 1:
            theta = math.atan2(math.sqrt(1-U[2][2]**2),U[2][2])
            phi = math.atan2(U[1][2],U[0][2])
            psi = math.atan2(U[2][1],-U[2][0])
        elif i == 0:
            theta = math.atan2(-math.sqrt(1-U[2][2]**2),U[2][2])
            phi = math.atan2(-U[1][2],-U[0][2])
            psi = math.atan2(-U[2][1],U[2][0])
        else:
            print("error")
    else:
        if U[2][2] == 1:
            theta = 0 # from U[2][2] = 1 
            phi_plus_psi = math.atan2(U[1][0],U[0][0]) # infinite solutions exist
            phi = 0; # in order to get a unique solution
            psi = phi_plus_psi - phi
        elif U[2][2] == -1:
            theta = math.pi # from U[2][2] = -1
            phi_minus_psi = math.atan2(-U[0][1],-U[0][0]) # infinite solutions exist
            phi = 0 # in order to get a unique solution
            psi = phi - phi_minus_psi
        else:
            print("error")
    euler = (theta*180/math.pi,phi*180/math.pi,psi*180/math.pi) # in degrees
    print("Euler Angles:",euler)