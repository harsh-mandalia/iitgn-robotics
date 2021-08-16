from render import Renderer
import cv2
import numpy as np
from scipy.integrate import ode

class TwoR3(Renderer):
    def __init__(self):
        super().__init__()
        self.ode = ode(self.func).set_integrator('vode', nsteps=500, method='bdf')
        self.l1=100
        self.l2=100
        self.q1=-(np.pi/2+np.pi/3)
        self.q2=-(np.pi/6)
        self.q1_1=0
        self.q2_1=0
        self.m1=1
        self.m2=1
        self.g=9.81
        self.x0=0
        self.y0=-105


        self.x_sum=0
        self.y_sum=0
        self.count=0
        self.x_av=0
        self.y_av=0

        self.k=10
        self.tau1s=0
        self.tau2s=0
        self.first=True


    def func(self,t,y):
        q1=y[0]
        q1_1=y[1]
        q2=y[2]
        q2_1=y[3]
        
        m1, m2, l1, l2, g = self.m1, self.m2, self.l1, self.l2, self.g
        a = -0.5*m2*l1*l2*q2_1*(q2_1-q1_1)*np.sin(q2-q1) - 0.5*m2*l1*l2*q1_1*q2_1*np.sin(q2-q1) + 0.5*m1*g*l1*np.cos(q1) + m2*g*l1*np.cos(q1) - self.tau1s
        b = -0.5*m2*l1*l2*q1_1*(q2_1-q1_1)*np.sin(q2-q1) + 0.5*m2*l1*l2*q1_1*q2_1*np.sin(q2-q1) + 0.5*m2*g*l2*np.cos(q2) - self.tau2s
        c = 1/3*m1*l1**2 + m2*l1**2
        d = 0.5*m2*l1*l2*np.cos(q2-q1)
        e = 1/12*m2*l1**2 + 1/4*m2*l2**2
        f = 0.5*m2*l1*l2*np.cos(q2-q1)

        q2_2 = (a/c-b/f)/(e/f-d/c)
        q1_2 = -a/c - d*q2_2/c
        return [q1_1,q1_2,q2_1,q2_2]

    def step(self,dt=0.01):
        x = self.l1*np.cos(self.q1)+self.l2*np.cos(self.q2)
        y = self.l1*np.sin(self.q1)+self.l2*np.sin(self.q2)
        
        self.x_sum+=x
        self.y_sum+=y
        self.count+=1
        self.x_av=self.x_sum/self.count
        self.y_av=self.y_sum/self.count

        fx = -self.k*(x-self.x0)
        fy = -self.k*(y-self.y0)
        self.tau1s = fy*self.l1*np.cos(self.q1) - fx*self.l1*np.sin(self.q1)   + 0.5*self.m2*self.g*self.l2*np.cos(self.q2)
        self.tau2s = fy*self.l2*np.cos(self.q2) - fx*self.l2*np.sin(self.q2)   + 0.5*self.m1*self.g*self.l1*np.cos(self.q1)
        
        if(self.first):
            state=[self.q1, self.q1_1, self.q2, self.q2_1]
            self.ode.set_initial_value(state,0)
            self.first=False
        newstate=self.ode.integrate(self.ode.t+dt)
        self.q1 = newstate[0]
        self.q2 = newstate[2]

    def draw(self, image):
        line1=(int(300+self.l1*np.cos(-self.q1)),int(300+self.l1*np.sin(-self.q1)))
        cv2.line(image, (300,300), line1, (255,0,0), 1)
        line2=(int(line1[0]+self.l2*np.cos(-self.q2)),int(line1[1]+self.l2*np.sin(-self.q2)))
        cv2.line(image, line1, line2, (0,255,0), 1)

        cv2.circle(image, (300+self.x0,-self.y0+300),2,(0,0,0),-1)

        return image

    def getInfo(self):
        info={"time":round(self.ode.t,3),"x average":round(self.x_av,3),"y average":round(self.y_av,3)}
        return info

mybot=TwoR3()
while(True):
    mybot.step(0.01)
    mybot.render()