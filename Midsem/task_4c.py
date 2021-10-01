from render import Renderer
import cv2
import numpy as np
from scipy.integrate import ode

class TwoR3(Renderer):
    def __init__(self):
        super().__init__()
        self.ode = ode(self.func).set_integrator('vode', nsteps=500, method='bdf')
        self.l1=100
        self.q1=0.1
        self.q1_1=0
        self.m1=1
        self.g=9.81
        self.q0=0

        self.x0 = self.l1*np.cos(self.q0)
        self.y0 = self.l1*np.sin(self.q0)

        self.count=0
        self.q_sum=0
        self.q_av=0
        self.e_tot=0

        self.k=10
        self.tau1s=0
        self.first=True


    def func(self,t,y):
        q1=y[0]
        q1_1=y[1]
        
        m1, l1, g = self.m1, self.l1, self.g


        # temp=self.tau1s+m1*g*l1*(np.sin(q1)-0.5*np.cos(q1))
        # q1_2=3*temp/(3*l1**1)
        q1_2=self.tau1s/(m1*l1**2)
        
        return [q1_1,q1_2]

    def step(self,dt=0.01):
        x = self.l1*np.cos(self.q1)
        y = self.l1*np.sin(self.q1)
        
        self.q_sum+=self.q1
        self.count+=1
        self.q_av=self.q_sum/self.count
        self.e_tot=self.m1*self.g*self.l1*np.cos(self.q1)+1/6*self.m1*self.l1**2*self.q1_1**2
        # self.e_tot=0.5*self.m1*self.g*self.l1*np.sin(self.q1)+(1/6)*self.m1*self.l1**2*self.q1_1**2

        self.tau1s = -self.k*(self.q1-self.q0)
        
        if(self.first):
            state=[self.q1, self.q1_1]
            self.ode.set_initial_value(state,0)
            self.first=False
        newstate=self.ode.integrate(self.ode.t+dt)
        self.q1 = newstate[0]
        self.q1_1 = newstate[1]

    def draw(self, image):
        line1=(int(300+self.l1*np.cos(-self.q1)),int(300+self.l1*np.sin(-self.q1)))
        cv2.line(image, (300,300), line1, (255,0,0), 1)

        point=(int(300+self.x0),int(-self.y0+300))
        cv2.circle(image, point,2,(0,0,0),-1)

        return image

    def getInfo(self):
        info={"time":round(self.ode.t,3),"q average":round(self.q_av,3),"Total energy":round(self.e_tot,3)}
        return info

mybot=TwoR3()
while(True):
    mybot.step(0.1)
    mybot.render()