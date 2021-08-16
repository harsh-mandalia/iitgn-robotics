from render import Renderer
import cv2
import numpy as np
from scipy.optimize import fsolve

class TwoR2(Renderer):
    def __init__(self):
        super().__init__()
        self.l1=100
        self.l2=100
        self.q1=np.pi/6
        self.q2=np.pi/2
        self.wall1=(100,10)
        self.wall2=(10,200)
        self.force=10

        self.reached=0
        self.x0=100
        self.y0=100
        self.tau1=0
        self.tau2=0
        
    def equations(self, p, q2,x1,x2,y1,y2):
        x0,y0,q1=p
        eq1=self.l1*np.cos(q1)+self.l2*np.cos(q2)-x0
        eq2=self.l1*np.sin(q1)+self.l2*np.sin(q2)-y0
         
        eq3=(y0-y1)*(x2-x1)-(x0-x1)*(y2-y1)
        return [eq1,eq2,eq3]

    def getInfo(self):
        info={"tau1":self.tau1,"tau2":self.tau2}
        return info

    def step(self,steps=100):
        x1,y1=self.wall1
        x2,y2=self.wall2
        if x2!=x1:
            wall_solpe=(y2-y1)/(x2-x1)
            q2=-np.arctan(1/wall_solpe)
        else:
            q2=0           
        fx=self.force*np.cos(q2)
        fy=self.force*np.sin(q2)

        self.x0,self.y0,q1=fsolve(self.equations, (x1,x2,self.q1),(q2,x1,x2,y1,y2))
        d1=q1-self.q1
        d2=q2-self.q2
        while(round(self.q1,2)!=round(q1,2) or round(self.q2,2)!=round(q2,2)):
            self.q1=self.q1+d1/steps
            self.q2=self.q2+d2/steps
            self.render()

        self.tau1 = fy*self.l1*np.cos(self.q1) - fx*self.l1*np.sin(self.q1)
        self.tau2 = fy*self.l2*np.cos(self.q2) - fx*self.l2*np.sin(self.q2)

    def draw(self, image):
        wall1=(self.wall1[0]+300,-self.wall1[1]+300)
        wall2=(self.wall2[0]+300,-self.wall2[1]+300)
        cv2.line(image, wall1, wall2, (0,0,0),2)

        line1=(int(300+self.l1*np.cos(-self.q1)),int(300+self.l1*np.sin(-self.q1)))
        cv2.line(image, (300,300), line1, (255,0,0), 1)
        line2=(int(line1[0]+self.l2*np.cos(-self.q2)),int(line1[1]+self.l2*np.sin(-self.q2)))
        cv2.line(image, line1, line2, (0,255,0), 1)

        center=(int(self.x0+300),int(-self.y0+300))
        cv2.circle(image, center, 1, (0,0,0),2)
        return image


mybot=TwoR2()
while(True):
    mybot.step(200)
    mybot.render()