from render import Renderer
import cv2
import numpy as np
from scipy.optimize import fsolve

class TwoR1(Renderer):
    def __init__(self):
        super().__init__()
        self.l1=100
        self.l2=100
        self.q1=np.pi/6
        self.q2=np.pi/2
        self.x=300
        self.y=300
        self.path=[]

    
    def equations(self, p):
        q1,q2=p
        eq1=self.l1*np.cos(q1)+self.l2*np.cos(q2)-self.x
        eq2=self.l1*np.sin(q1)+self.l2*np.sin(q2)-self.y
        return [eq1,eq2]

    def callback(self, event,x,y,flags,param):
        if event == cv2.EVENT_MOUSEMOVE:
            self.x=x-300
            self.y=-y+300
            q1,q2=fsolve(self.equations, (self.q1, self.q2))
            min=35*np.pi/180
            max=145*np.pi/180
            if(self.x**2+self.y**2)<=(self.l1+self.l2)**2:
                self.q1=q1
                self.q2=q2
                point=(x,y)
                self.path.append(point)
            else:
                print("I can't reach ",end="")
            print("({}, {})".format(self.x,self.y))
            
    def step(self,steps=10):
        cv2.setMouseCallback("window",self.callback)       

    def getInfo(self):
        info={}
        return info

    def draw(self, image):
        line1=(int(300+self.l1*np.cos(-self.q1)),int(300+self.l1*np.sin(-self.q1)))
        cv2.line(image, (300,300), line1, (255,0,0), 1)
        line2=(int(line1[0]+self.l2*np.cos(-self.q2)),int(line1[1]+self.l2*np.sin(-self.q2)))
        cv2.line(image, line1, line2, (0,255,0), 1)
        if(len(self.path)>200):
            path=self.path[len(self.path)-200:]
        else:
            path=self.path
        for center in path:
            cv2.circle(image, center, 1, (0,0,0),-1)
        return image
    
mybot=TwoR1()
while(True):
    cv2.namedWindow("window")
    mybot.step()
    mybot.render()