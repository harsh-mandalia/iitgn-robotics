from render import Renderer
import cv2
import numpy as np
from scipy.optimize import fsolve
import time

class TwoR1(Renderer):
    def __init__(self):
        super().__init__()
        self.l1=100
        self.l2=100
        self.m1=1
        self.m2=1
        self.q1=np.pi/6
        self.q2=np.pi/2
        self.x=300
        self.y=300
        self.g=9.81
        self.speed=200
        
        self._q1=self.q1
        self._q2=self.q2
        self.have_path=False
        self.Ldown=False
        self.done=False
        self.path=[]
        self.qpath=[]
        self.fast_qpath=[]
        self.info={}
    
    def equations(self, p):
        q1,q2=p
        eq1=self.l1*np.cos(q1)+self.l2*np.cos(q2)-self.x
        eq2=self.l1*np.sin(q1)+self.l2*np.sin(q2)-self.y
        return [eq1,eq2]

    def callback(self, event,x,y,flags,param):
        if event == cv2.EVENT_LBUTTONDOWN:
            self.Ldown=True

        elif self.done==False and event == cv2.EVENT_LBUTTONUP:
            self.have_path=True

            for i in range(len(self.path)):
                x=self.path[i][0]
                y=self.path[i][1]
                t=self.path[i][2]
                self.x=x
                self.y=y

                if i>2 and (self.x**2+self.y**2)<=(self.l1+self.l2)**2:
                    q1_0=self.qpath[i-2][0]
                    q1_1=self.qpath[i-1][0]
                    q1_2=self.qpath[i][0]
                    q2_0=self.qpath[i-2][1]
                    q2_1=self.qpath[i-1][1]
                    q2_2=self.qpath[i][1]
                    t0=self.path[i-2][2]
                    t1=self.path[i-1][2]
                    t2=self.path[i][2]

                    dt=(t2-t0)/2

                    q1_dot=(q1_2-q1_1)/dt
                    q2_dot=(q2_2-q2_1)/dt
                    q1_dotdot=(q1_0+q1_2-2*q1_1)/(2*t*t)
                    q2_dotdot=(q2_0+q2_2-2*q2_1)/(2*t*t)
                    
                    q1=q1_2
                    q2=q2_2
                    q1_1=q1_dot
                    q1_2=q1_dotdot
                    q2_1=q2_dot
                    q2_2=q2_dotdot

                    m1, m2, l1, l2, g = self.m1, self.m2, self.l1, self.l2, self.g
                    a = -0.5*m2*l1*l2*q2_1*(q2_1-q1_1)*np.sin(q2-q1) - 0.5*m2*l1*l2*q1_1*q2_1*np.sin(q2-q1) + 0.5*m1*g*l1*np.cos(q1) + m2*g*l1*np.cos(q1)
                    b = -0.5*m2*l1*l2*q1_1*(q2_1-q1_1)*np.sin(q2-q1) + 0.5*m2*l1*l2*q1_1*q2_1*np.sin(q2-q1) + 0.5*m2*g*l2*np.cos(q2)
                    c = 1/3*m1*l1**2 + m2*l1**2
                    d = 0.5*m2*l1*l2*np.cos(q2-q1)
                    e = 1/12*m2*l1**2 + 1/4*m2*l2**2
                    f = 0.5*m2*l1*l2*np.cos(q2-q1)

                    tau1=c*q1_2+d*q2_2+a
                    tau2=f*q1_2+e*q2_2+b
                    self.info={"tau1":tau1,"tau2":tau2}

                    self.q1=q1
                    self.q2=q2
                self.render()

            self.done=True

        elif self.Ldown and self.have_path==False and event == cv2.EVENT_MOUSEMOVE:
            # print(a)
            self.x=x-300
            self.y=-y+300
            if(self.x**2+self.y**2)<=(self.l1+self.l2)**2:
                self.path.append((self.x,self.y,time.time()))
                self._q1,self._q2=fsolve(self.equations, (self._q1, self._q2))
                self.qpath.append((self._q1,self._q2))
            else:
                print("Path out of workspace, run againg with proper path :(")
                exit()

    def step(self,steps=10):
        cv2.setMouseCallback("window",self.callback)        

    def getInfo(self):
        return self.info

    def draw(self, image):
        line1=(int(300+self.l1*np.cos(-self.q1)),int(300+self.l1*np.sin(-self.q1)))
        cv2.line(image, (300,300), line1, (255,0,0), 1)
        line2=(int(line1[0]+self.l2*np.cos(-self.q2)),int(line1[1]+self.l2*np.sin(-self.q2)))
        cv2.line(image, line1, line2, (0,255,0), 1)
        for x,y,t in self.path:
            cv2.circle(image, (x+300,-y+300), 1, (0,0,0),-1)
        return image
    
mybot=TwoR1()
print("press \"q\" to exit")
while(True):
    cv2.namedWindow("window")
    mybot.step()
    mybot.render()
    if(mybot.done):
        cv2.waitKey(0)
        break