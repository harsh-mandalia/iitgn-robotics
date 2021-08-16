from render import Renderer
import cv2
import numpy as np

class TwoR1(Renderer):
    def __init__(self):
        super().__init__()
        self.l1=100
        self.l2=100
        self.q1=np.pi/6
        self.q2=np.pi/2
        self.x=300
        self.y=300
        self.points=[]
    
    def equations(self, p):
        q1,q2=p
        eq1=self.l1*np.cos(q1)+self.l2*np.cos(q2)-self.x
        eq2=self.l1*np.sin(q1)+self.l2*np.sin(q2)-self.y
        return [eq1,eq2]

    def trace(self,speed=3):
        min=35
        max=145
        for q1 in range(min,max+1, speed):
            for q2 in range(min,max+1, speed):
                self.q1=q1*np.pi/180
                self.q2=q2*np.pi/180
                self.render()


    def getInfo(self):
        info={}
        return info

    def draw(self, image):
        line1=(int(300+self.l1*np.cos(-self.q1)),int(300+self.l1*np.sin(-self.q1)))
        cv2.line(image, (300,300), line1, (255,0,0), 1)
        line2=(int(line1[0]+self.l2*np.cos(-self.q2)),int(line1[1]+self.l2*np.sin(-self.q2)))
        cv2.line(image, line1, line2, (0,255,0), 1)
        self.points.append(line2)
        for i in self.points:
            cv2.circle(image,i,1,(0,0,0),1)
        return image
    
mybot=TwoR1()
mybot.trace(5)
cv2.waitKey(5000)