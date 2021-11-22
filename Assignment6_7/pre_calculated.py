import numpy as np
import matplotlib.pyplot as plt
from time import time
import pickle

with open('data\\_6a.pkl', 'rb') as file:
    xs_d_a,ys_d_a,xs_a,ys_a,time_a = pickle.load(file)
with open('data\\_6b.pkl', 'rb') as file:
    xs_d_b,ys_d_b,xs_b,ys_b,time_b = pickle.load(file)
with open('data\\_6c.pkl', 'rb') as file:
    xs_d_c,ys_d_c,xs_c,ys_c,time_c = pickle.load(file)
with open('data\\_6d.pkl', 'rb') as file:
    xs_d_d,ys_d_d,xs_d,ys_d,time_d = pickle.load(file)

plt.plot(time_a,ys_d_a,"black")
plt.plot(time_a,ys_a,"red")
plt.plot(time_b,ys_b,"cyan")
plt.plot(time_c,ys_c,"yellow")
plt.plot(time_d,ys_d,"blue")

# plt.plot(xs_d_a,ys_d_a,"black")
# plt.plot(xs_a,ys_a,"red")
# plt.plot(xs_b,ys_b,"cyan")
# plt.plot(xs_c,ys_c,"yellow")
# plt.plot(xs_d,ys_d,"blue")

plt.title('y vs time')
plt.xlabel('time')
plt.ylabel('y') 
# plt.legend(['x desired','xa','xb','xc','xd'])
# plt.legend(['y desired','ya','yb','yc','yd'])
plt.legend(['desired','a','b','c','d'])
# plt.xlim(0.399,0.401)
# plt.ylim(0,0.07)
plt.show()