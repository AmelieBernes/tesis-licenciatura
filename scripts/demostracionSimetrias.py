import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import flechas_2D

colores_ame=['hotpink', 'rebeccapurple']

fig, axis= plt.subplots()

axis.set_xlim(-1.3,1.3)
axis.set_ylim(-1.3,1.3)
flechas_2D.dibujar_flechas_2d(fig, axis)


x=[-1,-1,1]
y=[1,-1,1]
plt.scatter(x,y, color='black')

#para dimensión n=7.
#X=np.arange(-1.2, 1.2, 0.05)
#x=[-1,-2/3, -1/3, 0, 1/3, 2/3, 1]
#for i in range(7):
#    plt.plot(X, X**i, color=colores_ame[i%2])
#    y=[t**i for t in x]
#    plt.scatter(x,y, color=colores_ame[i%2])
#
##for abs in x:
##    plt.plot(0*X+abs, X, color='gray', linestyle=':')
#
#
#
#y=[0,0,0,0,0,0,0]
#plt.scatter(x,y, marker='x', color='black')


#para dimensión n=4
X=np.arange(-1.2, 1.2, 0.05)
x=[-1,-1/3,1/3,1]
for i in range(4):
    plt.plot(X, X**i, color=colores_ame[i%2])
    y=[t**i for t in x]
    plt.scatter(x,y, color=colores_ame[i%2])

#for abs in x:
#    plt.plot(0*X+abs, X, color='gray', linestyle=':')



y=[0,0,0,0]
plt.scatter(x,y, marker='x', color='black')

plt.grid()
plt.show()

