import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import flechas_2D

colores_ame=['hotpink','rebeccapurple','goldenrod','gray']

def f(x):
    return 0*x+1

def g(x):
    return 1-5*x*(x-0.4)*(x-0.8)*(x-1)*(x-2)

def h(x):
    return 1+3*x*(x-1)*(x-2)

def i(x):
    return 1-3*x*(x-1)*(x-2)


def p(x):
    return 1+x*(x-1)*(x-2.5)

#def q(x):
#    return -0.5*x**2+0.5*x+1 

def q(x):
    return -0.5*x**2 + x + 0.625

fig, axis= plt.subplots(1,2)
X=np.arange(-0.2,2.2,0.05)
dominio=[0,1,2]

senial=[1,1,1]

axis[0].plot(X,g(X), color=colores_ame[1], linestyle=':')
axis[0].plot(X,h(X), color=colores_ame[2], linestyle=':')

axis[0].plot(X,f(X), color=colores_ame[0])
axis[0].scatter(dominio,senial, color=colores_ame[0], label='$x=(1,1,1) \in \mathbb{R}^{3}$')

#--------------------------------------
senial=[1,1,0]
axis[1].plot(X,p(X) , color=colores_ame[1], linestyle=':')
axis[1].scatter(dominio,senial, color=colores_ame[0], label='$x=(1,1,0) \in \mathbb{R}^{3}$')

X=np.arange(-0.2,2.7,0.05)
axis[1].plot(X,q(X), color=colores_ame[0])
dominio= [0.5, 1.5, 2.5]
axis[1].scatter(dominio,senial, color=colores_ame[0], s = 100, label='$x=(1,1,0) \in \mathbb{R}^{3}$')




for i in range(2):
    flechas_2D.dibujar_flechas_2d(fig, axis[i])
    axis[i].grid(True)
    axis[i].legend()
plt.show()
