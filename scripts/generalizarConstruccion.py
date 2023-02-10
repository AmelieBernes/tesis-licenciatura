import numpy as np
import matplotlib.pyplot as plt
import pylab #Para usar LaTeX en captions
import math
import flechas_2D

colores_ame=['goldenrod','hotpink','mediumpurple']

fig, axis= plt.subplots(1,2)
axis[0].grid(True)
axis[1].grid(True)
dominio=[0,1,2]
X=np.arange(-1.2,2.2,0.05)

def f0(x):
    return 1
def f1(x):
    return x
def f2(x):
    return x**2

def g0(x):
    return 3
def g1(x):
    return 0.5*x+1
def g2(x):
    return x**2+2*x+3

funciones_f=[f0,f1,f2]
funciones_g=[g0,g1,g2] 
#--------------------------------------------
axis[0].plot(X,0*X+1, color=colores_ame[0], label='$f_{0}$')
axis[0].plot(X,X, color=colores_ame[1],label='$f_{1}$')
axis[0].plot(X,X**2, color=colores_ame[2],label='$f_{2}$')

axis[1].plot(X, 0*X+3, color=colores_ame[0],label='$g_{0}$')
axis[1].plot(X, 0.5*X+1, color=colores_ame[1],label='$g_{1}$')
axis[1].plot(X, X**2+2*X+3, color=colores_ame[2],label='$g_{2}$')

#--------------------------------------------
puntos1=[]
for i in range(3):
    lista=[funciones_f[i](x) for x in dominio]
    puntos1.append(lista)

puntos2=[]
for i in range(3):
    lista=[funciones_g[i](x) for x in dominio]
    puntos2.append(lista)
for i in range(3):
    axis[0].scatter(dominio, puntos1[i],color=colores_ame[i])
    axis[1].scatter(dominio, puntos2[i],color=colores_ame[i])

#--------------------------------------------

flechas_2D.dibujar_flechas_2d(fig, axis[0])
flechas_2D.dibujar_flechas_2d(fig, axis[1])

axis[0].legend()
axis[1].legend()
plt.show()
