import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import flechas_2D

colores_ame=['hotpink', 'rebeccapurple']

fig, axis= plt.subplots(1,2)

#--------------- Inicializamos las mallas. 
#--------------- Datos necesiarios: dimensión, punto inicial de P, paso de la malla P. 
n=11 #cantidad de puntos (dimensión)
t0=-2 #punto inicial de la malla P
h= 0.5 #paso de la malla P
malla_P=[]
malla_Pn=[]
for i in range(n):
    malla_Pn.append(i)

for i in range(n):
    malla_P.append(h*malla_Pn[i]+t0)

ceros=[0]*n
#------------------------------------
X=np.arange(-2.2,3.2,0.05)
Y=np.arange(-0.4,10.4,0.05)
valores_funcion=[t**3-3*t+1 for t in malla_P]
valores_phi=[h*t+t0 for t in malla_Pn]

axis[0].plot(X, X**3-3*X+1, color=colores_ame[0], label='$f(t)=t^{3}-3t+1$')
axis[0].scatter(malla_P, valores_funcion, color=colores_ame[0])
axis[0].scatter(malla_P,ceros, color='black', marker='|')

axis[1].plot(Y, h*Y+t0, color=colores_ame[1], label='$\phi=0.5t-2$')
axis[1].scatter(malla_Pn, valores_phi, color=colores_ame[1])
axis[1].scatter(malla_Pn,ceros, color='black', marker='|')



for i in range(2):
    flechas_2D.dibujar_flechas_2d(fig, axis[i])
    axis[i].grid(True)
    axis[i].legend()
plt.show()
