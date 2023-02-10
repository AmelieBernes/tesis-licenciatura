import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import flechas_2D
import legendre


#-- Definición de objetos

colores_ame=['goldenrod', 'mediumpurple', 'darkturquoise', 'hotpink']

fig, axis=plt.subplots(2,2)
X=np.arange(-1.2,1.2,0.05)

def f0(x):
    return 0*x+1

def f1(x):
    return x

def f2(x):
    return x**2

def f3(x):
    return x**3

#funciones_f=[f0,f1,f2,f3]
mallaP=[-1,-1/3,1/3,1]

#Discretización------------------------------
##obtengamos las discretizaciones:
#discretizaciones=[]
#for f in funciones_f:
#    discretizaciones.append([f(p) for p in mallaP])
#
#for i in range(2):
#    for j in range(2):
#        axis[i][j].scatter(mallaP, discretizaciones[2*i+j], color=colores_ame[2*i+j], s=100)
#        axis[i][j].set_title('Grado '+str(2*i+j))
#        axis[i][j].set_xlim([-1.2,1.2])
#        axis[i][j].set_ylim([-1.2,1.2])
#        axis[i][j].grid("True")
#        axis[i][j].axhline(y=0, color='gray')
#        axis[i][j].axvline(x=0, color='gray')
#
#
#fig.suptitle('Discretización')



#Ortonormalización en Rn-------------------
for i in range(2):
    for j in range(2):
        axis[i][j].scatter(mallaP, legendre.BON_L[4][2*i+j], color=colores_ame[2*i+j], s=100)
        axis[i][j].set_title('Grado '+str(2*i+j))
        axis[i][j].set_xlim([-1.2,1.2])
        axis[i][j].set_ylim([-1.2,1.2])
        axis[i][j].grid("True")
        axis[i][j].axhline(y=0, color='gray')
        axis[i][j].axvline(x=0, color='gray')


fig.suptitle('Ortonormalización en $\mathbb{R}^{n}$')

plt.show()
