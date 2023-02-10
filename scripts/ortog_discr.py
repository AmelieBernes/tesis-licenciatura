import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import flechas_2D

#-- Definición de objetos

colores_ame=['goldenrod', 'mediumpurple', 'darkturquoise', 'hotpink']

fig, axis=plt.subplots(2,2)
X=np.arange(-1.2,1.2,0.05)

#def f0(x):
#    return 0*x+1
#
#def f1(x):
#    return x
#
#def f2(x):
#    return x**2
#
#def f3(x):
#    return x**3

def p0(x):
    return 0*x+1/math.sqrt(2)

def p1(x):
    return math.sqrt(3/2)*x

def p2(x):
    return math.sqrt(5/8)*(3*x**2-1)

def p3(x):
    return (5/2)*math.sqrt(7/2)*(x**3-(3/5)*x)

#funciones_f=[f0,f1,f2,f3]
funciones_p=[p0,p1,p2,p3]
mallaP=[-1,-1/3,1/3,1]


#--------------------------------------
#Graficar

#Para dibujar a los 4 primeros polinomios ortonormales de Legendre
#for i in range(2):
#    for j in range(2):
#        axis[i][j].plot(X, funciones_p[2*i+j](X), color=colores_ame[2*i+j])
#        axis[i][j].set_title('Grado '+str(2*i+j))
#        axis[i][j].set_xlim([-1.2,1.2])
#        axis[i][j].set_ylim([-2.2,2.2])
#        axis[i][j].grid("True")
#        axis[i][j].axhline(y=0, color='gray')
#        axis[i][j].axvline(x=0, color='gray')
#fig.suptitle('Ortonormalización en $L^{2}[-1,1]$')


#obtengamos las discretizaciones:
discretizaciones=[]
for p in funciones_p:
    discretizaciones.append([p(t) for t in mallaP])

for i in range(2):
    for j in range(2):
        axis[i][j].scatter(mallaP, discretizaciones[2*i+j], color=colores_ame[2*i+j], s=100)
        axis[i][j].set_title('Grado '+str(2*i+j))
        axis[i][j].set_xlim([-1.2,1.2])
        axis[i][j].set_ylim([-2.2,2.2])
        axis[i][j].grid("True")
        axis[i][j].axhline(y=0, color='gray')
        axis[i][j].axvline(x=0, color='gray')
fig.suptitle('Discretización')

plt.show()
