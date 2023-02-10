import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import flechas_2D
import legendre
import proyecciones as pr #para dibujar parábolas

colores_ame=['hotpink', 'gray']
eje=[0,1,2,3]
X=np.arange(0,3, 0.05)
plt.xlim([-0.1,3.1])
plt.ylim([-2.2,2.2])


def constante():
    #Constante
    plt.scatter(eje, legendre.BON_L[4][0], color=colores_ame[0], s=100)
    plt.plot(eje, legendre.BON_L[4][0], color=colores_ame[0], linestyle='dotted')

    plt.plot([-0.5,3.5], [0,0], color=colores_ame[1])
    plt.plot([0,0], [-2.5,3.5], color=colores_ame[1])
    plt.title("Gráfica de $\mathcal{L}^{4,0}$", fontsize=16)
    plt.grid()
    plt.show()



def lineal():
    #Constante
    plt.scatter(eje, legendre.BON_L[4][0], color=colores_ame[1], s=100)
    plt.plot(eje, legendre.BON_L[4][0], color=colores_ame[1], linestyle='dotted')
    #Lineal
    plt.scatter(eje, legendre.BON_L[4][1], color=colores_ame[0], s=100)
    plt.plot(eje, legendre.BON_L[4][1], color=colores_ame[0], linestyle='dotted')

    plt.plot([-0.5,3.5], [0,0], color=colores_ame[1])
    plt.plot([0,0], [-2.5,3.5], color=colores_ame[1])
    plt.title("Gráfica de $\mathcal{L}^{4,0}$ y $\mathcal{L}^{4,1}$", fontsize=16)
    plt.grid()
    plt.show()

def cuadratico():
    fig, axis=plt.subplots(1,2)
    axis[0].set_title("Gráficas de $\mathcal{L}^{4,0}, \mathcal{L}^{4,1}$ y algunos elementos de $W_{4,2}$",fontsize=12)
    axis[1].set_title("Gráficas de $\mathcal{L}^{4,0}, \mathcal{L}^{4,1} y \mathcal{L}^{4,2}$", fontsize=12) 
    #Constante
    for i in range(2):
        axis[i].scatter(eje, legendre.BON_L[4][0], color=colores_ame[1], s=100)
        axis[i].plot(eje, legendre.BON_L[4][0], color=colores_ame[1], linestyle='dotted')
    #Lineal
    for i in range(2):
        axis[i].scatter(eje, legendre.BON_L[4][1], color=colores_ame[1], s=100)
        axis[i].plot(eje, legendre.BON_L[4][1], color=colores_ame[1], linestyle='dotted')
    #Cuadráticas equivocadas
    axis[0].scatter(eje, [0.5*t**2-1.5*t-1.3 for t in range(4)] , color=colores_ame[0], s=100)
    axis[0].plot(X, 0.5*X**2-1.5*X-1.2, color=colores_ame[0], linestyle='dotted')

    axis[0].scatter(eje, [-(0.7*t-1)**2+1.7 for t in range(4)] , color=colores_ame[0], s=100)
    axis[0].plot(X, -(0.7*X-1)**2+1.7, color=colores_ame[0], linestyle='dotted')

    #Cuadrática Legendre
    axis[1].scatter(eje, legendre.BON_L[4][2], color=colores_ame[0], s=100)
    #coeficientes de la parábola que pasa por los puntos anteriores
    coef=pr.parab3Puntos([legendre.BON_L[4][2][i] for i in range(4) ])
    axis[1].plot(X,coef[0]*X**2+coef[1]*X+coef[2],color=colores_ame[0], linestyle='dotted')
    
    for i in range(2):
        axis[i].plot([-0.5,3.5], [0,0], color=colores_ame[1])
        axis[i].plot([0,0], [-2.5,3.5], color=colores_ame[1])
        axis[i].grid()
    plt.xlim([-0.1,3.1])
    plt.ylim([-2.2,2.2])
    return plt.show()


def cubico():
    fig, axis=plt.subplots(1,2)
    #Constante Legendre
    axis[0].set_title("Gráficas de $\mathcal{L}^{4,0}, \mathcal{L}^{4,1}, \mathcal{L}^{4,2}$ y un elemento de $W_{4,3}$",fontsize=12)
    axis[1].set_title("Gráficas de $\mathcal{L}^{4,0}, \mathcal{L}^{4,1}, \mathcal{L}^{4,2}$ y $\mathcal{L}^{4,3}$", fontsize=12) 
    for i in range(2):
        axis[i].scatter(eje, legendre.BON_L[4][0], color=colores_ame[1], s=100)
        axis[i].plot(eje, legendre.BON_L[4][0], color=colores_ame[1], linestyle='dotted')
    #Lineal Legendre
    for i in range(2):
        axis[i].scatter(eje, legendre.BON_L[4][1], color=colores_ame[1], s=100)
        axis[i].plot(eje, legendre.BON_L[4][1], color=colores_ame[1], linestyle='dotted')
    #Cuadrática Legendre
    coef=pr.parab3Puntos([legendre.BON_L[4][2][j] for j in range(4) ])
    for i in range(2):
        axis[i].scatter(eje, legendre.BON_L[4][2], color=colores_ame[1], s=100)
        #coeficientes de la parábola que pasa por los puntos anteriores
        axis[i].plot(X,coef[0]*X**2+coef[1]*X+coef[2],color=colores_ame[1], linestyle='dotted')
   
    #Cúbica equivocada
    axis[0].scatter(eje, [0.17*t**3-0.7*t**2+2*t-2 for t in range(4)],color=colores_ame[0], s=100)
    axis[0].plot(X, 0.17*X**3-0.7*X**2+2*X-2, color=colores_ame[0], linestyle='dotted')

    #Cúbica Legendre
    axis[1].scatter(eje, legendre.BON_L[4][3],color=colores_ame[0], s=100)
    axis[1].plot(X, 0.745*X**3-3.35*X**2+3.5*X-0.22, color=colores_ame[0], linestyle='dotted')


    for i in range(2):
        axis[i].plot([-0.5,3.5], [0,0], color=colores_ame[1])
        axis[i].plot([0,0], [-2.5,3.5], color=colores_ame[1])
        axis[i].grid()
    plt.xlim([-0.1,3.1])
    plt.ylim([-2.2,2.2])
    return plt.show()


#constante()
#lineal()
#cuadratico()
cubico()
