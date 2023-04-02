#Script en el que se guarda el código de todas las figuras incluidas en la tesis
#que fueron programadas en Python.


# Módulos

import numpy as np
import matplotlib.pyplot as plt
import pylab #Para poder usar LaTeX? TODO creo que no 
import math
import flechas_2D

import random
from random import uniform


#Código para cambiar el tipo de fuente

plt.style.use('seaborn-v0_8-poster') 
#plt.style.use('seaborn-v0_8-pastel') 
params = {"ytick.color" : "black",
          "xtick.color" : "black",
          "axes.labelcolor" : "black",
          "axes.edgecolor" : "black",
          "text.usetex" : True,
          "font.family" : "serif",
          "font.serif" : ["Computer Modern Serif"]}
plt.rcParams.update(params)

colores=['hotpink', 'rebeccapurple','goldenrod','gray']

def figura_introduccion():
    X=np.arange(0, 30.5, 0.05)
    dominio=[]
    for i in range(30):
        dominio.append(i)
   
   #usamos la función random.uniform para agregar ruido a la figura.
    mediciones_parabola=[]
    for i in range(30):
        x=dominio[i]
        mediciones_parabola.append((0.23*x-2)**2-10+random.uniform(-0.5, 0.5))
    
    mediciones_recta=[]
    for i in range(30):
        x=dominio[i]
        mediciones_recta.append(-0.7*x+5+random.uniform(-0.5, 0.5))
    
    plt.scatter(dominio, mediciones_parabola, color='hotpink')
    plt.scatter(dominio, mediciones_recta, color='hotpink', marker='x')
    plt.grid()
    return plt.show()


def figura_discretizacionIntegral():

    fig, axis=plt.subplots(2,2)

    X=np.arange(-1.2, 1.2, 0.05)
    x=[-1,-0.5,0,0.5,1]
    for i in range(2):
        for j in range(2):
            axis[i][j].set_xlim(-1.3,1.3)
            axis[i][j].set_ylim(-1.3,1.3)
            flechas_2D.dibujar_flechas_2d(fig, axis[i,j])
            axis[i,j].grid(True)
            axis[i,j].plot(X, X**(2*i+j), color=colores[0])
            y=[t**(2*i+j) for t in x]
            axis[i,j].scatter(x, y, color=colores[0])
            k=2*i+j
            axis[i,j].set_title("Grado "+str(k)+", $u_{}$".format(str(k))+"$\mapsto \mathcal{L}^{4,grado}.$".replace("grado", str(k)), fontname="Computer Modern Serif", size=16)
    return plt.show()


def figura_discretizacionPuntual():
    fig, axis=plt.subplots(2,2)
    X=np.arange(-1.2, 2.2, 0.05)
    def f1(x):
        """
        x un objeto de tipo float
        """
        return 0*x+0.3
    def f2(x):
        """
        x un objeto de tipo float
        """
        return -0.4*x+1
    def f3(x):
        """
        x un objeto de tipo float
        """
        return 0.5*x**2-0.5*x
    def f4(x):
        """
        x un objeto de tipo float
        """
        return -x**3+2*x**2+x-2

    funciones=[f1, f2,f3,f4]
    x=[-1,0,1,2]
    signos=["+","-","+","-"]
    for i in range(2):
        for j in range(2):
            axis[i][j].set_xlim(-1.3,2.3)
            axis[i][j].set_ylim(-2.3,2.3)
            flechas_2D.dibujar_flechas_2d(fig, axis[i,j])
            axis[i,j].grid(True)
            axis[i,j].plot(X, funciones[2*i+j](X), color=colores[0])
            y=[funciones[2*i+j](t) for t in x]
            axis[i,j].scatter(x, y, color=colores[0], s=100)
            k=2*i+j
            axis[i,j].set_title("Grado "+str(k)+", $u_{}$".format(str(k))+"$\mapsto$"+signos[k]+"$\mathcal{L}^{4,grado}.$".replace("grado", str(k)), fontname="Computer Modern Serif", size=16)
    return plt.show()

def figura_cambioDeMalla():
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
    
    axis[0].plot(X, X**3-3*X+1, color=colores[0], label='$f(t)=t^{3}-3t+1$')
    axis[0].scatter(malla_P, valores_funcion, color=colores[0])
    axis[0].scatter(malla_P,ceros, color='black', marker='|')
    
    axis[1].plot(Y, h*Y+t0, color=colores[1], label='$\phi=0.5t-2$')
    axis[1].scatter(malla_Pn, valores_phi, color=colores[1])
    axis[1].scatter(malla_Pn,ceros, color='black', marker='|')
    
    for i in range(2):
        flechas_2D.dibujar_flechas_2d(fig, axis[i])
        axis[i].grid(True)
        axis[i].legend()

    return plt.show()

def figura_coseno():

    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)
    
    X=np.linspace(0,3.15,100)
    plt.plot(X, np.cos(X), color="gray", linestyle="dotted", label="$y=cos(\\theta)$")
    
    plt.scatter(0, 0, marker="o", color="gray", s=200, label="$v$ y $w$ son paralelos")
    plt.scatter(3.14, 0, marker="o", color="gray", s=200)
    plt.scatter(3.14/2, 0, marker="s", color="gray", s=200, label="$v$ y $w$ son perpendiculares")
    plt.scatter(0.7, np.cos(0.7), marker="*", color="gray", s=300, label="$\\frac{\langle v, w \\rangle}{||v|| \cdot ||w||}$")
    
    
    plt.text(0,0.1,"0", fontsize=17)
    plt.text(3.14,0.1,"$\pi$", fontsize=17)
    plt.text(3.14/2,0.1,"$\\frac{\pi}{2}$", fontsize=17)
    
    flechas_2D.dibujar_flechas_2d(fig, ax)
    
    plt.grid()
    plt.legend()
    return plt.show()


def figura_demSimetrias(n=7):
    """
    'n' sólo puede tomar el valor 7 o 4.
    """
    fig, axis= plt.subplots()
    
    axis.set_xlim(-1.3,1.3)
    axis.set_ylim(-1.3,1.3)
    flechas_2D.dibujar_flechas_2d(fig, axis)
    
    
    x=[-1,-1,1]
    y=[1,-1,1]
    plt.scatter(x,y, color='black')
    
    if n==7:
        X=np.arange(-1.2, 1.2, 0.05)
        x=[-1,-2/3, -1/3, 0, 1/3, 2/3, 1]
        for i in range(7):
            plt.plot(X, X**i, color=colores[i%2])
            y=[t**i for t in x]
            plt.scatter(x,y, color=colores[i%2])
        
        y=[0,0,0,0,0,0,0]
        plt.scatter(x,y, marker='x', color='black')
    
    else: #n=4
        X=np.arange(-1.2, 1.2, 0.05)
        x=[-1,-1/3,1/3,1]
        for i in range(4):
            plt.plot(X, X**i, color=colores[i%2])
            y=[t**i for t in x]
            plt.scatter(x,y, color=colores[i%2])
        
        y=[0,0,0,0]
        plt.scatter(x,y, marker='x', color='black')
    
    plt.grid()
    return plt.show()

def figura_discret_ortog(disc=True):
    """
    disc==True: discretizacion
    disc==False: ortogonalización
    """
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
    
    mallaP=[-1,-1/3,1/3,1]
   
    if disc==True:
        fig.suptitle('Discretización')
        discretizaciones=[]
        for f in funciones_f:
            discretizaciones.append([f(p) for p in mallaP])
        
        for i in range(2):
            for j in range(2):
                axis[i][j].scatter(mallaP, discretizaciones[2*i+j], color=colores_ame[2*i+j], s=100)
                axis[i][j].set_title('Grado '+str(2*i+j))
                axis[i][j].set_xlim([-1.2,1.2])
                axis[i][j].set_ylim([-1.2,1.2])
                axis[i][j].grid("True")
                axis[i][j].axhline(y=0, color='gray')
                axis[i][j].axvline(x=0, color='gray')
    
    
    if disc==False
        fig.suptitle('Ortonormalización en $\mathbb{R}^{n}$')
        legendre_4=legendre.calculo_base(4) #calculamos la base de legendre discreta de dimensión 4
        for i in range(2):
            for j in range(2):
                axis[i][j].scatter(mallaP, legendre_4[2*i+j], color=colores_ame[2*i+j], s=100)
                axis[i][j].set_title('Grado '+str(2*i+j))
                axis[i][j].set_xlim([-1.2,1.2])
                axis[i][j].set_ylim([-1.2,1.2])
                axis[i][j].grid("True")
                axis[i][j].axhline(y=0, color='gray')
                axis[i][j].axvline(x=0, color='gray')

    return plt.show()


if __name__=='__main__':
    #figura_discretizacionPuntual()
    #figura_cambioDeMalla()
    #figura_introduccion()
    #figura_coseno()
    #figura_demSimetrias(4)
    #figura_discret_ortog(True) #no reconoce el True/False...
