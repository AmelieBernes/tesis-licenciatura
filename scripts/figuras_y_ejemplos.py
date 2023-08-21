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

import base_legendreDiscreta as legendre
import proyecciones as pr


#Código para cambiar el tipo de fuente
#TODO no sé por qué no siempre funciona, a veces da muchos problemas y no deja que se genere la imagen!

plt.style.use('seaborn-v0_8-poster') 
##plt.style.use('seaborn-v0_8-pastel') 
#params = {"ytick.color" : "black",
#          "xtick.color" : "black",
#          "axes.labelcolor" : "black",
#          "axes.edgecolor" : "black",
#          "text.usetex" : True,
#          "font.family" : "serif",
#          "font.serif" : ["Computer Modern Serif"]}
#plt.rcParams.update(params)


def formato_axis(axis):
  """
  Agregando los elementos que me gustan a un axis
  """
  axis.axhline(y=0, color='gray')
  axis.axvline(x=0, color='gray')
  axis.grid(True)
  axis.legend()

colores=['hotpink', 'rebeccapurple','goldenrod','gray']


def figura_introduccion():
    fig, axis = plt.subplots()
    X=np.arange(0, 30.5, 0.05)
    dominio=[i for i in range(30)]
   
    mediciones_parabola=[]
    for i in range(30):
        x=dominio[i]
        mediciones_parabola.append((0.23*x-2)**2-10+random.uniform(-0.5, 0.5))
    
    mediciones_recta=[]
    for i in range(30):
        x=dominio[i]
        mediciones_recta.append(-0.7*x+5+random.uniform(-0.5, 0.5))
    
    axis.scatter(dominio, mediciones_parabola, color='hotpink')
    axis.scatter(dominio, mediciones_recta, color='hotpink', marker='x')
    axis.grid()
    return plt.show()


def figura_rectaParabola():
    fig, axis = plt.subplots(2,2)
    X = np.arange(-0.5, 15.5, 0.05)
    def recta(t):
        return 2*t+1
    def parabola(t):
        return (0.5*t-3)**2

    dominio = [t for t in range(15)]

    cte_disc = [20 for t in dominio]
    cte_discNoise = [20 + random.uniform(-1.5,1.5) for t in dominio]

    recta_disc = [recta(t) for t in dominio]
    recta_discNoise = [t + random.uniform(-1.5,1.5) for t in recta_disc]

    parab_disc = [parabola(t) for t in dominio]
    parab_discNoise = [t + random.uniform(-1.5,1.5) for t in parab_disc]

    axis[0,0].scatter(dominio, recta_disc, color = 'hotpink')
    axis[0,0].scatter(dominio, cte_disc, color = 'hotpink', marker = 'x')
    axis[0,0].plot(X, recta(X), linestyle = 'dotted', color = 'hotpink')
    axis[0,0].set_title("Señales afines")

    axis[0,1].scatter(dominio, recta_discNoise, color = 'hotpink')
    axis[0,1].scatter(dominio, cte_discNoise, color = 'hotpink', marker = 'x')
    axis[0,1].set_title("Señales casi afines")

    axis[1,0].plot(X, parabola(X), linestyle = 'dotted', color = 'hotpink')
    axis[1,0].scatter(dominio, parab_disc, color = 'hotpink')
    axis[1,0].set_title("Señal cuadrática")

    axis[1,1].scatter(dominio, parab_discNoise, color = 'hotpink')
    axis[1,1].set_title("Señal casi cuadrática")

    for i in range(2):
        for j in range(2):
            axis[i, j].grid()

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


def figura_generalizacion_operadorOmega():
    def g0(t):
        return 0*t+3
    def g1(t):
        return 0.5*t + 1
    def g2(t):
        return t**2 + 2*t +3

    fig, axis = plt.subplots(1,2)
    X = np.arange(-1.3, 2.3, 0.001)
    P = [0, 1, 2]
    axis[0].plot(X, 0*X+1, color = 'mediumturquoise', label = '$f_{0}$')
    axis[0].scatter(P, [1,1,1], color = 'mediumturquoise')
    axis[0].plot(X, X, color = 'rebeccapurple', label = '$f_{1}$')
    axis[0].scatter(P, [0, 1, 2], color = 'rebeccapurple')
    axis[0].plot(X, X**2, color = 'goldenrod', label = '$f_{2}$')
    axis[0].scatter(P, [0, 1, 4], color = 'goldenrod')

    X = np.arange(-2.3, 1.3, 0.001)
    P = [-2, -0.5, 1]
    axis[1].plot(X, g0(X), color = 'mediumturquoise', label = '$g_{0}$')
    axis[1].scatter(P, [3, 3, 3], color = 'mediumturquoise')
    axis[1].plot(X, g1(X), color = 'rebeccapurple', label = '$g_{1}$')
    axis[1].scatter(P, [g1(p) for p in P], color = 'rebeccapurple')
    axis[1].plot(X, g2(X), color = 'goldenrod', label = '$g_{2}$')
    axis[1].scatter(P, [g2(p) for p in P], color = 'goldenrod')

    for i in range(2):
        axis[i].grid()
        axis[i].legend()
        axis[i].axhline(y = 0, color = 'gray')
        axis[i].axvline(x = 0, color = 'gray')

    return plt.show()


def figura_coseno():

    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)
    
    X=np.linspace(0,3.15,100)
    plt.plot(X, np.cos(X), color="black", linestyle="dotted", label="$y=cos(\\theta)$")
    
    plt.scatter(0, 1, marker="o", color="black", s=200, label="$v$ y $w$ son paralelos")
    plt.scatter(3.14, -1, marker="o", color="black", s=200)
    plt.scatter(3.14/2, 0, marker="s", color="black", s=200, label="$v$ y $w$ son perpendiculares")
    plt.scatter(0.7, np.cos(0.7), marker="*", color="black", s=300, label="$\\frac{\langle v, w \\rangle}{||v|| \cdot ||w||}$")
    
    
    plt.text(0,0.1,"0", fontsize=17)
    plt.text(3.14,0.1,"$\pi$", fontsize=17)
    plt.text(3.14/2,0.1,"$\\frac{\pi}{2}$", fontsize=17)
   
    ax.axhline(y=0, color = 'gray')
    ax.axvline(x=0, color = 'gray')
    
    plt.grid()
    plt.legend()
    return plt.show()


def figura_demSimetrias():
    colores = ['red', 'blue']

    fig, axis= plt.subplots(1,2)
    for i in range(2):
        axis[i].set_xlim(-1.3,1.3)  
        axis[i].set_ylim(-1.3,1.3)
        flechas_2D.dibujar_flechas_2d(fig, axis[i])
        x=[-1,-1,1]
        y=[1,-1,1]
        axis[i].scatter(x,y, color='black')
        axis[i].grid(True)
   
    X=np.arange(-1.2, 1.2, 0.05)
    x=[-1,-2/3, -1/3, 0, 1/3, 2/3, 1]
    for i in range(7):
        axis[1].plot(X, X**i, color = colores[i%2])
        y=[t**i for t in x]
        axis[1].scatter(x,y, color =  colores[i%2])
    
    y=[0,0,0,0,0,0,0]
    axis[1].scatter(x,y, marker='x', color = 'black')
    
    x=[-1,-1/3,1/3,1]
    for i in range(4):
        axis[0].plot(X, X**i, color = colores[i%2])
        y=[t**i for t in x]
        axis[0].scatter(x,y, color = colores[i%2])
    
    y=[0,0,0,0]
    axis[0].scatter(x,y, marker='x', color='black')
    
    return plt.show()

def figura_discret_ortog(disc=True):
    """
    disc==True: discretizacion
    disc==False: ortogonalización
    """
    colores = ['goldenrod', 'mediumpurple', 'darkturquoise', 'limegreen']
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
    funciones_f=[f0, f1, f2, f3]
   
    if disc==True:
        fig.suptitle('Discretización')
        discretizaciones=[]
        for f in funciones_f:
            discretizaciones.append([f(p) for p in mallaP])
        
        for i in range(2):
            for j in range(2):
                axis[i][j].scatter(mallaP, discretizaciones[2*i+j], color=colores[2*i+j], s=100)
                axis[i][j].set_title('Grado '+str(2*i+j))
                axis[i][j].set_xlim([-1.2,1.2])
                axis[i][j].set_ylim([-1.2,1.2])
                axis[i][j].grid("True")
                axis[i][j].axhline(y=0, color='gray')
                axis[i][j].axvline(x=0, color='gray')
    
    
    if disc==False:
        fig.suptitle('Ortonormalización') #en \IRn
        legendre_4=legendre.calculo_base(4) #calculamos la base de legendre discreta de dimensión 4
        for i in range(2):
            for j in range(2):
                axis[i][j].scatter(mallaP, legendre_4[2*i+j], color=colores[2*i+j], s=100)
                axis[i][j].set_title('Grado '+str(2*i+j))
                axis[i][j].set_xlim([-1.2,1.2])
                axis[i][j].set_ylim([-1.2,1.2])
                axis[i][j].grid("True")
                axis[i][j].axhline(y=0, color='gray')
                axis[i][j].axvline(x=0, color='gray')

    return plt.show()


def figura_ortogYoscil(k):
    """
    Figuras para la discusión de relación entre ortogonalidad y oscilaciones.
    'k' puede valer 0,1,2 o 3; indica qué gráfica se quiere.
    """
    colores_ame=['hotpink', 'gray']
    X=np.arange(0,3, 0.05)
    eje = [m for m in range(4)] 

    base_legendre = legendre.calculo_base(4) #Calculamos la base de Legendre discreta de dim. 4
    
    if k == 0: 
         #Constante
         plt.scatter(eje, base_legendre[0], color=colores_ame[0], s=200)
         plt.plot(eje, base_legendre[0], color=colores_ame[0], linestyle='dotted')
     
         plt.plot([-0.5,3.5], [0,0], color=colores_ame[1])
         plt.plot([0,0], [-2.5,3.5], color=colores_ame[1])
         plt.title("Gráfica de $\mathcal{L}^{4,0}$", fontsize=16)
         plt.ylim([-1,1])
         plt.grid()
         return plt.show()
     
     
    if k ==1: 
         #Constante
         plt.scatter(eje, base_legendre[0], color=colores_ame[1], s=200)
         plt.plot(eje, base_legendre[0], color=colores_ame[1], linestyle='dotted')
         #Lineal
         plt.scatter(eje, base_legendre[1], color=colores_ame[0], s=200)
         plt.plot(eje, base_legendre[1], color=colores_ame[0], linestyle='dotted')
     
         plt.plot([-0.5,3.5], [0,0], color=colores_ame[1])
         plt.plot([0,0], [-2.5,3.5], color=colores_ame[1])
         plt.title("Gráficas de $\mathcal{L}^{4,0}$ y $\mathcal{L}^{4,1}$", fontsize=16)
         plt.ylim([-1.2,1.2])
         plt.xlim([-.2,3.2])
         plt.grid()
         return plt.show()
     
    if k ==2:
        fig, axis=plt.subplots(1,2)
        axis[0].set_title("Gráficas de $\mathcal{L}^{4,0}, \mathcal{L}^{4,1}$ y algunos elementos de $W_{4,2}$",fontsize=12)
        axis[1].set_title("Gráficas de $\mathcal{L}^{4,0}, \mathcal{L}^{4,1} y \mathcal{L}^{4,2}$", fontsize=12) 
        #Constante
        for i in range(2):
            axis[i].scatter(eje, base_legendre[0], color=colores_ame[1], s=100)
            axis[i].plot(eje, base_legendre[0], color=colores_ame[1], linestyle='dotted')
        #Lineal
        for i in range(2):
            axis[i].scatter(eje, base_legendre[1], color=colores_ame[1], s=100)
            axis[i].plot(eje, base_legendre[1], color=colores_ame[1], linestyle='dotted')
        #Cuadráticas equivocadas
        axis[0].scatter(eje, [0.5*t**2-1.5*t-1.3 for t in range(4)] , color=colores_ame[0], s=100)
        axis[0].plot(X, 0.5*X**2-1.5*X-1.2, color=colores_ame[0], linestyle='dotted')
    
        axis[0].scatter(eje, [-(0.7*t-1)**2+1.7 for t in range(4)] , color=colores_ame[0], s=100)
        axis[0].plot(X, -(0.7*X-1)**2+1.7, color=colores_ame[0], linestyle='dotted')
    
        #Cuadrática Legendre
        axis[1].scatter(eje, base_legendre[2], color=colores_ame[0], s=100)
        #coeficientes de la parábola que pasa por los puntos anteriores
        coef=pr.parab3Puntos([base_legendre[2][i] for i in range(4) ])
        axis[1].plot(X,coef[0]*X**2+coef[1]*X+coef[2],color=colores_ame[0], linestyle='dotted')
        
        for i in range(2):
            axis[i].plot([-0.5,3.5], [0,0], color=colores_ame[1])
            axis[i].plot([0,0], [-2.5,3.5], color=colores_ame[1])
            axis[i].grid()
        plt.xlim([-0.1,3.1])
        plt.ylim([-2.2,2.2])
        return plt.show()
    
    
    if k == 3:
         fig, axis=plt.subplots(1,2)
         #Constante Legendre
         axis[0].set_title("Gráficas de $\mathcal{L}^{4,0}, \mathcal{L}^{4,1}, \mathcal{L}^{4,2}$ y un elemento de $W_{4,3}$",fontsize=12)
         axis[1].set_title("Gráficas de $\mathcal{L}^{4,0}, \mathcal{L}^{4,1}, \mathcal{L}^{4,2}$ y $\mathcal{L}^{4,3}$", fontsize=12) 
         for i in range(2):
             axis[i].scatter(eje, base_legendre[0], color=colores_ame[1], s=100)
             axis[i].plot(eje, base_legendre[0], color=colores_ame[1], linestyle='dotted')
         #Lineal Legendre
         for i in range(2):
             axis[i].scatter(eje, base_legendre[1], color=colores_ame[1], s=100)
             axis[i].plot(eje, base_legendre[1], color=colores_ame[1], linestyle='dotted')
         #Cuadrática Legendre
         coef=pr.parab3Puntos([base_legendre[2][j] for j in range(4) ])
         for i in range(2):
             axis[i].scatter(eje, base_legendre[2], color=colores_ame[1], s=100)
             #coeficientes de la parábola que pasa por los puntos anteriores
             axis[i].plot(X,coef[0]*X**2+coef[1]*X+coef[2],color=colores_ame[1], linestyle='dotted')
        
         #Cúbica equivocada
         axis[0].scatter(eje, [0.17*t**3-0.7*t**2+2*t-2 for t in range(4)],color=colores_ame[0], s=100)
         axis[0].plot(X, 0.17*X**3-0.7*X**2+2*X-2, color=colores_ame[0], linestyle='dotted')
     
         #Cúbica Legendre
         axis[1].scatter(eje, base_legendre[3],color=colores_ame[0], s=100)
         axis[1].plot(X, 0.745*X**3-3.35*X**2+3.5*X-0.22, color=colores_ame[0], linestyle='dotted')
         #TODO cómo había encontrado esos coeficientes???
     
     
         for i in range(2):
             axis[i].plot([-0.5,3.5], [0,0], color=colores_ame[1])
             axis[i].plot([0,0], [-2.5,3.5], color=colores_ame[1])
             axis[i].grid()
         plt.xlim([-0.1,3.1])
         plt.ylim([-2.2,2.2])
         return plt.show()


def figura_defGrado():

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
    
    def q(x):
        return -0.5*x**2+0.5*x+1 
    
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
    
    axis[1].plot(X,q(X), color=colores_ame[0])
    axis[1].scatter(dominio,senial, color=colores_ame[0], label='$x=(1,1,0) \in \mathbb{R}^{3}$')
    
    
    
    for i in range(2):
        flechas_2D.dibujar_flechas_2d(fig, axis[i])
        axis[i].grid(True)
        axis[i].legend()
    return plt.show()



def figura_cosenoMuestreo(w, phi, n):
    """
    'w' es un real mayor a cero (la frecuencia), 'phi' \\in[0,1] es el desfase normalizado.
    'n', entero mayor a uno, indica la cantidad de muestras a tomar'
    """
    pi = math.pi 

    dominio = [t/n for t in range(n)]
    valores = [np.cos(2 * pi * w * t+ 2*pi*phi) for t in dominio]

    fig, axis = plt.subplots(1,1)
    X1 = np.arange(-0.3,0,0.01)
    X2 = np.arange(0,1,0.01)
    X3 = np.arange(1,1.3,0.01)

    axis.plot(X1, np.cos(2 * pi * w * X1 + 2*pi*phi), color = 'black', linestyle = 'dotted')
    axis.plot(X2, np.cos(2 * pi * w * X2+ 2*pi*phi), color = 'red', linestyle = 'dotted')
    axis.scatter(dominio, valores, color = 'red')
    axis.plot(X3, np.cos(2 * pi * w * X3+ 2*pi*phi), color = 'black', linestyle = 'dotted')
    axis.axhline(y=0, color = 'gray')
    axis.axvline(x=0, color = 'gray')
    axis.grid(True)

    #dominio = [t/n for t in range(n)]
    #valores = [np.sin(2 * pi * w * t+ 2*pi*phi) for t in dominio]
    #axis.plot(X1, np.sin(2 * pi * w * X1 + 2*pi*phi), color = 'black', linestyle = 'dotted')
    #axis.plot(X2, np.sin(2 * pi * w * X2+ 2*pi*phi), color = 'blue', linestyle = 'dotted')
    #axis.scatter(dominio, valores, color = 'blue')
    #axis.plot(X3, np.sin(2 * pi * w * X3+ 2*pi*phi), color = 'black', linestyle = 'dotted')
    
    return plt.show()


def figura_raicesUnidad(n1, n2):
    pi = math.pi
    fig, axis =plt.subplots(1,2)

    dominio_1 = [k/n1 for k in range(n1)]
    dominio_2 = [k/n2 for k in range(n2)]

    cosenos_puntos_1=[math.cos(2*pi*t) for t in dominio_1]
    senos_puntos_1=[math.sin(2*pi*t) for t in dominio_1]

    cosenos_puntos_2=[math.cos(2*pi*t) for t in dominio_2]
    senos_puntos_2=[math.sin(2*pi*t) for t in dominio_2]

    t = np.linspace(0, 2*np.pi, 500)

    for i in range(2):
        formato_axis(axis[i])
        axis[i].plot(np.cos(t), np.sin(t), color='black', linestyle = 'dotted')
    for i in range(n1):
        axis[0].scatter(cosenos_puntos_1[i], senos_puntos_1[i], color = colores[1])
    for i in range(n2):
        axis[1].scatter(cosenos_puntos_2[i], senos_puntos_2[i], color = colores[1])

    axis[0].set_title('$n={{ {0} }}$'.format(str(n1)))
    axis[1].set_title('$n={{ {0} }}$'.format(str(n2)))
    return plt.show()


def coeficientes_interpolacion(vector_abscisas, vector_ordenadas, n):
    """
    TODO actualiza
    dado un vector de ordenadas de dimensión n...
    Se pide n como argumento para no calcularlo.
    """
    A = np.empty((n,n), int)
    for i in range(n):
        entrada = vector_abscisas[i]
        A[i] = [entrada**j for j in range(n)]
    return np.linalg.solve(A, vector_ordenadas)


def PDL_grafica_versionContinua(n,k):
    fig, axis = plt.subplots(1,1)
    dominio = [m for m in range(n)]
    vector_legendre = legendre.calculo_base(n)[k]
    coeficientes =coeficientes_interpolacion(dominio, vector_legendre, n) 
    X = np.arange(0,n-1, 0.0001)
    vector_legendre = legendre.calculo_base(n)[k] 
    axis.scatter(dominio, vector_legendre, color = 'hotpink')
    
    def PDL_continuo(t):
        resultado = coeficientes[0]
        for i in range(1,n):
            resultado += coeficientes[i]*t**i
        return resultado
    
    axis.plot(X, PDL_continuo(X), color = '#09DEB4')
    formato_axis(axis)
    
    return plt.show()

def PDL_grafica_versionContinua_axis(vector_legendre, n, k, dominio, axis, Color):
    coeficientes =coeficientes_interpolacion(dominio, vector_legendre, n) 
    X = np.arange(0,n-1, 0.0001)
    axis.scatter(dominio, vector_legendre, color = Color, label = str(k))
    
    def PDL_continuo(t):
        resultado = coeficientes[0]
        for i in range(1,n):
            resultado += coeficientes[i]*t**i
        return resultado
    
    axis.plot(X, PDL_continuo(X), color = Color, linestyle = 'dotted')

def PDLdim_n_graficas_continuas(n, lista_colores):
    fig, axis = plt.subplots(1,1)
    dominio = [m for m in range(n)]
    base_legendre = legendre.calculo_base(n)
    for k in range(n):
        vector_legendre = base_legendre[k]
        PDL_grafica_versionContinua_axis(vector_legendre,n,k, dominio, axis, lista_colores[k])
    
    formato_axis(axis)
    plt.title('Gráficas de los polinomios de Legendre discretos de dimensión '+str(n))
    return plt.show()

def ejemplificando_operador_alternancia(x, nombre):
    """
    x es un array.
    Se grafican la gráfica de x y la de su alternado.
    """
    n = len(x)
    Ax = []
    for m in range(n):
        Ax.append((-1)**m * x[m])
    
    fig, axis = plt.subplots(1,2)
    X = np.arange(0,1,0.05)
    dominio = [m/n for m in range(n)]
    axis[0].scatter(dominio, x, color = 'hotpink')
    axis[1].scatter(dominio, Ax, color = 'hotpink')
    axis[0].plot(dominio, x, color = 'hotpink', linestyle = 'dotted')
    axis[1].plot(dominio, Ax, color = 'hotpink', linestyle = 'dotted')

    formato_axis(axis[0])
    axis[0].set_title('Gráfica de '+ '${0}$'.format(nombre))
    formato_axis(axis[1])
    axis[1].set_title('Gráfica de '+ r'$A_{{ {0}  }}($'.format(str(n)) + '${0}$'.format(nombre) + ')')

    return plt.show()

if __name__=='__main__':
    #n, k = 15,1
    #x = legendre.calculo_base(n)[k]
    #nombre = r'\mathcal{{L}}^{{{0}}}'.format(str(n) + ', ' + str(k))
    #ejemplificando_operador_alternancia(x, nombre)
    #figura_discretizacionPuntual()
    #figura_cambioDeMalla()
    figura_rectaParabola()
    #figura_introduccion()
    #figura_coseno()
    #figura_demSimetrias()
    #figura_ortogYoscil(2)  
    #figura_defGrado() #problemas
    #figura_cosenoMuestreo(0.0001,0.3, 20)
    figura_cosenoMuestreo(4.6, 0.3, 18)
    #figura_cosenoMuestreo(12, 0, 24)
    #figura_raicesUnidad(5, 8)
    #PDL_grafica_versionContinua(8,5)  
    #figura_generalizacion_operadorOmega()
    #lista_colores =['#fe00fa', '#feaf16', '#782ab6', '#1c8356','#3283fe', '#16ff32', '#dea0fd',  '#325a9b', '#feaf16', '#feaf16', '#1cffce'] 
    #for n in range(2,9):
    #    PDLdim_n_graficas_continuas(n, lista_colores)

