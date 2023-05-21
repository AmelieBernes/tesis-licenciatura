import numpy as np
from numpy.linalg import norm
from mpl_toolkits import mplot3d
import math
import pylab
import matplotlib.pyplot as plt

import proyecciones as proy
import flechas_2D
import base_legendreDiscreta as legendre
import funciones_figuras3d

#Script en el que se dan los códigos para las figuras del ejemplo
#en 3D de la tesis.

"""
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Computer Modern Serif'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Computer Modern Serif:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Computer Modern Serif:bold'
"""

colores=['goldenrod', 'hotpink', 'rebeccapurple', 'blueviolet', 'lawngreen', 'darkturquoise', 'mediumpurple']
pi=3.141592

def dibujar_W31_de_R3(axis, n, Color=colores[6]):
    xx, yy=np.meshgrid(range(-n,n),range(-n,n))
    espacio_W2=2*yy-xx

    #https://stackoverflow.com/questions/55531760/is-there-a-way-to-label-multiple-3d-surfaces-in-matplotlib
    c1= axis.plot_surface(xx,yy, espacio_W2, color=Color, alpha=0.6, label='$W_{3,1}= \{ (u,v, 2v-u) \in \mathbb{R}^{3} :u, v, \in \mathbb{R}  \} $')
    c1._facecolors2d=c1._facecolor3d
    c1._edgecolors2d=c1._edgecolor3d
    axis.plot_wireframe(xx,yy, espacio_W2, color='white', alpha=0.4)

def dibujar_W30_de_R3(axis, n, Color=colores[5]):
    X=np.arange(-n,n,0.5)
    axis.plot(X,X,X, color=Color, label= '$W_{3,0}=\{ (u, u, u) \in \mathbb{R}^{3} : u \in \mathbb{R} \}$')


#TODO nuevo:
def dibujar_plano_por_tres_puntos(p0, p1, p2, axis, titulo, Color):
    """
    p0, p1 y p2 son los puntos por los que pasa el plano.
    axis es el axis de una figura en el que vamos a dibujar tal plano.
    'titulo' es el string del título.
    'Color' es una string con el color del plano.
    """
    u, v = [], [] #puntos del plano
    for i in range(3):
        u.append(p1[i]-p0[i]) 
        v.append(p2[i]-p0[i]) 
    n = np.cross(u, v) # un vector normal al plano.
    
    #TODO incluir luego el caso n[2]==0.
    coef_x = -n[0]/n[2]
    coef_y = -n[1]/n[2]
    coef_cte = p0[2] + p0[0]*n[0]/n[2] + p0[1]*n[1]/n[2]
    
    xx, yy = np.meshgrid(range(-10, 10), range(-10, 10))
    zz = coef_x*xx - coef_y *yy + coef_cte
    axis.plot_surface(xx, yy, zz, color = Color, alpha = 0.5)
    axis.set_title(titulo)


def dibujar_espaciosLegendre_R2_R3(m, n):
    fig = plt.figure()
    
    # --- Graficando en R2 a W_{2,1}
    axis = fig.add_subplot(1,2,1) 
    X=np.arange(-1.0,1.0,0.05)
    axis.plot(X, X, color=colores[5], label='$W_{2,0}= \{(u,u): u \in \mathbb{R}  \}$') #TODO después cambia este color
    flechas_2D.dibujar_flechas_2d(fig, axis)
    axis.set_title("$\mathbb{R}^{2}$")
    plt.legend()
    plt.grid()
    
    # --- Graficando en R3 a W_{3,1} y a W_{3,2}
    axis = fig.add_subplot(1,2,2, projection='3d')

    dibujar_W30_de_R3(axis, n)
    dibujar_W31_de_R3(axis, n)
    funciones_figuras3d.dibuja_ejes(axis,m-1)
    axis.invert_zaxis() #para que la orientación del espacio sea la correcta!! No sé por qué esto no es así por default
    axis.set_title("$\mathbb{R}^{3}$")
    plt.legend()
    fig.suptitle('Subespacios de Legendre de dimensiones 2 y 3')
    plt.show()




#----------------------------------------------------------------------------------------

#Funciones para graficar las regiones I, II y III de R3 en las que el hiperplano
#W3,2 \subseteq \IR3 divide a este espacio.
def tres_regiones():
    fig=plt.figure()
    axis = fig.add_subplot(1,3,1,projection='3d')
    dibujar_W31_de_R3(axis, 4)
    funciones_figuras3d.dibuja_ejes(axis, 6)
    axis.quiver(0,0,0,1,-2,1,color=colores[0]) #vector que apuntan en la dirección L_32
    axis.scatter(1,-2,1,color=colores[0])
    
    axis.quiver(0,0,0,1,-2,4,color=colores[1]) #vector x
    axis.scatter(1,-2,4, color=colores[1])
    axis.quiver(0,0,0,-0.5,1,2.5,color=colores[3]) #proyección al espacio W_{3,1}
    axis.scatter(-0.5,1,2.5,color=colores[3])
    axis.plot([-0.5,1],[1,-2],[2.5,4], linestyle='--',color=colores[1])
    axis.set_title("Región I")
    axis.invert_zaxis()

    axis = fig.add_subplot(1,3,2,projection='3d')
    dibujar_W31_de_R3(axis, 5)
    funciones_figuras3d.dibuja_ejes(axis, 6)
    axis.quiver(0,0,0,1,-2,1,color=colores[0]) #vector que apuntan en la dirección L_32
    axis.scatter(1,-2,1,color=colores[0])
    axis.quiver(0,0,0,-1,1,3,color=colores[1]) #vector x
    axis.scatter(-1,1,3, color=colores[1])
    axis.set_title("Región II")
    axis.invert_zaxis()


    axis = fig.add_subplot(1,3,3,projection='3d')

    dibujar_W31_de_R3(axis, 5)
    funciones_figuras3d.dibuja_ejes(axis, 6)
    axis.quiver(0,0,0,1,-2,1,color=colores[0]) #vector que apuntan en la dirección L_32
    axis.scatter(1,-2,1,color=colores[0])
    axis.quiver(0,0,0,-3,3,-3,color=colores[1]) #vector x
    axis.scatter(-3,3,-3,color=colores[1])
    axis.quiver(0,0,0,-1,-1,-1,color=colores[3]) #proyección al espacio W_{1}
    axis.scatter(-1,-1,-1,color=colores[1])
    axis.plot([-3,-1],[3,-1],[-3,-1], linestyle='--',color=colores[1])
    axis.set_title("Región III")
    axis.invert_zaxis()

    fig.suptitle("Las tres regiones en las que el hiperplano $W_{3,1}$ divide a $\mathbb{R}^{3}$")

    return plt.show()


#-----------------------------------------------------------------------------------------

def generador_senAfinesR3(a0, a1, legendre_3):
    """
    'a0' y 'a1' son floats, son los coeficientes del vector afin
    (el output) respecto a la BON de Legendre discreta $\mathcal{L}^{3}$
    """
    return a0*legendre_3[0] + a1*legendre_3[1] 


def grafica_vector_3D_y_2D(parte_afin, alpha, signo):

    """
    'parte afin' es un np.array de dimensión 3, que representa una señal afín de dimensión 3 (i.e. un vector de R3 que pertenece al espacio W_{3,2}). Se interpreta como la proyección del vector 'x' que se planea estudiar.

    'alpha' es un float entre 0 y pi/2 que representa el ángulo que forma 'x' con el plano W_{3,2}

    'signo' es 1, 0 o -1, y se usa para determinar la región en la que yace 'x'.
    """
    c= norm(parte_afin)**2
    a2 = signo*math.sqrt(c*(math.tan(alpha))**2) #tercer coeficiente de 'x' respecto a la BON de legendre discreta
    vector_x = parte_afin+ a2 * legendre_3[2] 
    
    fig = plt.figure()
    
    # --------------------- Plot en 3D -----------------------------
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    dibujar_W31_de_R3(ax, int(math.sqrt(c))+1)
    #funciones_figuras3d.dibuja_plano(ax,-1,2)
    funciones_figuras3d.dibuja_ejes(ax, math.sqrt(c)+1)

    ax.quiver(0,0,0,1,-2,1,color=colores[0]) #vector que apuntan en la dirección L_32
    ax.scatter(1,-2,1,color=colores[0])
    
    ax.quiver(0,0,0,float(parte_afin[0]) ,float (parte_afin[1]), float(parte_afin[2]),color='black') #vector parte afin  
    ax.scatter(parte_afin[0],parte_afin[1],parte_afin[2], color='black', label='$\Pi_{W_{3,1}}(x)=$'+str(np.around(parte_afin, 2)))

    X=np.arange(-3,3,0.05)
    ax.plot(parte_afin[0]+X, parte_afin[1]-2*X, parte_afin[2]+X, color='gray', linestyle=':') #recta
    ax.invert_zaxis() #para que la orientación del espacio sea la correcta!
    
    ax.quiver(0,0,0,vector_x[0], vector_x[1], vector_x[2], color=colores[1])
    ax.scatter(vector_x[0], vector_x[1],vector_x[2], color=colores[1], s=100, label='Vector $x=$'+str(np.around(vector_x, 2)))
    ax.legend()
    plt.legend()
    
    
    # --------------------- Plot en 2D ----------------------------
    axis = fig.add_subplot(1, 2, 2)
    axis.set_xlim(-0.2,2.2)
    axis.set_ylim(-8.2,8.2)
    flechas_2D.dibujar_flechas_2d(fig, axis)
    
    X=np.arange(0,2.2,0.05)
    dominio=[0,1,2]
    plt.scatter(dominio, vector_x, s=130, color=colores[1], label='Gráfica del vector $x \in \mathbb{R}^{3}$')
    proyW2=proy.proyeccion(vector_x,2)
    parab=proy.parab3Puntos(proyW2)
    axis.plot(X, parab[0]*X**2+parab[1]*X+parab[2], linestyle='dotted', color=colores[1])
    plt.legend()
    
    #proy.graficas_senial_parteCuadrV2(x)
    axis.grid(True)

    return plt.show()


if __name__=='__main__':
    p0 = [0,0,0]
    p1= [0, 1,2]
    p2= [-7,5,3]
    fig = plt.figure()
    axis = fig.add_subplot(1,1,1, projection='3d')
    titulo = 'probando'
    dibujar_plano_por_tres_puntos(p0, p1, p2, axis, titulo, 'magenta')
    plt.show()

    #legendre_3=legendre.calculo_base(3)
    #legendre_3=[np.asarray(vector) for vector in legendre_3] #convertimos a numpy para poder hacer operaciones como suma o multiplicación por escalar
    #
    #parte_afin= math.sqrt(3)*legendre_3[0] +math.sqrt(2)*legendre_3[1] #es [0,1,2]  
    #print(parte_afin)
    #print(parte_afin[0])
    #alphas=[pi/3, pi/4, 0, pi/3, 6*pi/17]
    #signos=[1,1,0,-1,-1] #lista de signos (para determinar la región.)
    #
    #for i in range(5):
    #s    grafica_vector_3D_y_2D(parte_afin, alphas[i], signos[i])


    #tres_regiones()
    #dibujar_espaciosLegendre_R2_R3(7,5)
    #parte_afin=generador_senAfinesR3(math.sqrt(3), math.sqrt(2), legendre_3) #es (0,1,2)
    #parte_afin=generador_senAfinesR3(-4.5, 0.7, legendre_3)
    #grafica_vector_3D_y_2D(parte_afin, 0.9, 1)
