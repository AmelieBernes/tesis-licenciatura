import numpy as np
from mpl_toolkits import mplot3d
import math
import pylab
import matplotlib.pyplot as plt
import ame_flecha
import flechas_2D
import base_legendreDiscreta as legendre
import proyecciones as proy
import figuras3d

#Script en el que se dan los códigos para las figuras del ejemplo
#en 3D de la tesis.

colores_ame=['goldenrod', 'hotpink', 'rebeccapurple', 'blueviolet', 'lawngreen']
pi=3.141592

legendre_3=legendre.calculo_base(3)
legendre_3=[np.asarray(vector) for vector in legendre_3] #convertimos a numpy para poder hacer operaciones como suma o multiplicación por escalar


def dibuja_W2_de_R3(axis, n, Color='mediumpurple'):
    xx, yy=np.meshgrid(range(-n,n),range(-n,n))
    espacio_W2=2*yy-xx
    axis.plot_surface(xx,yy, espacio_W2, color=Color, alpha=0.6)
    axis.plot_wireframe(xx,yy, espacio_W2, color='white', alpha=0.4)

def dibuja_W1_de_R3(axis, n, Color='hotpink'):
    X=np.arange(-n,n,0.5)
    axis.plot(X,X,X, color='hotpink')



def dibuja_espaciosLegendre_R3(m, n):
    global fig
    fig=plt.figure()
    global axis
    axis=plt.axes(projection='3d')
    dibuja_W1_de_R3(axis, n)
    dibuja_W2_de_R3(axis, n)
    dibuja_ejes(axis,m)
    
    axis.invert_zaxis() #para que la orientación del espacio sea la correcta!! No sé por qué esto no es así por default
    plt.show()

#dibuja_espaciosLegendre_R3(7,5)

#----------------------------------------------------------------------------------------

#Funciones para graficar las regiones I, II y III de R3 en las que el hiperplano
#W3,2 \subseteq \IR3 divide a este espacio.
def regionI():
    fig=plt.figure()
    axis=plt.axes(projection='3d')
    dibuja_W2_de_R3(axis, 3)
    dibuja_ejes(axis, 5)
    axis.quiver(0,0,0,1,-2,1,color=colores_ame[0]) #vector que apuntan en la dirección L_32
    axis.scatter(1,-2,1,color=colores_ame[0])
    
    axis.quiver(0,0,0,1,-2,4,color=colores_ame[1]) #vector x
    axis.scatter(1,-2,4, color=colores_ame[1])
    axis.quiver(0,0,0,-0.5,1,2.5,color=colores_ame[3]) #proyección al espacio W_{1}
    axis.scatter(-0.5,1,2.5,color=colores_ame[1])
    axis.plot([-0.5,1],[1,-2],[2.5,4], linestyle='--',color=colores_ame[1])
    axis.invert_zaxis()

    return plt.show()

#segunda sección
def regionII():
    fig=plt.figure()
    axis=plt.axes(projection='3d')
    dibuja_W2_de_R3(axis, 3)
    dibuja_ejes(axis, 5)
    axis.quiver(0,0,0,1,-2,1,color=colores_ame[0]) #vector que apuntan en la dirección L_32
    axis.scatter(1,-2,1,color=colores_ame[0])
    axis.quiver(0,0,0,-1,1,3,color=colores_ame[1]) #vector x
    axis.scatter(-1,1,3, color=colores_ame[1])
    axis.invert_zaxis()

    return plt.show()

#tercera sección
def regionIII():
    fig=plt.figure()
    axis=plt.axes(projection='3d')
    dibuja_W2_de_R3(axis, 3)
    dibuja_ejes(axis, 5)
    axis.quiver(0,0,0,1,-2,1,color=colores_ame[0]) #vector que apuntan en la dirección L_32
    axis.scatter(1,-2,1,color=colores_ame[0])
    axis.quiver(0,0,0,-3,3,-3,color=colores_ame[1]) #vector x
    axis.scatter(-3,3,-3,color=colores_ame[1])
    axis.quiver(0,0,0,-1,-1,-1,color=colores_ame[3]) #proyección al espacio W_{1}
    axis.scatter(-1,-1,-1,color=colores_ame[1])
    axis.plot([-3,-1],[3,-1],[-3,-1], linestyle='--',color=colores_ame[1])
    axis.invert_zaxis()

    return plt.show()

#regionIII()

#-----------------------------------------------------------------------------------------


##Lista con todos los ángulos que voy a considerar.

parte_afin= math.sqrt(3)*legendre_3[0] +math.sqrt(2)*legendre_3[1] #es [0,1,2]  

alphas=[pi/3, pi/4, 0, pi/3, 6*pi/17]
signos=[1,1,0,-1,-1] #lista de signos (para determinar la región.)

#TODO parece qie hay un error en la toeria; calculas mal a2 a partir de a1 y  a0.
def grafica_vector_3D_y_2D(parte_afin, alpha, signo):

    """
    'parte afin' es un np.array de dimensión 3, que representa una señal afín de dimensión 3 (i.e. un vector de R3 que pertenece al espacio W_{3,2}). Se interpreta como la proyección del vector 'x' que se planea estudiar.

    'alpha' es un float entre 0 y pi/2 que representa el ángulo que forma 'x' con el plano W_{3,2}

    'signo' es 1, 0 o -1, y se usa para determinar la región en la que yace 'x'.
    """

    a2 = signo*math.sqrt(5*(math.tan(alpha))**2) #tercer coeficiente de 'x' respecto a la BON de legendre discreta
    vector_x = parte_afin+ a2 * legendre_3[2] 
    
    fig = plt.figure()
    
    # --------------------- Plot en 3D -----------------------------
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    #dibuja espacio W2 de R3, con parámetro 3
    figuras3d.dibuja_plano(ax,-1,2)
    
    #Dibuja ejes, con parámetro 5
    figuras3d.dibuja_ejes(ax, 5)
    
    ax.quiver(0,0,0,0,1,2,color='black') #vector parte afin  
    ax.scatter(0,1,2, color='black', label='Vector $(0,1,2)$')
    X=np.arange(-3,3,0.05)
    ax.plot(X, 1-2*X, 2+X, color='gray', linestyle=':', label='Recta $l_{0,1,2}$') #recta
    ax.invert_zaxis() #para que la orientación del espacio sea la correcta!
    
    ax.quiver(0,0,0,vector_x[0], vector_x[1], vector_x[2], color=colores_ame[1])
    ax.scatter(vector_x[0], vector_x[1],vector_x[2], color=colores_ame[1], s=100, label='Vector $x$')
    ax.legend()
    plt.legend()
    
    
    # --------------------- Plot en 2D ----------------------------
    axis = fig.add_subplot(1, 2, 2)
    axis.set_xlim(-0.2,2.2)
    axis.set_ylim(-8.2,8.2)
    flechas_2D.dibujar_flechas_2d(fig, axis)
    
    X=np.arange(0,2.2,0.05)
    dominio=[0,1,2]
    plt.scatter(dominio, vector_x, s=130, color=colores_ame[1], label='Gráfica del vector $x \in \mathbb{R}^{3}$')
    proyW2=proy.proyeccion(vector_x,2)
    parab=proy.parab3Puntos(proyW2)
    axis.plot(X, parab[0]*X**2+parab[1]*X+parab[2], linestyle='dotted', color=colores_ame[1])
    plt.legend()
    
    #proy.graficas_senial_parteCuadrV2(x)
    axis.grid(True)

    return plt.show()


grafica_vector_3D_y_2D(parte_afin, 0.9, signos[3])
