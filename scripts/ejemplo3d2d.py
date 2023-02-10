import numpy as np
from mpl_toolkits import mplot3d
import math
import pylab
import matplotlib.pyplot as plt
import ame_flecha
import flechas_2D
import legendre
import proyecciones as proy
import figuras3d

colores_ame=['goldenrod', 'hotpink', 'rebeccapurple', 'blueviolet', 'lawngreen']
pi=3.141592


#input ------------

##Lista con todos los ángulos que voy a considerar.
parte_afin=math.sqrt(3)*legendre.BON_L[3][0]+math.sqrt(2)*legendre.BON_L[3][1] #es [0,1,2]
alpha=[pi/3, pi/4, 0, pi/3, 6*pi/17]
signo=[1,1,0,-1,-1] #lista de signos (para determinar la región.)

a2=[signo[i]*math.sqrt(5*(math.tan(alpha[i]))**2) for i in range(len(alpha))] #lista de coeficientes a_2, que dependen de los ángulos.
vectores_x=[parte_afin+a2[i]*legendre.BON_L[3][2] for i in range(len(alpha))]
#fin del input -----------

fig = plt.figure()


# --------------------- Plot en 3D -----------------------------
ax = fig.add_subplot(1, 2, 1, projection='3d')
#dibuja espacio W2 de R3, con parámetro 3
figuras3d.dibuja_plano(ax,-1,2)

#Dibuja ejes, con parámetro 5
figuras3d.dibuja_ejes(ax, 5)

ax.quiver(0,0,0,0,1,2,color='black') #vector parte afin común. 
ax.scatter(0,1,2, color='black', label='Vector $(0,1,2)$')
X=np.arange(-3,3,0.05)
ax.plot(X, 1-2*X, 2+X, color='gray', linestyle=':', label='Recta $l_{0,1,2}$') #recta
ax.invert_zaxis() #para que la orientación del espacio sea la correcta!


def display_3d(i):
    """
    i es un entero entre cero y len(alpha).
    """
    ax.quiver(0,0,0,vectores_x[i][0],vectores_x[i][1],vectores_x[i][2],color=colores_ame[1])
    ax.scatter(vectores_x[i][0],vectores_x[i][1],vectores_x[i][2],color=colores_ame[1], s=100, label='$Vector x$')
    ax.legend()
    plt.legend()


# Plot en 2D
axis = fig.add_subplot(1, 2, 2)
axis.set_xlim(-0.2,2.2)
axis.set_ylim(-8.2,8.2)
flechas_2D.dibujar_flechas_2d(fig, axis)


def display_graphSignal(Axis, vect_x):
    X=np.arange(0,2.2,0.05)
    dominio=[0,1,2]
    plt.scatter(dominio, vect_x, s=130, color=colores_ame[1], label='Gráfica del vector $x \in \mathbb{R}^{3}$')
    proyW2=proy.proyeccion(vect_x,2)
    parab=proy.parab3Puntos(proyW2)
    Axis.plot(X, parab[0]*X**2+parab[1]*X+parab[2], linestyle='dotted', color=colores_ame[1])
    plt.legend()

#proy.graficas_senial_parteCuadrV2(x)
axis.grid(True)



#output 
def output_ame(i):
    """
    i es un entero entre 0 y len(alpha)
    """
    display_graphSignal(axis, vectores_x[i])
    display_3d(i)
    plt.show()

output_ame(4)
