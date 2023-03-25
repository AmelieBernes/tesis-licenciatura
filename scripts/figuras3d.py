from mpl_toolkits import mplot3d
import numpy as np
import math
import pylab
import matplotlib.pyplot as plt
import ame_flecha
import flechas_2D
import proyecciones as proy

#Recuerda que si quieres llamar a este script en otro para usar funciones aquí definidas NO es bueno que este scritp, en sí mismo, haga algo. 

colores_ame=['goldenrod', 'hotpink', 'rebeccapurple', 'blueviolet']

#estoy trabajando este script de manera que sólo tenga funciones útiles a la hora de graficar cosas en tercera dimensión.

def dibuja_ejes(axis, m):
    """
    Función que dibuja los ejes x, y, z del espacio tridimensional, con flechas para las
    ramas positivas y etiquetas. 
    m es un real positivo que da el máximo valor positivo para los tres ejes.
    """
    axis.set_xlim(m, -m)
    axis.set_ylim(m, -m)
    axis.set_zlim(m, -m)
    X=np.arange(-m,m, 0.05)
    Y=np.arange(-m,m, 0.05)
    Z=np.arange(-m,m, 0.05)
    axis.plot(X, 0*X, 0*X, color='gray')
    axis.plot(0*X,X, 0*X, color='gray')
    axis.plot(0*X, 0*X,X, lw=2, color='gray')
    #Dibujando las flechas de los ejes
    a=ame_flecha.Arrow3D([m,m+0.5], [0,0],[0,0], arrowstyle='-|>', lw=2, color='gray', mutation_scale=20)
    b=ame_flecha.Arrow3D([0,0],[m, m+0.5], [0,0],  arrowstyle='-|>', lw=2, color='gray', mutation_scale=20)
    c=ame_flecha.Arrow3D([0,0], [0,0], [m,m+0.5], arrowstyle='-|>', lw=2, color='gray', mutation_scale=20)
    axis.add_artist(a)
    axis.add_artist(b)
    axis.add_artist(c)
    axis.text(m+0.8,0,0.3, 'x', color='gray')
    axis.text(0,m+0.8,0.3, 'y', color='gray')
    axis.text(0,0,m+1, 'z', color='gray')
    xx, yy,=np.meshgrid(range(-5,5),range(-5,5))
    planoXY=0*xx+0*yy
    axis.plot_surface(xx,yy,planoXY, color='gray', alpha=0.1)

def dibuja_plano(axis,a,b,Color='mediumpurple',n=3):
    """
    Función que grafica el plano de ecuación cartesiana z=a*x+y.
    """
    xx, yy=np.meshgrid(range(-n,n),range(-n,n))
    coord_z=a*xx+b*yy
    axis.plot_surface(xx,yy,coord_z, color=Color, alpha=0.6)
    axis.plot_wireframe(xx,yy,coord_z, color='white', alpha=0.4)
#creo que las siguientes funciones, por ser tan específicas, son inútiles.


def dibuja_W2_de_R3(n, Color='mediumpurple'):
    xx, yy=np.meshgrid(range(-n,n),range(-n,n))
    espacio_W2=2*yy-xx
    axis.plot_surface(xx,yy, espacio_W2, color=Color, alpha=0.6)
    axis.plot_wireframe(xx,yy, espacio_W2, color='white', alpha=0.4)

def dibuja_W1_de_R3(n, Color='hotpink'):
    X=np.arange(-n,n,0.5)
    axis.plot(X,X,X, color='hotpink')


#Después deberías quitar estos comentarios.

#def dibuja_espaciosLegendre_R3(m, n):
#    global fig
#    fig=plt.figure()
#    global axis
#    axis=plt.axes(projection='3d')
#    dibuja_W2_de_R3(n)
#    dibuja_ejes(m)
#    
#    axis.invert_zaxis() #para que la orientación del espacio sea la correcta!! No sé por qué esto no es así por default
#    plt.show()


#----------------------------------------------------------------------------------------

#figuras en 3d para algunos ejemplos de tal sección.


fig=plt.figure()
axis=plt.axes(projection='3d')
dibuja_W2_de_R3(3)
dibuja_ejes(5)
axis.quiver(0,0,0,1,-2,1,color=colores_ame[0]) #vector que apuntan en la dirección L_32
axis.scatter(1,-2,1,color=colores_ame[0])

axis.quiver(0,0,0,1,-2,4,color=colores_ame[1]) #vector x
axis.scatter(1,-2,4, color=colores_ame[1])
axis.quiver(0,0,0,-0.5,1,2.5,color=colores_ame[3]) #proyección al espacio W_{1}
axis.scatter(-0.5,1,2.5,color=colores_ame[1])
axis.plot([-0.5,1],[1,-2],[2.5,4], linestyle='--',color=colores_ame[1])

#segunda sección
axis.quiver(0,0,0,-1,1,3,color=colores_ame[1]) #vector x
axis.scatter(-1,1,3, color=colores_ame[1])

#tercera sección
axis.quiver(0,0,0,-3,3,-3,color=colores_ame[1]) #vector x
axis.scatter(-3,3,-3,color=colores_ame[1])
axis.quiver(0,0,0,-1,-1,-1,color=colores_ame[3]) #proyección al espacio W_{1}
axis.scatter(-1,-1,-1,color=colores_ame[1])
axis.plot([-3,-1],[3,-1],[-3,-1], linestyle='--',color=colores_ame[1])


