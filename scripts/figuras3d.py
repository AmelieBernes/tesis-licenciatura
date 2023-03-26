from mpl_toolkits import mplot3d
import numpy as np
import math
import pylab
import matplotlib.pyplot as plt
import ame_flecha
import flechas_2D
import proyecciones as proy

#Script con funciones útiles para dibujar gráficas en 3d

colores_ame=['goldenrod', 'hotpink', 'rebeccapurple', 'blueviolet']

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


