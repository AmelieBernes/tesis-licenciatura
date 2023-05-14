from mpl_toolkits import mplot3d
import numpy as np
import math
import pylab
import matplotlib.pyplot as plt
import proyecciones as proy



from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))

        return np.min(zs)


#Script con funciones útiles para dibujar gráficas en 3d

colores_ame=['goldenrod', 'hotpink', 'rebeccapurple', 'blueviolet']

def dibuja_ejes_labelsPersonalizados(axis, m, labelx, labely, labelz):
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
    a=Arrow3D([m,m+0.5], [0,0],[0,0], arrowstyle='-|>', lw=2, color='gray', mutation_scale=20)
    b=Arrow3D([0,0],[m, m+0.5], [0,0],  arrowstyle='-|>', lw=2, color='gray', mutation_scale=20)
    c=Arrow3D([0,0], [0,0], [m,m+0.5], arrowstyle='-|>', lw=2, color='gray', mutation_scale=20)
    axis.add_artist(a)
    axis.add_artist(b)
    axis.add_artist(c)
    axis.text(m+0.8,0,0.3, labelx, color='gray')
    axis.text(0,m+0.8,0.3, labely, color='gray')
    axis.text(0,0,m+1, labelz, color='gray')
    xx, yy,=np.meshgrid(range(-5,5),range(-5,5))
    planoXY=0*xx+0*yy
    axis.plot_surface(xx,yy,planoXY, color='gray', alpha=0.1)

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
    a=Arrow3D([m,m+0.5], [0,0],[0,0], arrowstyle='-|>', lw=2, color='gray', mutation_scale=20)
    b=Arrow3D([0,0],[m, m+0.5], [0,0],  arrowstyle='-|>', lw=2, color='gray', mutation_scale=20)
    c=Arrow3D([0,0], [0,0], [m,m+0.5], arrowstyle='-|>', lw=2, color='gray', mutation_scale=20)
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

