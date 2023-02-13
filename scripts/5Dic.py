#Gráficas rosas para el final de la discusión de las formas alternativas para la definición de la BON discreta de Legendre.


import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import flechas_2D
#import formato

#-----------------------------------

#Intento de cambiar el tipo de fuente

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


#-----------------------------------------------------

#Primera Figura (discretización integral)
colores_ame=['hotpink', 'rebeccapurple']
fig, axis=plt.subplots(2,2)
X=np.arange(-1.2, 1.2, 0.05)
x=[-1,-0.5,0,0.5,1]
for i in range(2):
    for j in range(2):
        axis[i][j].set_xlim(-1.3,1.3)
        axis[i][j].set_ylim(-1.3,1.3)
        flechas_2D.dibujar_flechas_2d(fig, axis[i,j])
        axis[i,j].grid(True)
        axis[i,j].plot(X, X**(2*i+j), color=colores_ame[0])
        y=[t**(2*i+j) for t in x]
        axis[i,j].scatter(x, y, color=colores_ame[0])
        k=2*i+j
        axis[i,j].set_title("Grado "+str(k)+", $u_{}$".format(str(k))+"$\mapsto \mathcal{L}^{4,grado}.$".replace("grado", str(k)), fontname="Computer Modern Serif", size=16)


#-----------------------------------------------------

#Segunda Figura (discretización puntual)

#X=np.arange(-1.2, 2.2, 0.05)
#def f1(x):
#    """
#    x un objeto de tipo float
#    """
#    return 0*x+0.3
#def f2(x):
#    """
#    x un objeto de tipo float
#    """
#    return -0.4*x+1
#def f3(x):
#    """
#    x un objeto de tipo float
#    """
#    return 0.5*x**2-0.5*x
#def f4(x):
#    """
#    x un objeto de tipo float
#    """
#    return -x**3+2*x**2+x-2
#funciones=[f1, f2,f3,f4]
#x=[-1,0,1,2]
#signos=["+","-","+","-"]
#for i in range(2):
#    for j in range(2):
#        axis[i][j].set_xlim(-1.3,2.3)
#        axis[i][j].set_ylim(-2.3,2.3)
#        flechas_2D.dibujar_flechas_2d(fig, axis[i,j])
#        axis[i,j].grid(True)
#        axis[i,j].plot(X, funciones[2*i+j](X), color=colores_ame[0])
#        y=[funciones[2*i+j](t) for t in x]
#        axis[i,j].scatter(x, y, color=colores_ame[0], s=100)
#        k=2*i+j
#        axis[i,j].set_title("Grado "+str(k)+", $u_{}$".format(str(k))+"$\mapsto$"+signos[k]+"$\mathcal{L}^{4,grado}.$".replace("grado", str(k)), fontname="Computer Modern Serif", size=16)

plt.show()
