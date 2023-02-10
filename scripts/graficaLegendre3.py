import numpy as np
import matplotlib.pyplot as plt
import pylab #Para usar LaTeX en captions
import math
import legendre #En este archivo de python se guardan las BON de Legendre discretas de dimensión desde 2 hasta 6
import flechas_2D
import proyecciones as proy
import graficasLegendre as graficas


colores_ame=['goldenrod','hotpink','mediumpurple']



#Dimensión 3-----------------------------------------------
fig, axis= plt.subplots(1,3)
fig.suptitle('Gráficas de los elementos de la base $\mathcal{L}^{3}$')
dominio=np.array([0,1,2])
X=np.arange(0,2.05, 0.05) #para trazar la gráfica de la parábola
axis[0].scatter(dominio, legendre.BON_L[3][0], color=colores_ame[0], s=100,label='$\mathcal{L}^{4,{0}}$')
axis[1].scatter(dominio, legendre.BON_L[3][1], color=colores_ame[0], s=100,label='$\mathcal{L}^{4,{2}}$')
axis[2].scatter(dominio, legendre.BON_L[3][2], color=colores_ame[0], s=100,label='$\mathcal{L}^{4,{3}}$')
for i in range(3):
    axis[i].set_xlim(-0.5,2.5)
    axis[i].set_ylim(-1.2,1.2)
    graficas.dibujos_axis(axis[i], dominio)
    flechas_2D.dibujar_flechas_2d(fig, axis[i])
    axis[i].grid(True)


#Graficando las rectas/parábolas de mínimos cuadrados para las primeras tres gráficas
axis[0].plot(dominio, legendre.BON_L[3][0], color=colores_ame[0], linestyle=':')
axis[1].plot(dominio, legendre.BON_L[3][1], color=colores_ame[0], linestyle=':')
parab=proy.parab3Puntos(legendre.BON_L[3][2])
axis[2].plot(X, parab[0]*X**2+parab[1]*X+parab[2], color=colores_ame[0], linestyle=':')

plt.show()
