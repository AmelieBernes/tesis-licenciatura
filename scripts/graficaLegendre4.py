import numpy as np
import matplotlib.pyplot as plt
import pylab #Para usar LaTeX en captions
import math
import legendre #En este archivo de python se guardan las BON de Legendre discretas de dimensión desde 2 hasta 6
import flechas_2D
import proyecciones as proy
import graficasLegendre as graficas


colores_ame=['goldenrod','hotpink','mediumpurple']



#Dimensión 4-----------------------------------------------
fig, axis= plt.subplots(2,2)
fig.suptitle('Gráficas de los elementos de la base $\mathcal{L}^{4}$')
dominio=np.array([0,1,2,3])
X=np.arange(0,3.05, 0.05) #para trazar la gráfica de la parábola
axis[0,0].scatter(dominio, legendre.BON_L[4][0], color=colores_ame[0], s=100,label='$\mathcal{L}^{4,{0}}$')
axis[0,1].scatter(dominio, legendre.BON_L[4][1], color=colores_ame[0], s=100,label='$\mathcal{L}^{4,{1}}$')
axis[1,0].scatter(dominio, legendre.BON_L[4][2], color=colores_ame[0], s=100,label='$\mathcal{L}^{4,{2}}$')
axis[1,1].scatter(dominio, legendre.BON_L[4][3], color=colores_ame[0], s=100,label='$\mathcal{L}^{4,{3}}$')
for i in range(2):
    for j in range(2):
        axis[i][j].set_xlim(-0.5,3.5)
        axis[i][j].set_ylim(-1.2,1.2)
        graficas.dibujos_axis(axis[i,j], dominio)
        flechas_2D.dibujar_flechas_2d(fig, axis[i,j])
        axis[i,j].grid(True)


#Graficando las rectas/parábolas de mínimos cuadrados para las primeras tres gráficas
axis[0][0].plot(dominio, legendre.BON_L[4][0], color=colores_ame[0], linestyle=':')
axis[0][1].plot(dominio, legendre.BON_L[4][1], color=colores_ame[0], linestyle=':')
parab=proy.parab3Puntos(legendre.BON_L[4][2])
axis[1][0].plot(X, parab[0]*X**2+parab[1]*X+parab[2], color=colores_ame[0], linestyle=':')




plt.show()
