import numpy as np
import matplotlib.pyplot as plt
import pylab #Para usar LaTeX en captions
import math
import legendre #En este archivo de python se guardan las BON de Legendre discretas de dimensión desde 2 hasta 6
import flechas_2D
import proyecciones as proy
import graficasLegendre as graficas


colores_ame=['goldenrod','hotpink','mediumpurple']



#Dimensión 6-----------------------------------------------
fig, axis= plt.subplots(2,3)
fig.suptitle('Gráficas de los elementos de la base $\mathcal{L}^{6}$')
dominio=np.array([0,1,2,3,4,5])
X=np.arange(0,5.05, 0.05) #para trazar la gráfica de la parábola
axis[0,0].scatter(dominio, legendre.BON_L[6][0], color=colores_ame[0], s=100,label='$\mathcal{L}^{6,{0}}$')
axis[0,1].scatter(dominio, legendre.BON_L[6][1], color=colores_ame[0], s=100,label='$\mathcal{L}^{6,{1}}$')
axis[0,2].scatter(dominio, legendre.BON_L[6][2], color=colores_ame[0], s=100,label='$\mathcal{L}^{6,{2}}$')
axis[1,0].scatter(dominio, legendre.BON_L[6][3], color=colores_ame[0], s=100,label='$\mathcal{L}^{6,{3}}$')
axis[1,1].scatter(dominio, legendre.BON_L[6][4], color=colores_ame[0], s=100,label='$\mathcal{L}^{6,{4}}$')
axis[1,2].scatter(dominio, legendre.BON_L[6][5], color=colores_ame[0], s=100,label='$\mathcal{L}^{6,{5}}$')
for i in range(2):
    for j in range(3):
        axis[i][j].set_xlim(-0.5,6.5)
        axis[i][j].set_ylim(-1.2,1.2)
        graficas.dibujos_axis(axis[i,j], dominio)
        flechas_2D.dibujar_flechas_2d(fig, axis[i,j])
        axis[i,j].grid(True)


#Graficando las rectas/parábolas de mínimos cuadrados para las primeras tres gráficas
axis[0][0].plot(dominio, legendre.BON_L[6][0], color=colores_ame[0], linestyle=':')
axis[0][1].plot(dominio, legendre.BON_L[6][1], color=colores_ame[0], linestyle=':')
parab=proy.parab3Puntos(legendre.BON_L[6][2])
axis[0][2].plot(X, parab[0]*X**2+parab[1]*X+parab[2], color=colores_ame[0], linestyle=':')




plt.show()
