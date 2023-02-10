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
fig, axis= plt.subplots(1,2)
fig.suptitle('Gráficas de los elementos de la base $\mathcal{L}^{2}$')
dominio=np.array([0,1])
X=np.arange(0,1.05, 0.05) #para trazar la gráfica de la parábola
axis[0].scatter(dominio, legendre.BON_L[2][0], color=colores_ame[0], s=100,label='$\mathcal{L}^{2,{0}}$')
axis[1].scatter(dominio, legendre.BON_L[2][1], color=colores_ame[0], s=100,label='$\mathcal{L}^{2,{1}}$')
for i in range(2):
    axis[i].set_xlim(-0.5,1.5)
    axis[i].set_ylim(-1.2,1.2)
    graficas.dibujos_axis(axis[i], dominio)
    flechas_2D.dibujar_flechas_2d(fig, axis[i])
    axis[i].grid(True)


#Graficando las rectas/parábolas de mínimos cuadrados para las primeras tres gráficas
axis[0].plot(dominio, legendre.BON_L[2][0], color=colores_ame[0], linestyle=':')
axis[1].plot(dominio, legendre.BON_L[2][1], color=colores_ame[0], linestyle=':')



#axis[0].plot(dominio, legendre.BON_L[2][0], color='goldenrod', linestyle=':', label='Gráfica de $\mathcal{L}^{2,0}$')
#axis[1].plot(dominio, legendre.BON_L[2][1], color='goldenrod', linestyle=':', label='Gráfica de $\mathcal{L}^{2,1}$')

plt.show()
