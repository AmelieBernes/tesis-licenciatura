import matplotlib.pyplot as plt
import numpy as np
import pylab
import math
import legendre
import proyecciones as proy

#(no muy importante) script para dibujar algunas de las gráficas del Ejemplo 1 de la tesis.
#caso particular para dibujar las BON de legendre; quiero dejarlas porque aquí no dibujo los ejes x y y.

#--DIMENSION 2

fig, axis =plt.subplots(1,2)
fig.suptitle('Gráficas de los elementos de $\mathcal{L}^{2}$')
dominio=np.array([0,1])
X=np.arange(0, 1.05, 0.05)
for i in range(2):
    axis[i].set_xlim(-0.5,1.5)
    axis[i].set_ylim(-1.5,1.5)
    axis[i].grid(True)
    mediciones=legendre.BON_L[2][i]
    axis[i].scatter(dominio, mediciones, s=100, color='hotpink')
    axis[i].set_title('Grado '+str(i))

axis[0].plot(dominio, legendre.BON_L[2][0], color='hotpink', linestyle=':', label='Gráfica de $\mathcal{L}^{2,0}$')
axis[1].plot(dominio, legendre.BON_L[2][1], color='hotpink', linestyle=':', label='Gráfica de $\mathcal{L}^{2,1}$')

#--DIMENSION 3

#fig, axis =plt.subplots(1,3)
#fig.suptitle('Gráficas de los elementos de $\mathcal{L}^{3}$')
#dominio=np.array([0,1,2])
#X=np.arange(0, 2.05, 0.05)
#for i in range(3):
#    axis[i].set_xlim(-0.5,2.5)
#    axis[i].set_ylim(-1.5,1.5)
#    axis[i].grid(True)
#    mediciones=legendre.BON_L[3][i]
#    axis[i].scatter(dominio, mediciones, s=100, color='hotpink')
#    axis[i].set_title('Grado '+str(i))
#
#b=proy.coef_RMC(dominio, legendre.BON_L[3][0])
#y=b[0]+b[1]*dominio
#axis[0].plot(dominio, y, color='hotpink', linestyle=':', label='Gráfica de $\mathcal{L}^{3,0}$')
#
#
#b=proy.coef_RMC(dominio, legendre.BON_L[3][1])
#y=b[0]+b[1]*dominio
#axis[1].plot(dominio, y, color='hotpink', linestyle=':', label='Gráfica de $\mathcal{L}^{3,1}$')
#
#parab=proy.parab3Puntos(legendre.BON_L[3][2])
#axis[2].plot(X, parab[0]*X**2+parab[1]*X+parab[2], color='hotpink', linestyle=':', label='Gráfica de $\mathcal{L}^{3,2}$')


plt.show()



