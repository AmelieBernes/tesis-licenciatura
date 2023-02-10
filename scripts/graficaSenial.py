import numpy as np
import matplotlib.pyplot as plt
import pylab #Para usar LaTeX en captions
import math
import legendre #En este archivo de python se guardan las BON de Legendre discretas de dimensión desde 2 hasta 6
import flechas_2D
import proyecciones as proy
import graficasLegendre as graficas


"""
Script para dibujar la gráfica de una señal.
"""

#el input: un array de mediciones, de longitud cualquiera.
mediciones=[-0.5, 2.4, 1.6, 1.7, 2.3]
n=len(mediciones)
dominio=[]#inicializamos el dominio..
#...y lo llenamos
for i in range(n):
    dominio.append(i)


plt.scatter(dominio, mediciones, s=150, color='goldenrod')
plt.grid()
plt.show()
