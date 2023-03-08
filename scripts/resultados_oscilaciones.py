"""
Script para graficar la información del estudio de oscilaciones
de los polinomios discretos de Legendre via representación en sistemas de oscilaciones.

Input: n \in \IN, 0 \leq k leq n-1, ambos objetos de tipo 'int'. 'n' es la dimensión del polinomio,
'k' su grado.


La figura consta de cuatro axis:
    1.- axis[0][0]: En el que se grafica, en el intervalo [0,1], la gráfica de la señal \mathcal{L}^{n,k}. Se hace en 
    [0,1] en lugar de en [0, n-1] pues es en el primer intervalo en el que se grafican las funciones sinusoidales a partir
    de las cuales se obtienen los vectores de frecuencia.

    2.- axis[0][1] Gráfica de los coeficientes \sigma_{n,k} como se calculan en la función 'calculando_sigmasYesp'
    del script 'oscilaciones_legendre.py'
    TODO: puedes seguir graficando la esperanza si realizas el paso preliminar de normalizar el vector de coeficientes sigma.
    Esto no era necesario cuando usabas la base clásica de oscilaciones de Fourier, pero puesto que nuestras modificaciones
    pueden no ser BON's, ya no estas segura de obtener de inmediato con las sigmas una distribución discreta de probabilidad.

"""

import numpy as np
import matplotlib.pyplot as plt
import pylab #Para usar LaTeX en captions
import math


import legendre_discreto as legendre 
import oscilaciones_legendre as osc
import proyecciones as proy

colores_amelie=['hotpink', 'mediumpurple', 'darkgoldenrod']
n=50
k=9





fig, axis= plt.subplots(2,1)
fig.suptitle(r"Análisis oscilatorio del polinomio discreto $ \mathcal{{ L }}^{{ {0} , {1} }} \in \mathbb{{ R }}^{{ {0}  }}$".format(str(n), str(k)) )

for j in range(1):
    axis[j].axhline(y=0, color='gray')
    axis[j].axvline(x=0, color='gray')
    axis[j].grid(True)



#----------- Axis[0,0]: WIP -------------------------

dominio=[k/n for k in range(n)]
vector_legendre=legendre.base_Legendre(n)[k] #Calculamos el vector de legendre de interés para el análisis.
axis[0].scatter(dominio, vector_legendre, color=colores_amelie[0])



#----------- Axis[0,1]: WIP -------------------------
(sigma, esp)= osc.calculando_sigmasYesp(n, k)
dominio_sigma=[t for t in range(len(sigma))]
sigma_max=max(sigma)
print(sigma_max)
#TODO falta realizar normalizaciones en el vector de sigmas para que calcular la esperanza tenga sentido siempre.

axis[1].scatter(dominio_sigma, sigma, s=100, color=colores_amelie[1], marker="*")
axis[1].scatter(esp, 0, s=100, color=colores_amelie[2], marker="^", label='Esperanza: '+str(esp.round(4))) 




plt.grid()
plt.show()