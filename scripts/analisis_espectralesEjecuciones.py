"""
#  ---------------------------------------- -- ----------------------------------------

Script en el que ejecutamos las siguientes funciones definidas en el módulo
'analisis_espectrales' que importamos:

	1.- grafica_analisisEspectrales(x, nombre, graficar = True), es una función
	    para graficar o guardar el análisis espectral de una señal (array) x
	    
	2.- grafica_analisisEspectrales_PDL(n,k, graficar = True), es una función
	    para graficar el análisis espectral del polinomio discreto de Legendre (PDL)
	    de dimensión n y grado k
	    
	3.- grafica_analisisGlobal_k_fijo(k, graficar = True), es una función para dibujar
	    las gráficas de los puntos de la forma (n, FP0(\mathcal{L}^{n,k}))$
    	    y $(n, FP1(\mathcal{L}^{n,k}))$, con $k < n \leq 69$
    	
    	4.- grafica_3d_n_k_FP(N)
    	
    	5.- grafica_analisisGlobal_n_fijo(n, graficar = True), es una función para dibujar
	    las gráficas de los puntos de la forma $(k, FP0(\mathcal{L}^{n,k}))$
	    y $(k, FP1(\mathcal{L}^{n,k}))$, con $0 \leq k \leq n-1$
	    
	6.- grafica_nube_b0m0_b1m1(), para visualizar las gráficas de nube de los coeficientes
	    $(b_0n, m_0n)$ y $(b_1n, m_1n)$
				

    7.- grafica_pendientes_oOrigen_RMC(), para graficar los puntos
        (n, m0_n), (n, m1_n), (n, b0_n), (n, b1_n), con 3 \leq n \leq 69.
  
    8.- grafica_coefEspectrales_de_frecPrincipales(k), para graficar los puntos
        de la forma (n, tau_{n}(FP0(L^{n,k}))) y (n, sigma_{n}(FP1(L^{n,k}))) para toda
        k+1 \leq n \leq 69
#  ---------------------------------------- -- ----------------------------------------
"""



import numpy as np
from numpy.linalg import norm
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import pylab
from tqdm import tqdm
import random #para ejemplos

import pickle   #for python data persistency
import pandas as pd

#módulos personales
import base_legendreDiscreta as legendre
import proyecciones as proy #aquí tengo una función para hacer regresiones lineales
import funciones_figuras3d
import analisis_espectrales as ae #aquí tengo todas las funciones para realizar análisis espectrales


pi=math.pi
mpl.rcParams.update(mpl.rcParamsDefault)

"""
#  ---------------------------------------- -- ----------------------------------------

				  FUNCIONES PARA CREAR EJEMPLOS
				
#  ---------------------------------------- -- ----------------------------------------
"""


def cos_con_ruido(t, A, w, phi):
    """
    A: amplitud (número real)
    w: frecuencia (no negativa)
    phi: desfase (en [0,1[)
    t: argumento de la función
    """
    return A * math.cos(2*pi*w*t + 2*pi*phi) + np.random.uniform(-0.5, 0.5)

def cos_sin_ruido(t, A, w, phi):
    """
    A: amplitud (número real)
    w: frecuencia (no negativa)
    phi: desfase (en [0,1[)
    t: argumento de la función
    """
    return A * math.cos(2*pi*w*t + 2*pi*phi) 

def sinusoide_espectros(n, w, A, phi, nombre, ruido = True):
    """
    n: tamaño de la muestra
    w: frecuencia 
    A: amplitud,
    phi: (entre 0 y 1) desfase
    nombre: (str) nombre de la señal
    ruido == True sii se usa a la función cos_con_ruido para muestrear, en caso contrario, se usa cos_sin_ruido
    """
    frecuencias = [a/100 for a in range(int(n*100/2) + 1)]
    if ruido == True : 
        x = [cos_con_ruido(m/n, A, w, phi) for m in range(n)]
    else:
        x = [cos_sin_ruido(m/n, A, w, phi) for m in range(n)]

    return ae.analisis_espectrales_mostrarGrafica(x, nombre)
    
    """
#  ---------------------------------------- -- ----------------------------------------

				         EJECUCIONES

#  ---------------------------------------- -- ----------------------------------------
"""


def f(t):
      return 3* np.sin(2*np.pi*t) + np.sin(2*np.pi*4*t) + 0.5* np.cos(2*np.pi*7*t)

x = [f(t/25) for t in range(25)]
#Fs = 25
#ts = 1/Fs
#t = np.arange(0,1,ts)
#x = f(t) #muestreando y guardando los resultados, dando lugar a la señal 'x'.
#frecuencias = [a/100 for a in range(int(25*100/2) + 1)]
##ae.grafica_analisisEspectrales(x, 'x')
#
#
#n, w, A, phi, nombre = 16, 3, 2.3, 0, 'x'

#x = [cos_sin_ruido(t/36, -1.5, 3.4, 0.2) for t in range(36)]
#x = [cos_sin_ruido(t/36, -1.5, 5, 0.2) for t in range(36)]
#x = [1,2,-3,-6, -9, 5, 8, 7, 7, 2.3, 6.9]



# -------------- Ejemplo alternado 0 ----------------
#x = [cos_con_ruido(t/50,A = -1.5, w = 22.4, phi = 0.2) for t in range(50)]
#nombre = 'x'
#ae.grafica_analisisEspectrales(x, nombre, graficar = True) 
#x_alter = ae.alternado(x)
#nombre = r'A_{50}(x)'
#ae.grafica_analisisEspectrales(x_alter, nombre, graficar = True) 
# -------------- Fin ejemplo alternado 0 ----------------

# -------------- Ejemplo alternado 0 ----------------
x = [t/12 for t in range(12)] 
nombre = 'x'
#ae.grafica_analisisEspectrales(x, nombre, graficar = True)  #Sí funciona
x_alter = ae.alternado(x)
print(len(x_alter))
nombre = r'A_{12}(x)'
ae.grafica_analisisEspectrales(x_alter, nombre, graficar = True) #No funciona
# -------------- Fin ejemplo alternado 0 ----------------


#x = [1,2,3,4,6,5]
#ae.grafica_analisisEspectrales(x, nombre, graficar = True) 
#x = [1,-2,3,-4,6,-5]
#ae.grafica_analisisEspectrales(x, nombre, graficar = True) 
#
#x = np.array([12,4,3,-3.5,4.3,21,8, 9])
#norma = np.linalg.norm(x)
#x = (1/norma) * x


# -------------- Ejemplo alternado 1 ----------------
#nombre = r'\mathcal{L}^{28,5}'
#x = legendre.calculo_base(28)[5]
#ae.grafica_analisisEspectrales(x, nombre, graficar = True) 
#nombre = r'A_{28}(\mathcal{L}^{28,5})'
#x_alter = ae.alternado(x)
#ae.grafica_analisisEspectrales(x_alter, nombre, graficar = True) 
# -------------- Fin ejemplo alternado 1 ----------------



#
#x = legendre.calculo_base(12)[1]
#ae.grafica_analisisEspectrales(x, nombre, graficar = True) 
#
#x = legendre.calculo_base(35)[1]
#ae.grafica_analisisEspectrales(x, nombre, graficar = True) 

#ae.grafica_analisisEspectrales_PDL(7,5)

#ae.grafica_nube_b0m0_b1m1() 

#ae.grafica_pendientes_oOrigen_RMC()

"""
#No ejecutar más:

for n in range(2, 40):
    for k in range(n):
        ae.grafica_analisisEspectrales_PDL(n, k, graficar = False) 

for n in range(2,70):
    ae.grafica_analisisGlobal_n_fijo(n, graficar = False)

for k in range(2,70):
    ae.grafica_analisisGlobal_k_fijo(k, graficar = False)

for k in range(10, 70):
    ae.grafica_coefEspectrales_de_frecPrincipales(k, False)
"""
 
 
 
 
        
