"""
Script en el que se realizan todos los cálculos
para realizar análisis espectrales de señales finitas, 
en particular, análisis espectrales de los PDL.

Este es el script que se importa en GUI_PDL.py

Este script se apoya en las funciones escritas en 
analisis_espectrales; decidí no usar este último para
programar la GUI pues los formatos que necesito para las 
gráficas son ligeramente distintos.
"""


 # --------------------------- Módulos -----------------------------------------
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
import analisis_espectrales as ae #de aquí voy a llamar a las funciones para realizar cálculos.

# --------------------------- Algunos objetos -------------------------------------

pi=math.pi
colores=['#fe01b1', '#3b2747', '#feb308', '#6832e3', '#feb308', '#8f99fb', 'gray', '#8e82fe', '#fd183d']

"""
El orden de los colores es:
    0: color de la señal a analizar
    1: color para los puntos del espectro que no son los de FP (frecuencia máxima) YA NO
    2: Frecuencia máxima del estudio basado en espacios monofrecuenciales
    3: Frecuencia máxima del estudio basado en la DFT
    4: Recta de mínimos cuadrados del estudio global de las FP encontradas con espacios monofrecuenciales
    5: Recta de mínimos cuadrados del estudio global de las FP encontradas con la DFT
    6: Color para los puntos del espectro monofrecuencial que no son máximos
    7: Color para los puntos del espectro de la DFT que no son máximos
    8: Color para lineas de referencia
"""


# ------------------------- Funciones de cálculos -------------------------------------



def grafica_taus_axis(x, n, nombre, axis1, axis2, leyenda_interior = False):
  """
  'x' es un array de dimensión mayor a dos. 
  Esta función dibuja la gráfica de 'x' con dominio el tiempo, junto
  con la gráfica de los coeficientes tau de x.

  En 'axis1' se grafica la señal x, #TODO ver si uso lo mismo para grafica sigmas.
  en 'axis2' se grafica el espectro.
  """
  M = math.ceil(n/2) #TODO poner como argumento
  coef_cosenos, coef_senos = ae.coeficientes_base_Fourier(x)
  
  dominio=[m/n for m in range(n)]
  taus = ae.coeficientes_tau(x)
  cant_freq=len(taus)

  axis1.scatter(dominio, x, color= colores[0], s=60, label= "${0}$".format(nombre), zorder = 3)

  for i in range(cant_freq):
    axis2.scatter(i, taus[i], color=colores[7], s=100, marker = '*', zorder = 2)

  axis2.set_xlabel('Frecuencias enteras $\omega$')
  axis2.set_ylabel(r'$\tau_{{{0}}}($'.format(str(n))+"${0}$".format(nombre)+ r'$, \omega)$' ) 

  X=np.arange(0, 1, 0.0001)
  
  max_w = taus.index(max(taus)) #máxima frecuencia
  axis2.scatter(max_w, max(taus), color = colores[3], s = 70, label = '( ' + str(max_w) + ', ' + str(round(max(taus), 4))  + ' )', marker = '^', zorder = 3)
  stemlines = axis2.stem(max_w, max(taus), markerfmt = ' ', linefmt = '--')
  plt.setp(stemlines, color = colores[3]) #artista Python

  if n %2 == 0: 
    if max_w == 0 or max_w == M:
      coef_cos = coef_cosenos[max_w] * math.sqrt(1/n)
      axis1.scatter(dominio, [coef_cos * np.cos(2*pi*max_w*t) for t in dominio], zorder=2, color = colores[3], s= 80)
      label = r'${{{0}}} \cdot cos(2 \pi \cdot {{{1}}} t) $'.format(str(round(coef_cos,2)), str(max_w))
      axis1.plot(X, coef_cos * np.cos(2*pi*max_w*X), color = colores[3])
    else:
      coef_cos = coef_cosenos[max_w] * math.sqrt(2/n)
      coef_sen = coef_senos[max_w-1] * math.sqrt(2/n)
      muestreo = [coef_cos * np.cos(2*pi*max_w*t)+coef_sen* np.sin(2*pi*max_w*t) for t in dominio]
      axis1.scatter(dominio, muestreo, zorder = 2, color = colores[3], s=80) 
      label = r'${{{0}}} \cdot cos(2 \pi \cdot {{{1}}} t) + {{{2}}} \cdot sen(2 \pi \cdot {{{1}}} t) $'.format(str(round(coef_cos,2)), str(max_w), str(round(coef_sen,2)))
      axis1.plot(X, coef_cos * np.cos(2*pi*max_w*X)  + coef_sen * np.sin(2*pi*max_w*X), color = colores[3])
  else: #o sea, si n%2 == 1
    if max_w == 0:
      coef_cos = coef_cosenos[max_w] * math.sqrt(1/n)
      axis1.scatter(dominio, [coef_cos * np.cos(2*pi*max_w*t) for t in dominio], zorder = 2, color = colores[3], s=80)
      label = r'${{{0}}} \cdot cos(2 \pi \cdot {{{1}}} t) $'.format(str(round(coef_cos,2)), str(max_w))
      axis1.plot(X, coef_cos * np.cos(2*pi*max_w*X), color = colores[3]) #redundante
    else: 
      coef_cos = coef_cosenos[max_w] * math.sqrt(2/n)
      coef_sen = coef_senos[max_w-1] * math.sqrt(2/n)
      axis1.scatter(dominio, [coef_cos * np.cos(2*pi*max_w*t)+coef_sen* np.sin(2*pi*max_w*t) for t in dominio], zorder=2, color = colores[3], s=80)
      label = r'${{{0}}} \cdot cos(2 \pi \cdot {{{1}}} t) + {{{2}}} \cdot sen(2 \pi \cdot {{{1}}} t) $'.format(str(round(coef_cos,2)), str(max_w), str(round(coef_sen,2)))
      axis1.plot(X, coef_cos * np.cos(2*pi*max_w*X) + coef_sen * np.sin(2*pi*max_w*X), color = colores[3])
  

  ae.formato_axis(axis1, leyenda_interior)
  ae.formato_axis(axis2)
  return label




def grafica_sigma_amplDesfase_axis_casosNoExtremos(x, w, axis):

  #TODO cambiar nombre a grafica_amplDesfase_axis_caso1
  n = len(x)

  #sigma = sigma_caso1(x, w) #TODO por qué calculo la sigma? esto no es necesario.
  A, phi = ae.amplDesfase_caso1(x, w)

  def coseno_amplDes(t):
    return A * np.cos(2*pi*w*t-2*pi*phi)

  dominio=[k/n for k in range(n)]
  proyeccion_Pnw = [coseno_amplDes(m/n) for m in range(n)]

  label=r'Gráfica de $\Pi_{{P_{{ {0}, {1} }} }}(x)$'.format(str(n),str( round(w, 2) ) )
  axis.scatter(dominio, proyeccion_Pnw, color=colores[2], s=80)
  
  X=np.arange(0, 1, 0.0001)
  axis.plot(X, coseno_amplDes(X), color=colores[2])

  ae.formato_axis(axis, leyenda_interior = False)
  return label


def analisis_espectral_espaciosMonofrecuenciales(x, n, frecuencias, nombre, axis0, axis1, leyenda_interior  = False):
  """
  x es un array
  'frecuencias' es un vector de frecuencias.
  En 'axis0' se grafica la señal x, en el 'axis1' el espectro.
  """
  dominio_tiempo=[m/n for m in range(n)] 

  axis0.scatter(dominio_tiempo, x, color= colores[0], s= 50)

  #Calculando el espectro
  sigmas = []
  for w in frecuencias:
    sigmas.append(ae.sigma_caso1(x,w))
  axis1.scatter(frecuencias, sigmas, color=colores[1], s=5, marker = '*', zorder = 2)

  sigmas.insert(0, ae.limite_cero(x,n))
  axis1.scatter(0, sigmas[0],color=colores[1], s=5, marker = '*', zorder = 2)
  sigmas.append(ae.limite_n_medios(x,n))
  axis1.scatter(n/2, sigmas[-1], color=colores[1], s=5, marker = '*', zorder = 2)
  #limite_0 = limite_cero(x, n) #NUEVO 
  #limite_n_2 = limite_n_medios(x, n) #NUEVO 
  sigma_max = max(sigmas)

  
  #Se tiene que agregar un ciclo if-else porque acorté el vector de frecuencias.
  # posiciones de sigmas =[0, 1, 2, ... , n*100/2]
  # posiciones de frecuencias = [1, 2, ..., n*100/2 -1]
  indice_sigmaMax = sigmas.index(sigma_max) #Mayor frecuencia en la que el espectro tiene su máximo
  if indice_sigmaMax == 0:
      frec_max = 0
  elif indice_sigmaMax == n*100/2:
      frec_max = n/2
  else:
    frec_max = frecuencias[indice_sigmaMax - 1] #TODO esto está bien ??

  axis1.scatter(frec_max, sigma_max, s = 100, color = colores[2], label = '( ' + str(round(frec_max, 2)) + ', ' + str(round(sigma_max, 2)) + ' )', marker = 'v', zorder=4)
  stemlines = axis1.stem(frec_max, sigma_max, markerfmt = ' ', linefmt = '--')
  plt.setp(stemlines, color = colores[2]) 
  
  if frec_max == 0:
    X = np.arange(0,1, 0.0001)
    proyeccion_W_n1 = proy.proyeccion(x,1)
    b0, b1 = proy.coef_RMC(dominio_tiempo, proyeccion_W_n1)
    label = r'Gráfica de $\Pi_{{ W_{{ {0}, 1 }} }}(x)$'.format(str(n))
    axis0.scatter(dominio_tiempo, proyeccion_W_n1, color=colores[2], s=80)
    axis0.plot(X, b0 + b1*X, color = colores[2])
    ae.formato_axis(axis0, leyenda_interior = False)
  elif frec_max == n/2:
    X = np.arange(0,1, 0.0001)
    x_alter = alternado(x)
    proyeccion_W_n1 = proy.proyeccion(x_alter,1)
    label = r'Gráfica de $\Pi_{{ W_{{ {0}, 1 }} }}(A_{{ {0}  }} (x))$'.format(str(n))
    axis0.scatter(dominio_tiempo, proyeccion_W_n1, color=colores[2], s=80)
    axis0.plot(dominio_tiempo, x_alter, color = colores[2])
    ae.formato_axis(axis0, leyenda_interior = False)
  else:    
    label = grafica_sigma_amplDesfase_axis_casosNoExtremos(x, frec_max, axis0)


  axis1.set_xlabel('Frecuencias ' + r'$\omega$')
  axis1.set_ylabel(r'$\sigma_{{{0}}}($'.format(str(n))+"${0}$".format(nombre)+ r'$, \omega)$' + r', $\tau_{{{0}}}($'.format(str(n))+"${0}$".format(nombre)+ r'$, \omega)$') 


  ae.formato_axis(axis0, leyenda_interior)
  ae.formato_axis(axis1)
  return label



# -----------------------------------------------------------------------------------------------------------

def grafica_analisisEspectrales_PDL_GUI(fig, n,k):
  """
  
  Función para graficar el análisis espectral de un polinomio discreto de Legendre (PDL).
  
  'n' (entero mayor o igual a dos) es la dimensión del PDL, 
  'k' \in \{0, 1, ... , n-1 \} es su grado.
  
  #TODO dar las especificaciones de la figura y sus axis.
  """
  #TODO nombre antiguo: analisis_espectrales_PDL_guardarGrafica(n,k)
  #TODO poner a la ruta como argumento. Como es siempre la misma ruta, ponla como una variable global,
  #para no tener que pedirla siempre como argumento.
 

  x = legendre.calculo_base(n)[k]
  frecuencias = [a/100 for a in range(int(n*100/2) + 1)] 
  #Quitamos las frecuencias extremas, que son casos especiales.
  frecuencias.pop(0)
  frecuencias.pop(-1)
  nombre = r'\mathcal{{L}}^{{{0}}}'.format(str(n)+','+str(k)) 

  #Recuperamos y limpiamos los axis de la figura
  ax_list = fig.axes
  ax0, ax1, ax2 = ax_list[0], ax_list[1], ax_list[2]
  ax0.clear()
  ax1.clear()
  ax2.clear()

  ax2.set_title('Espectros')
  ax2.axvline(x = k/2, color = colores[8], linewidth = 2, label = "x = " + str(k/2))
  ax2.axhline(y = 1, color = colores[8], linewidth = 2, linestyle = 'dotted')

  label_TDF =grafica_taus_axis(x, n, nombre, ax0, ax2)
  label_EM = analisis_espectral_espaciosMonofrecuenciales(x, n, frecuencias, nombre, ax1, ax2)

  return label_TDF, label_EM
  #return fig