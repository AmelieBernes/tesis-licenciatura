"""

Sistema con frecuencias "personalizadas": aunque la señal a analizar tenga
dimensión n, no necesariamente se consideran n frecuencias para formar este sistema.
Además, las frecuencias usadas no necesariamente son enteras.

Se usa pues como argumento un array de frecuencias.

Estatus: EN CONSTRUCCIÓN.
"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import pylab

pi=math.pi

#Funciones coseno y seno a partir de las que se construye todo lo que sigue.

def c_w(n,t, w):
  return math.sqrt(1/n)*np.cos(2*pi*w*t)

def s_w(n,t, w):
  return math.sqrt(1/n)*np.sin(2*pi*w*t)	

#------------------------------------

def calculo_base(omegas): 
  """
  'omegas' es un array de frecuencias.
  """
  cant_frec=len(omegas) #cantidad de frecuencias a analizar
  dominio=[k/n for k in range(cant)]
  M=math.ceil(N/2) #cota superior de las frecuencias consideradas en la base
  

  base_F=[(1/math.sqrt(N))*np.ones([N])] #inicializamos la base, que será un array. Ya incluimos la primera entrada.

  for w in range(1,M): 
    f_w=[]
    g_w=[]
    for t in dominio:
      f_w.append(math.sqrt(2)*c_w(N, t, w))
      g_w.append(math.sqrt(2)*s_w(N, t, w))
    base_F.append(f_w)
    base_F.append(g_w)

  if N%2==1: #si N es impar, ya terminamos
    return base_F 
  else: #en caso contrario, falta agregar un vector con una frecuencia más alta
    f_w=[]
    for t in dominio:
      f_w.append(math.sqrt(2)*c_w(N, t, omega[cant_frec-1]))
    base_F.append(f_w) 
    return base_F
