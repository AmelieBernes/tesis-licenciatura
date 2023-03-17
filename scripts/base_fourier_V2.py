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

def calculo_base(n, omegas): 
  """
  'omegas' es un array de frecuencias.
  En esta función 'calculo_base' también se requiere dar como argumento a la dimensión del array x cuyas
  frecuencias nos interesa estudiar, pues
  los vectores de frecuencias que formemos deben tener la misma dimensión de x para que tenga sentido formar productos punto.
  """

  dominio=[k/n for k in range(n)]

  base_F=[(1/math.sqrt(N))*np.ones([N])] #inicializamos la base, que será un array. Ya incluimos la primera entrada.

  for w in omegas: 
    f_w=[]
    g_w=[]
    for t in dominio:
      f_w.append(math.sqrt(2)*c_w(n, t, w))
      g_w.append(math.sqrt(2)*s_w(n, t, w))
    base_F.append(f_w)
    base_F.append(g_w)
    
    return base_F 
