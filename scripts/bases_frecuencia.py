#Extendiendo el rango de la frecuencia omega (w).
#En realidad no es una base. Ver por qué los productos punto son una buena forma de conocer las frecuencias en la base de Legendre discreta.


import numpy as np #numpy tiene funciones seno y coseno.
import matplotlib.pyplot as plt
import math #para sqrt y ceiling function

pi=math.pi
# Definición de vectores ------------------------------------------------------

def vector_c_N_v(N, v):
  #El único caso especial: N par, v=N/2. TODO: comenta por qué se hace esto en la teoría.
  if N%2==0 and v==N/2:
    resultado= [(-1)**mu/math.sqrt(N) for mu in range(N)]

  if v==0:
    return (1/math.sqrt(N))*np.ones([N])
  else:
    resultado=np.empty(0) #inicializamos el vector
    for mu in range(N):
      resultado=np.append(resultado, [math.sqrt(2/N) * math.cos(2*pi*v*mu/N)])
  return resultado


def vector_s_N_v(N, v):
  resultado=np.empty(0) #inicializamos el vector
  for mu in range(N):
      resultado=np.append(resultado, [math.sqrt(2/N) * math.sin(2*pi*v*mu/N) ])
  return resultado