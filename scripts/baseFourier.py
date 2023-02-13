import numpy as np #numpy tiene funciones seno y coseno.
import math #para sqrt y ceiling function

pi=math.pi
# Definición de vectores ------------------------------------------------------

def vector_c_N_v(N, v):
  if v==0:
    return (1/math.sqrt(N))*np.ones([N])
  else:
    resultado=np.empty(0) #inicializamos el vector
    for mu in range(N):
      resultado=np.append( (math.sqrt(2/N)) * math.cos(2*pi*v*mu/N), resultado)
  return resultado


def vector_s_N_v(N, v):
  resultado=np.empty(0) #inicializamos el vector
  for mu in range(N):
      resultado=np.append( (math.sqrt(2/N)) * math.sin(2*pi*v*mu/N), resultado)
  return resultado


# Base de Fourier de dimensión N ------------------------------------------------

def base_Fourier(N):
  M=math.ceil(N/2)
  base_F=[vector_c_N_v(N, 0)] #inicializando la base, que será un array. Ya incluimos la primera entrada
  for v in range(1,M):
    base_F.append(vector_c_N_v(N, v))
    base_F.append(vector_s_N_v(N, v))
  if N%2==1: #si N es impar, ya terminamos
    return base_F
  else: #en caso contrario, falta agregar un último vector a la base.
    base_F.append(vector_c_N_v(N, M))
    return base_F



#Comprobando ortonormalidad ---------------------------------------------------------

def comprobar_BON(N, i, j):
  """
  N es un entero positivo. i y j son enteros no negativos menores a N.
  Se regresa el producto punto entre el i-ésimo y el j-ésimo elemento
  de la base de Fourier F de dimensión N.
  """
  base=base_Fourier(N)
  print(np.dot(base[i], base[j]))