import matplotlib.pyplot as plt
import numpy as np #numpy tiene funciones seno y coseno.
from numpy.linalg import norm
import math #para sqrt y ceiling function

pi=math.pi

"""
Aquí voy a definir la segunda base de Fourier (BON de oscilaciones). 
"""



def vector_g_w_N(N, w):
	"""
	'N' es la dimensión (int mayor a uno), 'w' es una frecuencia (float positivo).
	"""
	#Revisamos el único caso fácil de calcular:
	if w==0:
		return 1/math.sqrt(N)*np.ones(N)
	#Creamos la malla uniforme
	malla_uniforme=[ -1/2+mu/(N-1) for mu in range (N)]
	G_w_N=np.array([1/math.sqrt(N)*math.cos(2*pi*w*t) for t in malla_uniforme])
	norma=norm(G_w_N)
	return 1/norma*G_w_N

def vector_h_w_N(N, w):
	"""
	'N' es la dimensión (int mayor a uno), 'w' es una frecuencia (float positivo).
	"""
	#Creamos la malla uniforme
	malla_uniforme=[ -1/2+mu/(N-1) for mu in range (N)]

	G_w_N=np.array([1/math.sqrt(N)*math.sin(2*pi*w*t) for t in malla_uniforme])
	norma=norm(G_w_N)
	return 1/norma*G_w_N


def calculo_base(N):
  M=math.ceil(N/2)
  base_F=[vector_g_w_N(N, 0)] #inicializando la base, que será un array. Ya incluimos la primera entrada
  for w in range(1,M):
    base_F.append(vector_g_w_N(N, w))
    base_F.append(vector_h_w_N(N, w))
  if N%2==1: #si N es impar, ya terminamos
    return base_F
  else: #en caso contrario, falta agregar un último vector a la base.
    base_F.append(vector_g_w_N(N, M))
    return base_F


if __name__ == "__main__":
	print(calculo_base(4))
