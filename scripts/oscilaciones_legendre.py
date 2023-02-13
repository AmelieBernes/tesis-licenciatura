#amelie

import math 
import numpy as np
import matplotlib.pyplot as plt
import pylab

import calculo_basesLegendre_NUEVO as legendre
import baseFourier as fourier

#--------------------------------------------------------

def esperanza(dominio, mediciones):
	"""
	Función que calcula la esperanza de la variable aleatoria cuyo dominio es la lista
	'dominio' y que toma los valores del vector 'mediciones'.
	Claro que se supone que las listas 'dominio' y 'mediciones' tienen la misma longitud.
	"""
	esp=0 #inicializamos la esperanza
	for i in range(len(dominio)):
		esp+=dominio[i]*mediciones[i]
	return esp


def graficando_esperanzas(N):
	"""
	Función que calcula las esperanzas de las N distribuciones de los coeficientes al cuadrado
	correspondientes a cada grado 0 leq k leq N-1. 
	"""
	baseFourier=fourier.base_Fourier(N)
	baseLegendre=legendre.base_Legendre(N)
	dominio=[t for t in range(N)]

	esperanzas=[] #inicializando el vector de esperanzas
	for k in range(N): #iterando en la variable de grado
		vectorLegendre=baseLegendre[k]
		coeficientes_Fourier=[np.dot(vectorLegendre, baseFourier[v])**2 for v in range(N)]
		esperanza_coefFourier=esperanza(dominio, coeficientes_Fourier)
		esperanzas.append(esperanza_coefFourier)

	plt.scatter(dominio, esperanzas, s=100, color="mediumpurple", marker="*")
	plt.grid()
	plt.axhline(y=0, color='black')
	plt.axvline(x=0, color='black')

	plt.title("Esperanzas de las distribuciones para $N= $"+str(N))
	#Agregar un título!
	plt.show()


graficando_esperanzas(90)

#--------------------------------------

N=8
k=5

#Calculamos la base de Fourier de dimensión n
baseFourier=fourier.base_Fourier(N)
#Calculamos la base de Legendre discreta de dimensión n
baseLegendre=legendre.base_Legendre(N)
dominio=[t for t in range(N)]


def graficando_coefFourier(N,k):
	"""
	N es un entero positivo, k es un entero entre 0 y N-1 (inclusivo)
	Se calculan los cuadrados de los coeficientes de la señal discreta de legendre de dimensión N y grado k.
	"""
	vectorLegendre=baseLegendre[k] #extraemos el vector de nuestro interés; el de grado k.
	coeficientes_Fourier=[np.dot(vectorLegendre, baseFourier[v])**2 for v in range(N)]
	esperanza_coefFourier=esperanza(dominio, coeficientes_Fourier)
	return (coeficientes_Fourier, esperanza_coefFourier)


(coeficientes_Fourier, esperanza_coefFourier)=graficando_coefFourier(N,k)

plt.axhline(y=0, color='black')
plt.axvline(x=0, color='black')


plt.scatter(dominio, coeficientes_Fourier, s=200)
plt.scatter(esperanza_coefFourier, 0, s=200, color='red')
plt.axvline(x=esperanza_coefFourier, color='red', linestyle='dotted')

#plt.grid()
#plt.show()