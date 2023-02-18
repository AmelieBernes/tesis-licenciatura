import matplotlib.pyplot as plt
import numpy as np #numpy tiene funciones seno y coseno.
import math #para sqrt y ceiling function

"""
Aquí voy a definir la segunda base de Fourier (BON de oscilaciones). 
"""

#NOTA: El contenido de este script ya no coincide con su título

#TODO: graficar, para dimensiones pequeñas, las graficos
#Quiero graficar de cuatro en cuatro; cuando los cuatro axis de la figura se llenen, quiero guardarla
#y limpiarla para empezar de nuevo.


fig, axis=plt.subplots(2,2) 

def imprimiendo_graficas(N):
	"""
	N: tipo int. Dimensión de la base de Fourier que se va a graficar.
	Por el momento, sólo contemplo N mútliplo de 4.
	"""
	k=0
	for u in range(N//4): #Creando, de ser posible, imágenes con los cuatro axis llenos.
		for i in range(2):
			for j in range(2):
				baseFourier_V0.graficando_Fourier_v0(N, k, fig, axis[i,j])
				k+=1 #aumentamos en uno la variable de grado
		plt.legend()
		plt.show()
		plt.clf()

	#if N-N//4=!0: #si aún faltan grados que graficar, llenamos una última figura.


if __name__ == "__main__":
	imprimiendo_graficas(4)