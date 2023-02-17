import numpy as np #numpy tiene funciones seno y coseno.
import matplotlib.pyplot as plt
import math 
import baseFourier_V0

"""

Programa hecho para graficar, de cuatro en cuatro, los elementos de 
bases de R^{N}. Se planea usar este script para graficar bases
de frecuencias y las de Legendre discretas.

"""

#TODO: graficar, para dimensiones pequeñas, las graficos
#Quiero graficar de cuatro en cuatro; cuando los cuatro axis de la figura se llenen, quiero guardarla
#y limpiarla para empezar de nuevo.


#------------------------------------------------------------------------------

def graficando_en_axis(vector, fig, ax):
  """
  wip
  """
  dominio=[t for t in range(len(vector))]
  ax.grid()
  return ax.scatter(dominio, vector, color='hotpink')



#inicializamos la figura y cuatro axis, acomodados en un arreglo de 2x2
fig, axis=plt.subplots(2,2) 

#TODO: para que este script tenga sentido, debes de asegurarte que en TODOS los scripts en los que estas
#definiendo una base, la función que calcula la base se llame 'calculo_base'. 
def imprimiendo_graficas(N, modulo_base=baseFourier_V0):
	"""
	N: tipo int, mayor a uno. Representa una dimensión.

	'funcion_base' es un módulo en el que se ha definido una función llamada 'calculo_base', función
	que calcula la base cuyos elementos se quieren graficar (output: np.array o array. #TODO: homogeiniza esto). 
	"""
	base=modulo_base.calculo_base(N)
	k=0
	for u in range(N//4): #Creando, de ser posible, imágenes con los cuatro axis llenos.
		for i in range(2):
			for j in range(2):
				vector=base[k] #extrayendo el k-ésimo elemento de la base.
				graficando_en_axis(vector, fig, axis[i,j])
				axis[i,j].grid()
				axis[i,j].set_title("Dimensión " +str(N)+", grado "+str(k)) #TODO: pedir la primera parte del título (Fourier, Legendre) como input
				k+=1 #aumentamos en uno la variable de grado
		plt.legend()
		plt.show()
		plt.clf()

	#if N-N//4=!0: #si aún faltan grados que graficar, llenamos una última figura.


if __name__ == "__main__":
	imprimiendo_graficas(4)