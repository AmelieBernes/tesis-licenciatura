import numpy as np #numpy tiene funciones seno y coseno.
import matplotlib.pyplot as plt
import math 

"""

Programa hecho para graficar, de cuatro en cuatro, los elementos de 
bases de R^{n}; se planea usar este script para graficar

1.- Bases de Frecuencia (i.e. de Fourier),
2.- Bases de Legendre discretas

Estas bases se calculan en los scripts
1.0 base_fourier_v0
1.1 base_fourier_v1
2. base_legendreDiscreta

En todos estos scripts se ha definido una función llamada 'calculo_base', cuyo
entero es un int n mayor a uno y cuyo output es un array con n arrays, cada uno
conteniendo la información de los vectores de la base cuyo nombre se establece en
el título del script.


#Creé una rama en Github para seguir editando esto. Veamos si funciona como creo.
"""

import base_fourier_V0
import base_fourier_V1
import base_legendreDiscreta



#TODO: graficar, para dimensiones pequeñas, las graficos
#Quiero graficar de cuatro en cuatro; cuando los cuatro axis de la figura se llenen, quiero guardarla
#y limpiarla para empezar de nuevo.


#---------------------------------------------------------------------------------------

def graficando_en_axis(vector, fig, ax):
  """
  wip
  """
  dominio=[t for t in range(len(vector))]
  ax.grid()
  return ax.scatter(dominio, vector, color='hotpink')


#---------------------------------------------------------------------------------------

#inicializamos la figura y cuatro axis, acomodados en un arreglo de 2x2
fig, axis=plt.subplots(2,2) 

#TODO: para que este script tenga sentido, debes de asegurarte que en TODOS los scripts en los que estas
#definiendo una base, la función que calcula la base se llame 'calculo_base'. 
def guardando_graficas(N, modulo_base=base_fourier_V0, nombre_base= 'Fourier (v0)' ,ruta=None):
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
				axis[i,j].set_title( nombre_base +": dimensión " +str(N)+", grado "+str(k)) #TODO: pedir la primera parte del título (Fourier, Legendre) como input
				k+=1 #aumentamos en uno la variable de grado
		#plt.legend() #TODO: por qué no corre el script con esta linea?


		#Guardando la imagen
		if ruta==None: 
			plt.savefig("AMELIA_1.png") #El archivo se guarda en la carpeta en la que está este script.
		else:
			plt.savefig(ruta+"AMELIA_1.png")
		#plt.show()
		plt.clf()

	#if N-N//4=!0: #si aún faltan grados que graficar, llenamos una última figura.



#TODO: hacer otra función 'imprimiendo gráficas' :)

if __name__ == "__main__":
	guardando_graficas(4, base_legendreDiscreta, 'LegDis' ,"/home/ame/GitHub/tesis-licenciatura/imagenes/")