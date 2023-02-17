"""

Programa hecho para graficar, de cuatro en cuatro, los elementos de 
bases de R^{N}. Se planea usar este script para graficar bases
de frecuencias y las de Legendre discretas.

"""

#TODO: graficar, para dimensiones pequeñas, las graficos
#Quiero graficar de cuatro en cuatro; cuando los cuatro axis de la figura se llenen, quiero guardarla
#y limpiarla para empezar de nuevo.


#inicializamos la figura y cuatro axis, acomodados en un arreglo de 2x2
fig, axis=plt.subplots(2,2) 

def imprimiendo_graficas(N, generar_base=baseFourier_V0.graficando_Fourier_v0()):
	"""
	N: tipo int, mayor a uno. Representa una dimensión.

	'generar_base' es una función que calcula los elementos de la base
	que en esta función se desean graficar. Se supone que el único argumento de esta
	función es la dimensión N del espacio.
	"""
	base=generar_base(N)
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