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
#Sí funciona! La rama en la que estoy trabajando en la edición de este script se llama
#'script_graficas'.
"""

import base_fourier_V0
import base_fourier_V1
import base_legendreDiscreta



#TODO: graficar, para dimensiones pequeñas, las graficos
#Quiero graficar de cuatro en cuatro; cuando los cuatro axis de la figura se llenen, quiero guardarla
#y limpiarla para empezar de nuevo.

#TODO: study the difference between window, figure and axis. Recuerda que en matplotlib
#importa la jerarquía!


#Intento de cambiar el tipo de fuente

plt.style.use('seaborn-v0_8-poster') 
#plt.style.use('seaborn-v0_8-pastel') 
params = {"ytick.color" : "black",
          "xtick.color" : "black",
          "axes.labelcolor" : "black",
          "axes.edgecolor" : "black",
          "text.usetex" : True,
          "font.family" : "serif",
          "font.serif" : ["Computer Modern Serif"]}
plt.rcParams.update(params)

#---------------------------------------------------------------------------------------



def graficando_en_axis(vector, fig, ax):
  """
  wip
  """
  dominio=[t for t in range(len(vector))]
  ax.grid()
  return ax.scatter(dominio, vector, color='hotpink')

def guardando_cuatro_axis(fig, axis, N, k, u, base, nombre_base, ruta):
	for i in range(2):
		for j in range(2):
			vector=base[k]
			graficando_en_axis(vector, fig, axis[i,j])
			axis[i,j].set_title( nombre_base +": dimensión " +str(N)+", grado "+str(k)) #TODO: pedir la primera parte del título (Fourier, Legendre) como input
			axis[i,j].axhline(y=0, color='gray')
			axis[i,j].axvline(x=0, color='gray')
			k+=1 #aumentamos en uno la variable de grado
	#Guardando la imagen
	nombre_imagen= nombre_base +": dimensión"+str(N)+", figura "+str(u)
	fig.suptitle(nombre_imagen)
	if ruta==None: 
		plt.savefig(nombre_imagen) #El archivo se guarda en la carpeta en la que está este script.
	else:
		plt.savefig(ruta+nombre_imagen)

def guardando_cuatro_axis_finales(fig, axis, N, k, u, base, nombre_base, ruta):
	K=0
	for i in range(2):
		if K==N-1: #Si ya llegamos al último grado 
			break #salimos
		for j in range(2):
			vector=base[k]
			graficando_en_axis(vector, fig, axis[i,j])
			axis[i,j].set_title( nombre_base +": dimensión " +str(N)+", grado "+str(k)) #TODO: pedir la primera parte del título (Fourier, Legendre) como input
			axis[i,j].axhline(y=0, color='gray')
			axis[i,j].axvline(x=0, color='gray')
			K+=1 #aumentamos en uno la variable de grado
			if K==N-1: #si ya graficamos el último elemento de la base...
				break #salimos del ciclo for.
	#Guardando la imagen
	nombre_imagen= nombre_base +": dimensión"+str(N)+", figura "+str(u+1)
	fig.suptitle(nombre_imagen)
	if ruta==None: 
		plt.savefig(nombre_imagen) #El archivo se guarda en la carpeta en la que está este script.
	else:
		plt.savefig(ruta+nombre_imagen)


#---------------------------------------------------------------------------------------
 
def guardando_graficas(N, modulo_base=base_fourier_V0, nombre_base= 'Fourier (v0)' ,ruta=None):
	"""
	N: tipo int, mayor a uno. Representa una dimensión.

	'funcion_base' es un módulo en el que se ha definido una función llamada 'calculo_base', función
	que calcula la base cuyos elementos se quieren graficar (output: np.array o array. 
	#TODO: homogeiniza esto). 

	#NOTA: por el momento esta función sólo funciona para cuando N=4 (también cuando N es múltiplo de 4??)
	"""
	q=N//4
	r=N%4 #(por lo tanto, N=q*4+r)
	base=modulo_base.calculo_base(N)
	k=0
	for u in range(q): #Creando, de ser posible, imágenes con los cuatro axis llenos.

		#inicializamos la figura y cuatro axis, acomodados en un arreglo de 2x2
		fig, axis=plt.subplots(2,2) 
		fig.tight_layout(pad=3.0) #para cambiar el espacio entre los axis de la figura.
		guardando_cuatro_axis(fig, axis, N, k, u, base, nombre_base, ruta)
		k+=4
		plt.close() #para cerrar la VENTANA! Por lo tanto, la figura y los axis.

	if r!=0: #si aún faltan grados que graficar, llenamos una última figura. Nota que $r \in \{ 1,2,3 \}$.
		fig, axis=plt.subplots(2,2)  #inicializamos pues la figura y cuatro axis en ella. Seguro no los llenamos todos.
		guardando_cuatro_axis_finales(fig, axis, N, k, u, base, nombre_base, ruta)





#TODO: hacer otra función 'imprimiendo gráficas' :)

if __name__ == "__main__":
	guardando_graficas(15, base_fourier_V0, 'Fourier (v0)',"/home/ame/GitHub/tesis-licenciatura/imagenes/bases/")