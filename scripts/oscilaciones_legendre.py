import math 
import numpy as np
import matplotlib.pyplot as plt
import pylab

import base_legendreDiscreta as legendre
import base_fourier_V0 
import base_fourier_V1 
import proyecciones #para calcular rectas de mínimos cuadrados.

#TODO: recuerda que, para poder importar mejor, en un script en el que sólo quieres definir funciones debes
#evitar ejecutar cosas. UPDATE: acabo de descubrir el idiom if __name__ == "__main__".

#TODO: como ahora tengo más de una versión de base de Fourier (base de frecuencias), tengo que poner a esta
#como argumento en las funciones de abajo.

#-----------------------------------------------------------------------------------

#Función auxiliar
def esperanza(dominio, mediciones):
	"""
	Función que calcula la esperanza de la variable aleatoria cuyo dominio es la lista
	'dominio' y que toma los valores del vector 'mediciones'.
	Claro que se supone que las listas 'dominio' y 'mediciones' tienen la misma longitud,
	y que el i-ésimo valor de 'dominio' es mapeado al i-ésimo
	vector de 'mediciones'.'
	"""
	esp=0 #inicializamos la esperanza
	for i in range(len(dominio)):
		esp+=dominio[i]*mediciones[i]
	return esp

#-----------------------------------------------------------------------------------
#Funciones principales del script
def calculando_sigmasYesp(N,k, fourier=base_fourier_V1):
	"""
	Función que calcula los coeficientes sigma
	del polinomio discreto de Legendre de dimensión N
	y grado k, junto con la esperanza de la distribución
	que ellos forman.
	"fourier" es un script en el que se ha definido una base de fourier
	(i.e. una base ortonormal de frecuencias).
	"""

	M=math.ceil(N/2)

	#inicializamos las bases de Fourier y Legendre discreta de dim N
	baseFourier=fourier.calculo_base(N)
	baseLegendre=legendre.calculo_base(N)

	#Llamamos al vector de Legendre de dim N y grado k
	vectorLegendre= baseLegendre[k]

	#Creando el vector de productos punto al cuadrado.
	prod_punto_cuadrado=[np.dot(vectorLegendre, baseFourier[mu])**2 for mu in range(N)] 

	#a partir de los productos punto al cuadrado almacenados en la lista de arriba calculamos los coeficientes sigma
	sigma=[prod_punto_cuadrado[0]] #inicializamos la lista con la primera entrada
	for l in range(1,M):
		sigma.append(prod_punto_cuadrado[2*l-1]+prod_punto_cuadrado[2*l])

	if N%2==1: #si N es impar, calculamos la esperanza y ya terminamos.
		dominio=[t for t in range(M)] #TODO: No sería mejor empezar desde uno? Para que la primera frecuencia no
									  #se elimine al calcular la esperanza...
		esp= esperanza(dominio, sigma)
		return (sigma, esp) 
	else: #en caso contrario, falta agregar una sigma y calcular la esperanza
		sigma.append(prod_punto_cuadrado[N-1])
		dominio=[t for t in range(M+1)] 
		esp= esperanza(dominio, sigma)
		return (sigma, esp) #para calcular esp necesito a sigma, por lo que no creo poder separar estos outputs

def graficando_sigmasYesp(N,k, fourier=base_fourier_V1):
	"""
	Función que grafica los coeficientes sigma
	del polinomio discreto de Legendre de dimensión N
	y grado k, junto con la esperanza de la distribución
	que ellos forman.
	"""
	sigma, esp = calculando_sigmasYesp(N,k, fourier) #calculamos los datos
	dominio=[t for t in range(len(sigma))]  #calculamos len(sigma) para no tener que hacer un if-else con la paridad.

	#graficando las sigmas
	plt.scatter(dominio, sigma, s=100, color="mediumpurple", marker="*")

	#graficando la esperanza (un solo punto).
	plt.scatter(esp, 0, s=100, color="darkgoldenrod", marker="^", label='Esperanza: '+str(esp.round(4))) #elegí la forma de una cuña (esperanza como punto de equilibrio)

	plt.xlabel("Frecuencia $\\omega$")
	plt.ylabel(r"$\sigma_{{\omega}}^{{ {0} }}( \mathcal{{ L }}^{{ {0} , {1} }} )$".format(str(N), str(k)) )
	plt.grid()
	plt.legend()
	plt.axhline(y=0, color='gray')	
	plt.axvline(x=0, color='gray')
	plt.title("La distribución $\\sigma_{{ {0} , {1} }}$ y su esperanza".format(str(N), str(k)) )
	


def graficando_esperanzas(N, fourier=base_fourier_V1):
	"""
	Función que calcula las esperanzas de los coeficientes sigma de cada uno de los N polinomios
	discretos de Legendre de grado N.
	
	NOTA: Observa cómo el primer punto siempre es cero. Esto se corresponde con el hecho
	de que todo polinomio discreto de Legendre de grado cero no tiene oscilaciones.
	TODO: fue bueno usar numpy arrays en lugar de arrays. Tienes que hacer esto para las demás funciones también.
	"""

	baseFourier=fourier.calculo_base(N)
	baseLegendre=legendre.calculo_base(N)

	dominio=np.array([t for t in range(N)])
	esperanzas=np.array([calculando_sigmasYesp(N,k, fourier)[1] for k in range(N)]) #iteramos en la variable de grado 'k'.

	plt.scatter(dominio, esperanzas, s=100, color="darkgoldenrod", marker="^")

	#Graficando la recta f(k)=k/2
	X=np.linspace(0, N, 100)
	plt.plot(X, X/2, color="black", linestyle='dashed', label="Gráfica de la recta $y=\\frac{1}{2}k$")
	b0, b1= proyecciones.coef_RMC(dominio, esperanzas)
	plt.plot(X, b1*X+b0, color="mediumblue", linestyle='dashed', label='Ajuste lineal de mínimos cuadrados')

	plt.xlabel("Grado $k$")
	plt.legend()
	plt.grid()
	plt.axhline(y=0, color='gray')	
	plt.axvline(x=0, color='gray')

	plt.title("Esperanzas de las distribuciones sigma de los pol. de Legendre de dimensión {0}".format(str(N)))


if __name__ == "__main__":
	#graficando_esperanzas(90, base_fourier_V0) 
	graficando_sigmasYesp(100,15, base_fourier_V0) 
	plt.show()


#Intento fallido de guardar en lugar de mostrar gráficas.
#for k in range(10):
#	fig, ax =plt.subplots(1,1)
#	graficando_sigmasYesp(10,k, fig, ax)
#	plt.savefig("Desktop/Amelie/Dimension10_Grado{}".format(str(k))
#	plt.close()



#-------------------------------------------------------------------------------------------------------#

                          #Funciones viejas: No sé si sea buena idea usarlas!

#-------------------------------------------------------------------------------------------------------#


##Calculamos la base de Fourier de dimensión n
#baseFourier=fourier.base_Fourier(N)
##Calculamos la base de Legendre discreta de dimensión n
#baseLegendre=legendre.base_Legendre(N)
#dominio=[t for t in range(N)]


#def graficando_coefFourier(N,k):
#	"""
#	N es un entero positivo, k es un entero entre 0 y N-1 (inclusivo)
#	Se calculan los cuadrados de los coeficientes de la señal discreta de legendre de dimensión N y grado k.
#	"""
#	vectorLegendre=baseLegendre[k] #extraemos el vector de nuestro interés; el de grado k.
#	coeficientes_Fourier=[np.dot(vectorLegendre, baseFourier[v])**2 for v in range(N)]
#	esperanza_coefFourier=esperanza(dominio, coeficientes_Fourier)
#	return (coeficientes_Fourier, esperanza_coefFourier)


#(coeficientes_Fourier, esperanza_coefFourier)=graficando_coefFourier(N,k)

#plt.axhline(y=0, color='black')
#lt.axvline(x=0, color='black')


#plt.scatter(dominio, coeficientes_Fourier, s=200)
#plt.scatter(esperanza_coefFourier, 0, s=200, color='red')
#plt.axvline(x=esperanza_coefFourier, color='red', linestyle='dotted')

#plt.grid()
#plt.show()


def graficando_esperanzas(N):
	"""
	Función que calcula las esperanzas de las N distribuciones de los coeficientes al cuadrado
	correspondientes a cada grado 0 leq k leq N-1. 
	TODO: Esta es una función antigua, no sé qué tan bueno sea usarla.
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
	plt.title("Esperanzas de las distribuciones $\\sigma$ de los pol. de Legendre discretos de dimensión {dimen}".format(dimen=str(N)))
	plt.xlabel("Degree")
	plt.ylabel("Means")
	#plt.title("Esperanzas de las distribuciones para $N= $"+str(N))
    
	#Agregar un título!
	plt.show()


#graficando_esperanzas(15)
