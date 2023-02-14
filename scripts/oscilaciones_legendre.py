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
	Claro que se supone que las listas 'dominio' y 'mediciones' tienen la misma longitud,
	y que el i-ésimo valor de 'dominio' es mapeado al i-ésimo
	vector de 'mediciones'.'
	"""
	esp=0 #inicializamos la esperanza
	for i in range(len(dominio)):
		esp+=dominio[i]*mediciones[i]
	return esp

def calculando_sigmasYesp(N,k):
	"""
	Función que calcula los coeficientes sigma
	del polinomio discreto de Legendre de dimensión N
	y grado k, junto con la esperanza de la distribución
	que ellos forman.
	"""

	M=math.ceil(N/2)

	#inicializamos las bases de Fourier y Legendre discreta de dim N
	baseFourier=fourier.base_Fourier(N)
	baseLegendre=legendre.base_Legendre(N)

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

def graficando_sigmasYesp(N,k):
	"""
	Función que grafica los coeficientes sigma
	del polinomio discreto de Legendre de dimensión N
	y grado k, junto con la esperanza de la distribución
	que ellos forman.
	"""
	sigma, esp = calculando_sigmasYesp(N,k) #calculamos los datos
	dominio=[t for t in range(len(sigma))]  #calculamos len(sigma) para no tener que hacer un if-else con la paridad.

	#graficando las sigmas
	plt.scatter(dominio, sigma, s=100, color="mediumpurple", marker="*")

	#graficando la esperanza (un solo punto).
	plt.scatter(esp, 0, s=100, color="darkgoldenrod", marker="^") #elegí la forma de una cuña (esperanza como punto de equilibrio)

	plt.grid()
	plt.axhline(y=0, color='gray')	
	plt.axvline(x=0, color='gray')

	plt.title("Coeficientes sigma del pol. discrete de Legendre de dim. "+str(N)+" y grado "+str(k))


def graficando_esperanzas(N):
	"""
	Función que calcula las esperanzas de los coeficientes sigma de cada uno de los N polinomios
	discretos de Legendre de grado N.
	"""
	#Cambia el nombre del eje horizontal a k.
	baseFourier=fourier.base_Fourier(N)
	baseLegendre=legendre.base_Legendre(N)

	dominio=[t for t in range(N)] 
	esperanzas=[calculando_sigmasYesp(N,k)[1] for k in range(N)] #iteramos en la variable de grado 'k'.

	plt.scatter(dominio, esperanzas, s=100, color="darkgoldenrod", marker="^")

	X=np.linspace(0, N, 100)
	plt.plot(X, X/2, color="black", linestyle='dashed', label="Gráfica de la recta $y=\\frac{1}{2}k$")

	plt.xlabel("Grado k")
	plt.ylabel("Esperanza de la distribución NAME")
	plt.legend()
	plt.grid()
	plt.axhline(y=0, color='gray')	
	plt.axvline(x=0, color='gray')

	plt.title("Esperanzas de las distribuciones sigma de los pol. de Legendre de dim. "+str(N))

#graficando_sigmasYesp(50,25)
#graficando_esperanzas(20)
graficando_sigmasYesp(60,5)
#graficando_esperanzas(80)
plt.show()


#plt.clf()
#graficando_sigmasYesp(10,2)
#graficando_sigmasYesp(10,3)
#plt.show()

#for k in range(10):
#	fig, ax =plt.subplots(1,1)
#	graficando_sigmasYesp(10,k, fig, ax)
#	plt.savefig("Desktop/Amelie/Dimension10_Grado{}".format(str(k))
#	plt.close()


#--------------------------------------


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



#def graficando_esperanzas(N):
#	"""
#	Función que calcula las esperanzas de las N distribuciones de los coeficientes al cuadrado
#	correspondientes a cada grado 0 leq k leq N-1. 
#	"""
#	baseFourier=fourier.base_Fourier(N)
#	baseLegendre=legendre.base_Legendre(N)
#	dominio=[t for t in range(N)]

#	esperanzas=[] #inicializando el vector de esperanzas
#	for k in range(N): #iterando en la variable de grado
#		vectorLegendre=baseLegendre[k]
#		coeficientes_Fourier=[np.dot(vectorLegendre, baseFourier[v])**2 for v in range(N)]
#		esperanza_coefFourier=esperanza(dominio, coeficientes_Fourier)
#		esperanzas.append(esperanza_coefFourier)

#	plt.scatter(dominio, esperanzas, s=100, color="mediumpurple", marker="*")
#	plt.grid()
#	plt.axhline(y=0, color='black')
#	plt.axvline(x=0, color='black')

#	plt.title("Esperanzas de las distribuciones para $N= $"+str(N))
#	#Agregar un título!
#	plt.show()