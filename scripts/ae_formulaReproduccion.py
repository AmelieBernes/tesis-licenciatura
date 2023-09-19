#Planeo sustituir el script de analisis_espectrales (en el que las fórmulas aplicadas son correctas, pero pueden
#simplificarse con la teoría desarrollada en los últimos meses), pero por el momento en este script sólo voy a 
#escribir lo necesario para hacer los procesos de síntesis y análisis de una señal x in IRn 
#usando el espectro Tx de esta y la función de desfases Dx.


import math
import numpy as np
import matplotlib.pyplot as plt

pi = math.pi


def formato_axis(axis):
    """
    Agregando los elementos que me gustan a un axis
    """
    axis.axhline(y=0, color='gray', zorder = 0)
    axis.axvline(x=0, color='gray', zorder = 0)
    axis.grid(True)
    axis.legend(loc = 'best')



"""
#  ---------------------------------------- -- ----------------------------------------

				  CÁLCULO DE LA BASE DE FOURIER

#  ---------------------------------------- -- ----------------------------------------
"""


def calculo_baseFourier(n):
	"""
  	Cálculo de la base de Fourier de dimensión n.
	Resultados revisados, sí funciona correctamente.
	"""
	
	M = math.ceil((n-1)/2) #Cota superior de frecuencias
	rango_recuencias = [i for i in range(M+1)]
	dominio=[k/n for k in range(n)]

	#TODO el convertir los elementos de base_fourier es np.arrays me arruinó todo.
	base_fourier = [ [math.sqrt(1/n)*np.ones([n])] ] #inicializamos la base de Fourier con el primer vector

	for w in range(1,M): 
		c_w = np.array([math.sqrt(2/n)*np.cos(2*pi*w*t) for t in dominio])
		s_w = np.array([math.sqrt(2/n)*np.sin(2*pi*w*t) for t in dominio])
		base_fourier.append(c_w)
		base_fourier.append(s_w)

	
	if n %2 == 1: #n impar
		base_fourier.append(np.array([math.sqrt(2/n)*np.cos(2*pi*M*t) for t in dominio]))
		base_fourier.append(np.array([math.sqrt(2/n)*np.sin(2*pi*M*t) for t in dominio]))

	else: #n par
		base_fourier.append(np.array([math.sqrt(1/n)*np.cos(2*pi*M*t) for t in dominio]))


	return base_fourier



"""
#  ---------------------------------------- -- ----------------------------------------

				  PROCESO DE ANÁLISIS

#  ---------------------------------------- -- ----------------------------------------
"""


#Se ve bien.
def elementos_sintesis(x):

	n = len(x)
	M = math.ceil((n-1)/2)
	x = np.array(x) #Convirtiendo a 'x' (array) a np.array en caso de ser necesario.
	base_fourier = calculo_baseFourier(n)

	alphas = [np.dot(x, base_fourier[0])] #productos punto de x con c_nw
	betas = [0] #productos punto de x con s_nw
	gammas = [0] #productos punto de c_nw con s_nw

	for i in range(1, 2*M-2, 2):
		alphas.append(np.dot(x, base_fourier[i]))
		betas.append(np.dot(x, base_fourier[i+1]))
		gammas.append(np.dot(base_fourier[i], base_fourier[i+1]))


	if n % 2 == 1: #si n es impar
		alphas.append(np.dot(x, base_fourier[2*M-1]))
		betas.append(np.dot(x, base_fourier[2*M]))
		gammas.append(np.dot(base_fourier[2*M-1], base_fourier[2*M]))

	else: #si n es par
		alphas.append(np.dot(x, base_fourier[2*M-1]))
		betas.append(0)
		gammas.append(0)

	return alphas, betas, gammas



def espectro_y_desfases_fourier(x, n, M):
	alphas, betas, gammas = elementos_sintesis(x)
	norma_cuad = np.dot(x, x)
	espectro = [( (alphas[w]**2 + betas[w]**2)/ norma_cuad )**(1/2) for w in range(M+1) ]

	deltas = [(alphas[w] - gammas[w] * betas[w] ) / (1 - gammas[w]**2) for w in range(M+1)]
	epsilons = [(betas[w] - gammas[w] * alphas[w] ) / (1 - gammas[w]**2) for w in range(M+1)]

	desfases = [] #inicializamos el array de desfases
	for w in range(M+1): #iterando en las frecuencias
		c = deltas[w]
		d = epsilons[w]

		tan = math.atan(d/c)
		print("tan  " +  str(tan))

		if c >= 0 and d > 0:
			phi = tan/(2*pi)
		elif (c>0 and d < 0) or (c < 0 and d < 0):
			phi = (tan + pi)/ (2 * pi)
		elif (c == 0):
			if d ==1:
				phi = 1/4
			else:
				phi = 3/4
		elif (d == 0):
			if c == 1:
				phi =0
			else:
				phi = 1/2
		else:
			phi = (tan + 2 * pi) * (2 * pi) 
		print("phi = " + str(phi))
		desfases.append(phi)

	return (espectro, desfases)



def grafica_espectro_fourier(x):
	n = len(x)
	M = math.ceil((n-1)/2)

	dominio=[m/n for m in range(n)]
	fig, axis = plt.subplots(3,1)
	axis[0].scatter(dominio, x, color = 'green', label = "Gráfica de la señal")

	frecuencias = [w for w in range(M+1)]

	espectro = espectro_y_desfases_fourier(x, n, M)[0]
	desfases = espectro_y_desfases_fourier(x, n, M)[1] 

	print(desfases)

	axis[1].scatter(frecuencias, espectro, color = 'purple', label = "Espectro de la señal")
	axis[2].scatter(frecuencias, desfases, color = 'pink', label = "Función de desfases")


	formato_axis(axis[0])
	formato_axis(axis[1])
	formato_axis(axis[2])

	return plt.show()




"""
#  ---------------------------------------- -- ----------------------------------------

				  PROCESO DE SÍNTESIS

Objetivo: dadas las funciones T: DOM_n -----> [0,1] y D: DOM_n -----> [0,1[, encontrar
a la señal x in IRn de norma uno tal que T_x = T y D_x = D.

Tales funciones están dadas en formas de n-arrays, donde la k-ésima entrada es igual
al valor de la respectiva función en la frecuencia k.

#  ---------------------------------------- -- ----------------------------------------
"""

def formando_proyecciones(T, D):
	M = len(T)
	#TODO: código para checar que T y D tienen la misma longitud.
	proyecciones = [] 
	for w in range(M):
		A = T[w]
		phi = D[w]
		proyecciones.append([A * np.cos(2*pi*(w * m/n - phi)) for m in range(n)])
	return proyecciones


def sintesis(T, D, n):
	"""
	Nota que se necesita la dimensión n para recuperar a x (o, al menos, su paridad).
	"""
	proyecciones = formando_proyecciones(T, D)
	base_fourier = calculo_baseFourier(n)
	M = math.ceil((n-1)/2)

	#inicializamos la síntesis con el primer sumando
	x_sintesis = np.dot(proyecciones[0], base_fourier[0]) * base_fourier[0]
	for w in range(1, M):
		cos_w = np.dot(proyecciones[w], base_fourier[2*w-1]) * base_fourier[2*w-1]
		sen_w = np.dot(proyecciones[w], base_fourier[2*w]) * base_fourier[2*w]
		x_sintesis += cos_w + sen_w

	if n % 2 ==1: #n impar
		cos_w = np.dot(proyecciones[M], base_fourier[2*M-1]) * base_fourier[2*M-1]
		sen_w = np.dot(proyecciones[M], base_fourier[2*M]) * base_fourier[2*M]
		x_sintesis += cos_w + sen_w

	else: #n par
		cos_w = np.dot(proyecciones[M], base_fourier[2*M-1]) * base_fourier[2*M-1]
		x_sintesis += cos_w

	return x_sintesis








if __name__ == "__main__":
	print(calculo_baseFourier(5))

	x = [12, 5, 8, 9, 7, 4, -1.3, 12]
	n = len(x)
	M = math.ceil((n-1)/2)
	T, D = espectro_y_desfases_fourier(x, n, M)
	sintesis = sintesis(T, D, n)
	print(sintesis)
	#x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
	#grafica_espectro_fourier(x)
