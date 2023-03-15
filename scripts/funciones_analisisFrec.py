"""
Funciones usadas para hacer análisis de frecuencia. 
Se supone que en otros scipts se definen bases/sistemas de frecuencia
via una función llamada 'calculo_base'



Aquí ejecuta todos los análisis, no en los scripts en los que
defines las bases!

"""
	
#------------------------------------------------------------------------------------


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import pylab
from tqdm import tqdm
from numpy.linalg import norm

import base_fourier_V0
import base_fourier_V1

pi=math.pi
mpl.rcParams.update(mpl.rcParamsDefault)
colores=['red', 'blue', 'darkviolet', 'gray', 'hotpink']


#Funciones auxiliares

def escXarr(a,x):
  """
  'x' es un array, 'n' es de tipo int o float.
  Se regresa un array cuya i-ésima entrada es n*x[i]
  """
  resultado=[] #inicializamos el array
  for entrada in x: #iteramos en el array x
    resultado.append(a*entrada)
  return resultado

def sumando_arrays(x,y):
  """
  'x' y 'y' son ambos de tipo array de la misma longitud.
  Esta función regresa un array cuya i-ésima entrada es X[i]+y[i].
  """
  resultado=[]
  for i in range(len(x)):
    resultado.append(x[i]+y[i])
  return resultado
  
  
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
	
	
	
#------------------------------------------------------------------------------------

def coeficientes_espectrales(x, script_sistFreq=base_fourier_V0):
  """
  'x' es un array de dimensión (digamos, 'n') mayor a dos
  'script_sistFreq' es un módulo de Python en el que se ha definido una función llamada
  'calculo_base' que calcula una base o un sistema de frecuencias.
  
  Se regresan los coeficientes de x respecto al sistema de frecuencias (i.e.
  sus productos puntos con esta) separados en dos listas: 
  una correspondiente a los vectores cosenos, y otra a los vectores senos.
  """
  n=len(x)
  M=math.ceil(n/2)
  base_frecuencias=script_sistFreq.calculo_base(n)

  coef_cosenos=[np.dot(x, base_frecuencias[0])] #agregamos el primer coeficiente, que corresponde a un coseno.
  coef_senos=[]

  for i in range(M-1):
    coef_senos.append(np.dot(x, base_frecuencias[2*i+2]))
    coef_cosenos.append(np.dot(x, base_frecuencias[2*i+1]))

  if n%2==0:
    coef_cosenos.append(np.dot(x, base_frecuencias[n-1]))

  return (coef_cosenos, coef_senos)



def reconstruirS_coefFrec(x, script_sistFreq=base_fourier_V0):
  """
  No funciona.
  Sólo tiene sentido cuando el sistema de frecuencias de hecho es una BON de Rn.
  
  Reconstruyendo a x (un array de dimensión n mayor a uno) a partir de sus coeficientes de frecuencia
  """
  n=len(x)
  M=math.ceil(n/2)

  base_frecuencias=calculo_base(n)
  coef_cosenos, coef_senos= coeficientes_espectrales(x, script_sistFreq)

  sintesis=escXarr(coef_cosenos[0],base_frecuencias[0]) #inicializando la síntesis
  
  #sumamos las contribuciones de los senos
  for i in range(M-1):
    sintesis=sumando_arrays( sintesis,escXarr(coef_senos[i],base_frecuencias[2*i+2]) )

  for i in range(1, M-1):
    sintesis=sumando_arrays( sintesis,escXarr(coef_cosenos[i],base_frecuencias[2*i+1]) )

  if n%2==0:
    sintesis=sumando_arrays( sintesis,escXarr(coef_cosenos[M],base_frecuencias[n-1]) )

  return sintesis




def grafica_dominio_frecuencia(x, script_sistFreq=base_fourier_V0):
  """
  'x' es un array de dimensión mayor a dos. 
  Esta función dibuja la gráfica de 'x' (como se definió en ??), que es una cuyo dominio es el tiempo, junto
  con la gráfica de los coeficientes de 'x' respecto a la BON de frecuencias Fn, que puede pensarse como una gráfica
  de 'x' cuyo dominio es la frecuencia.

  """
  fig, axis= plt.subplots(1,2)
  fig.set_size_inches(13, 5.5)
  n=len(x)
  M=math.ceil(n/2)
  dominio_tiempo=[t for t in range(n)]
  coef_cosenos, coef_senos=coeficientes_espectrales(x, script_sistFreq)

  axis[0].scatter(dominio_tiempo, x, color=colores[4], s=150, label='Gráfica de $x$')
  axis[0].set_title('Dominio: tiempo')


  axis[1].set_title('Dominio: frecuencia')

  for i in range(M): #iteramos en las frecuencias
    axis[1].scatter(i, coef_cosenos[i], s=150, color='red')
    
  for i in range(M-1): #iteramos en las frecuencias
    axis[1].scatter(i+1, coef_senos[i], s=150, color='blue')

  #Dibujando etiquetas
  axis[1].scatter(0, coef_cosenos[0], s=150, color='red', label='Amplitudes asociadas a cosenos')
  axis[1].scatter(1, coef_senos[0], s=150, color='blue', label='Amplitudes asociadas a senos')

  if n%2==0:
    axis[1].scatter(M, coef_cosenos[M], s=150, color='red')

  axis[0].set_xlabel('Tiempo')
  axis[0].set_ylabel('Coeficiente respecto a la base canónica')
  axis[1].set_xlabel('Frecuencias enteras')
  axis[1].set_ylabel('Amplitudes')

  for i in range(2):
    axis[i].axhline(y=0, color='gray')
    axis[i].axvline(x=0, color='gray')
    axis[i].grid(True)
    axis[i].legend()

  return plt.show()
  
  
 
  
def coeficientes_sigma(x, script_sistFreq=base_fourier_V0):
  """
  En construcción. Estos coeficientes juntan la información de cosenos y senos de la misma
  frecuencia en un solo coeficiente.
  """
  n=len(x)
  M=math.ceil(n/2)

  coef_cosenos, coef_senos=coeficientes_espectrales(x, script_sistFreq)
  sigmas=[coef_cosenos[0]**2] #inicializamos la lista de sigmas con la primera entrada
  coef_cosenos.pop(0)
  
  if n%2==1:
    cuadrados_cosenos=np.square(coef_cosenos)
    cuadrados_senos=np.square(coef_senos)
    for i in range(M-1):
      sigmas.append(cuadrados_cosenos[i]+cuadrados_senos[i])
    return sigmas
  
  else:
    sigma_final=coef_cosenos[M-1]**2 #guardamos el último sigma
    coef_cosenos.pop(M-1)
    cuadrados_cosenos=np.square(coef_cosenos)
    cuadrados_senos=np.square(coef_senos)
    for i in range(M-1):
      sigmas.append(cuadrados_cosenos[i]+cuadrados_senos[i])
    sigmas.append(sigma_final)
    return sigmas
    
    
    
    
def grafica_sigmas(x, script_sistFreq=base_fourier_V0):
  """
  'x' es un array de dimensión mayor a dos. 
  Esta función dibuja la gráfica de 'x' (como se definió en ??), que es una cuyo dominio es el tiempo, junto
  con la gráfica de los coeficientes sigma de x.

  """
  fig, axis= plt.subplots(1,2)
  fig.set_size_inches(13, 5.5)
  n=len(x)

  dominio_tiempo=[t for t in range(n)]
  sigmas=coeficientes_sigma(x, script_sistFreq)
  cant_freq=len(sigmas)


  axis[0].scatter(dominio_tiempo, x, color=colores[4], s=150, label='Gráfica de $x$')
  axis[0].set_title('Gráfica de $x$')

  axis[1].set_title('Gráfica de los coeficientes sigma de x')
  for i in range(cant_freq):
    axis[1].scatter(i, sigmas[i], color='purple', s=150)


  axis[0].set_xlabel('Tiempo')
  axis[0].set_ylabel('Coeficiente respecto a la base canónica')
  axis[1].set_xlabel('Frecuencias enteras')
  axis[1].set_ylabel('Coef. sigma')

  for i in range(2):
    axis[i].axhline(y=0, color='gray')
    axis[i].axvline(x=0, color='gray')
    axis[i].grid(True)
    axis[i].legend()

  return plt.show()
  

if __name__ == "__main__":
	x=[1,-3.5,4,0,4,0,12,20,8]
	grafica_dominio_frecuencia(x, base_fourier_V1)
	#grafica_sigmas(x, base_fourier_V1)
	

