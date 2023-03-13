#Nueva versión del script 'base_fourier_V0'


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import pylab

pi=math.pi


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



#Funciones coseno y coseno a partir de las que se construye todo lo que sigue.

def c_w(n,t, w):
  return math.sqrt(1/n)*np.cos(2*pi*w*t)

def s_w(n,t, w):
  return math.sqrt(1/n)*np.sin(2*pi*w*t)


def calculo_base(n): 
  dominio=[k/n for k in range(n)]
  M=math.ceil(n/2) #cota superior de las frecuencias consideradas en la base
  

  base_F=[(1/math.sqrt(n))*np.ones([n])] #inicializamos la base, que será un array. Ya incluimos la primera entrada.

  for w in range(1,M): #Nota que, a diferencia del caso complejo, el rango de frecuencia no tiene a 'n-1' como cota superior!
    f_w=[]
    g_w=[]
    for t in dominio:
      f_w.append(math.sqrt(2)*c_w(n, t, w))
      g_w.append(math.sqrt(2)*s_w(n, t, w))
    base_F.append(f_w)
    base_F.append(g_w)

  if n%2==1: #si N es impar, ya terminamos
    return base_F #Debemos multiplicar por \sqrt{2} para obtener elementos de Rn de norma uno
  else: #en caso contrario, falta agregar un vector con una frecuencia más alta
    f_w=[]
    for i in range(M):
      f_w.append(1/math.sqrt(n))
      f_w.append(-1/math.sqrt(n))
    base_F.append(f_w) #Nota que aquí no multiplicamos por \sqrt{2}
    return base_F



def coeficientes_espectrales(x):
  """
  'x' es un array de dimensión (digamos, 'n') mayor a dos
  Se regresan los coeficientes de x respecto a la BON Fn de Rn separados en dos listas: 
  una correspondiente a los vectores cosenos, y otra a los vectores senos.
  """
  n=len(x)
  M=math.ceil(n/2)
  base_frecuencias=calculo_base(n)

  coef_cosenos=[np.dot(x, base_frecuencias[0])] #agregamos el primer coeficiente, que corresponde a un coseno.
  coef_senos=[]

  for i in range(M-1):
    coef_senos.append(np.dot(x, base_frecuencias[2*i+2]))
    coef_cosenos.append(np.dot(x, base_frecuencias[2*i+1]))

  if n%2==0:
    coef_cosenos.append(np.dot(x, base_frecuencias[n-1]))

  return (coef_cosenos, coef_senos)


def coeficientes_sigma(x):
  """
  En construcción. Estos coeficientes juntan la información de cosenos y senos de la misma
  frecuencia en un solo coeficiente.
  """
  n=len(x)
  M=math.ceil(n/2)

  coef_cosenos, coef_senos=coeficientes_espectrales(x)
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


