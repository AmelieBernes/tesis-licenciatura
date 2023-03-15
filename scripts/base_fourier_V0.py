#Nueva versión del script 'base_fourier_V0'


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import pylab

pi=math.pi

#Funciones coseno y seno a partir de las que se construye todo lo que sigue.

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

