#script para realizar los dos análisis espectrales finales.

#librerías de python

import numpy as np
from numpy.linalg import norm
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import pylab
from tqdm import tqdm
import random #para ejemplos

import pickle   #for python data persistency
import pandas as pd

#módulos personales
import base_legendreDiscreta as legendre
import proyecciones as proy #aquí tengo una función para hacer regresiones lineales
import funciones_figuras3d

pi=math.pi
mpl.rcParams.update(mpl.rcParamsDefault)
colores=['#fe01b1', 'gray', '#feb308', '#5f34e7', '#feb308', '#8f99fb', 'gray', '#8e82fe']

"""
El orden de los colores es:
    0: color de la señal a analizar
    1: color para los puntos del espectro que no son los de FP (frecuencia máxima) YA NO
    2: Frecuencia máxima del estudio basado en espacios monofrecuenciales
    3: Frecuencia máxima del estudio basado en la DFT
    4: Recta de mínimos cuadrados del estudio global de las FP encontradas con espacios monofrecuenciales
    5: Recta de mínimos cuadrados del estudio global de las FP encontradas con la DFT
    6: Color para los puntos del espectro monofrecuencial que no son máximos
    7: Color para los puntos del espectro de la DFT que no son máximos
"""


def formato_axis_derecha(axis):
  """
  Agregando los elementos que me gustan a un axis
  """
  axis.axhline(y=0, color='gray')
  axis.axvline(x=0, color='gray')
  axis.grid(True)
  axis.legend()
  #axis.legend(loc = 'lower right')

def formato_axis_izquierda(axis):
  """
  Agregando los elementos que me gustan a un axis
  """
  axis.axhline(y=0, color='gray')
  axis.axvline(x=0, color='gray')
  axis.grid(True)
  axis.legend(loc = 'upper left')
  
#def formato_axis_legendaAfuera_2col(axis):
#  """
#  Agregando los elementos que me gustan a un axis
#  """
#  axis.axhline(y=0, color='gray')
#  axis.axvline(x=0, color='gray')
#  axis.grid(True)
#
#  #To shrink the axis, but I don't like it.
#  #box = axis.get_position()
#  #axis.set_position([box.x0, box.y0 + box.height *0.1,
#  #                   box.width, box.height*0.9])
#  #axis.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.15), fancybox = True)
#  axis.legend(bbox_to_anchor = (1.1, 1.05), loc = 'center left')


#  ---------------------------------------- TDF ----------------------------------------

#Funciones coseno y coseno a partir de las que se construye todo lo que sigue.
def c_w(n,t, w):
  return math.sqrt(1/n)*np.cos(2*pi*w*t)

def s_w(n,t, w):
  return math.sqrt(1/n)*np.sin(2*pi*w*t)
  
def calculo_base(n): 
  """
  función que calcula la BON de Fourier (versión real) de dimensión n.
  """
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
    f_w = [ c_w(n, t, M) for t in dominio ]
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
    
def coeficientes_tau(x):
  """
  Estos coeficientes juntan la información de cosenos y senos de la misma
  frecuencia en un solo coeficiente.
  """
  n = len(x)
  M = math.ceil(n/2)
  coef_cosenos, coef_senos = coeficientes_espectrales(x) 
  norma = np.linalg.norm(x)
  if norma == 0:
      return None
  if norma == 1:
    taus = [abs(coef_cosenos[0])] #inicializamos la lista de sigmas con la primera entrada
    coef_cosenos.pop(0) #quitamos el coseno que ya agregamos a la lista de taus.
    
    if n%2==1:
      cuadrados_cosenos = np.square(coef_cosenos)
      cuadrados_senos = np.square(coef_senos)
      for i in range(M-1):
        taus.append(math.sqrt(cuadrados_cosenos[i]+cuadrados_senos[i]))
      return taus
    
    else:
      sigma_final = abs(coef_cosenos[M-1]) #guardamos el último sigma
      coef_cosenos.pop(M-1) #lo quitamos de la lista
      cuadrados_cosenos = np.square(coef_cosenos)
      cuadrados_senos = np.square(coef_senos)
      for i in range(M-1):
        taus.append(math.sqrt(cuadrados_cosenos[i]+cuadrados_senos[i]))
      taus.append(sigma_final)
      
      return taus

  else:  
    taus = [abs(coef_cosenos[0])/norma] #inicializamos la lista de sigmas con la primera entrada
    coef_cosenos.pop(0) #quitamos el coseno que ya agregamos a la lista de taus.
    
    if n%2==1:
      cuadrados_cosenos = np.square(coef_cosenos)
      cuadrados_senos = np.square(coef_senos)
      for i in range(M-1):
        taus.append(math.sqrt(cuadrados_cosenos[i]+cuadrados_senos[i])/norma)

      return taus
    
    else:
      sigma_final = abs(coef_cosenos[M-1])/norma #guardamos el último sigma
      coef_cosenos.pop(M-1) #lo quitamos de la lista
      cuadrados_cosenos = np.square(coef_cosenos)
      cuadrados_senos = np.square(coef_senos)
      for i in range(M-1):
        taus.append(math.sqrt(cuadrados_cosenos[i]+cuadrados_senos[i])/norma)
      taus.append(sigma_final)
      return taus


def grafica_taus_axis(x, n, nombre, axis1, axis2, legend_derecha = True):
  """
  'x' es un array de dimensión mayor a dos. 
  Esta función dibuja la gráfica de 'x' con dominio el tiempo, junto
  con la gráfica de los coeficientes tau de x.

  En 'axis1' se grafica la señal x,
  en 'axis2' se grafica el espectro.
  """
  M = math.ceil(n/2)

  dominio=[m/n for m in range(n)]
  taus = coeficientes_tau(x)
  cant_freq=len(taus)

  axis1.scatter(dominio, x, color= colores[0], s=50, label= "${0}$".format(nombre), zorder = 3)
  axis1.set_title('Gráfica de '+ "${0}$".format(nombre))

  for i in range(cant_freq):
    axis2.scatter(i, taus[i], color=colores[7], s=100, marker = '*', zorder = 2)

  axis1.set_ylabel('Coeficiente respecto a la base canónica')
  axis1.set_xlabel('Tiempo')

  axis2.set_xlabel('Frecuencias enteras $\omega$')
  axis2.set_ylabel(r'$\tau_{{{0}}}($'.format(str(n))+"${0}$".format(nombre)+ r'$, \omega)$' ) 

  X=np.arange(0, 1, 0.0001)

  coef_cosenos, coef_senos = coeficientes_espectrales(x) #TODO :( código con muchas partes repetidas !
  max_w = taus.index(max(taus)) #máxima frecuencia
  axis2.scatter(max_w, max(taus), color = colores[3], s = 70, label = '( ' + str(max_w) + ', ' + str(round(max(taus), 4))  + ' )', marker = '^', zorder = 3)

  if n %2 == 0: 
    if max_w == 0 or max_w == M:
      coef_cos = coef_cosenos[max_w] * math.sqrt(1/n)
      axis1.scatter(dominio, [coef_cos * np.cos(2*pi*max_w*t) for t in dominio], zorder=2, color = colores[3], s= 80)
      axis1.plot(X, coef_cos * np.cos(2*pi*max_w*X), color = colores[3], label = r'${{{0}}} \cdot cos(2 \pi \cdot {{{1}}} t) $'.format(str(round(coef_cos,2)), str(max_w)))
    else:
      coef_cos = coef_cosenos[max_w] * math.sqrt(2/n)
      coef_sen = coef_senos[max_w-1] * math.sqrt(2/n)
      muestreo = [coef_cos * np.cos(2*pi*max_w*t)+coef_sen* np.sin(2*pi*max_w*t) for t in dominio]
      axis1.scatter(dominio, muestreo, zorder = 2, color = colores[3], s=80) 
      axis1.plot(X, coef_cos * np.cos(2*pi*max_w*X)  + coef_sen * np.sin(2*pi*max_w*X), color = colores[3], label = r'${{{0}}} \cdot cos(2 \pi \cdot {{{1}}} t) + {{{2}}} \cdot sen(2 \pi \cdot {{{1}}} t) $'.format(str(round(coef_cos,2)), str(max_w), str(round(coef_sen,2))))
  else: #o sea, si n%2 == 1
    if max_w == 0:
      coef_cos = coef_cosenos[max_w] * math.sqrt(1/n)
      axis1.scatter(dominio, [coef_cos * np.cos(2*pi*max_w*t) for t in dominio], zorder = 2, color = colores[3], s=80)
      axis1.plot(X, coef_cos * np.cos(2*pi*max_w*X), color = colores[3], label = r'${{{0}}} \cdot cos(2 \pi \cdot {{{1}}} t) $'.format(str(round(coef_cos,2)), str(max_w))) #redundante
    else: 
      coef_cos = coef_cosenos[max_w] * math.sqrt(2/n)
      coef_sen = coef_senos[max_w-1] * math.sqrt(2/n)
      axis1.scatter(dominio, [coef_cos * np.cos(2*pi*max_w*t)+coef_sen* np.sin(2*pi*max_w*t) for t in dominio], zorder=2, color = colores[3], s=80)
      axis1.plot(X, coef_cos * np.cos(2*pi*max_w*X) + coef_sen * np.sin(2*pi*max_w*X), color = colores[3], label = r'${{{0}}} \cdot cos(2 \pi \cdot {{{1}}} t) + {{{2}}} \cdot sen(2 \pi \cdot {{{1}}} t) $'.format(str(round(coef_cos,2)), str(max_w), str(round(coef_sen,2))))
  
  if legend_derecha == True:
    formato_axis_derecha(axis1)
    formato_axis_derecha(axis2)
  else:
    formato_axis_izquierda(axis1)
    formato_axis_izquierda(axis2)



#  ------------------------------------- espacios monofrecuenciales ----------------------------------------

def vectores_frecuencia(n, w): #TODO quita a n del argumento.
  """
  n y w ambas de tipo int, n mayor o igual a dos, w no negativa.
  n es la dimensión, w la frecuencia 
  
  """
  c_nw = np.array([math.cos(2*pi*w*m/n) for m in range(n) ])

  if w % (n/2) == 0: # caso 2
    f_nw= 1/math.sqrt(n) * c_nw
    return f_nw
  
  else: #caso 1
    a= math.sin(2*pi*w) 
    b= math.cos(2*pi*w*(1-1/n))
    c= math.sin(2*pi*w/n)

    xi_nw=math.sqrt(2) * (n+a*b/c)**(-1/2)
    eta_nw=math.sqrt(2) * (n-a*b/c)**(-1/2)

    s_nw = np.array([math.sin(2*pi*w*m/n) for m in range(n) ])

    f_nw= xi_nw * c_nw
    g_nw= eta_nw * s_nw

  return (f_nw, g_nw, xi_nw, eta_nw)

#Funciones para el caso 1

def elementos_basicos_caso1(x, w):
  x = np.array(x) #Convirtiendo a 'x' (array) a np.array en caso de ser necesario.
  n = len(x)
  f_nw, g_nw, xi_nw, eta_nw = vectores_frecuencia(n, w)

  p = np.dot(f_nw, g_nw) #TODO no estoy usando la fórmula que calculé, no importa?
  q = np.dot(x, f_nw)
  r = np.dot(x, g_nw) 
  s = np.dot(x, x)

  return f_nw, g_nw, xi_nw, eta_nw, p, q, r, s


def sigma_caso1(x, w):
  """
  "x" es un np.array, w es un float mayor o igual a cero tal que
  n/2 (donde n= len(x)) no divide a w.
  """

  f_nw, g_nw, xi_nw, eta_nw, p, q, r, s = elementos_basicos_caso1(x, w)
  sigma = math.sqrt( (q**2 + r**2 -2*p*q*r)/ (s*(1-p**2)) )
  return sigma

def amplDesfase_caso1(x,w):
  f_nw, g_nw, xi_nw, eta_nw, p, q, r, s = elementos_basicos_caso1(x, w)
  den = 1 - p**2
  c = (q -p*r)*xi_nw/ den
  d = (r -p*q)*eta_nw/ den
  A = math.sqrt(c**2 + d**2) #Amplitud de la proyección
  
  tan = math.atan(d/c)
  if c>0 and d>0:
    phi = tan/(2*pi)
  elif (c<0 and d>0) or (c<0 and d<0): #tenía mal el último caso; esto era el error que buscabas:)
    phi = (tan + pi)/(2*pi)
  else:
    phi = (tan + 2*pi)/(2*pi)
  return A, phi

def grafica_sigma_amplDesfase_axis_caso1(x, w, axis, label_derecha = True):
  """
  Se dibujan la gráfica de x junto con la de su proyección al espacio de frecuencias P_{w}.

  """
  n = len(x)

  sigma = sigma_caso1(x, w)
  A, phi = amplDesfase_caso1(x, w)

  def coseno_amplDes(t):
    return A * np.cos(2*pi*w*t-2*pi*phi)

  dominio=[k/n for k in range(n)]
  proyeccion_Pnw = [coseno_amplDes(m/n) for m in range(n)]

  #axis.scatter(dominio, x, color=colores[0], s=80)
  axis.scatter(dominio, proyeccion_Pnw, color=colores[2], s=80, label = r'Gráfica de $\Pi_{{P_{{ {0}, {1} }} }}(x)$'.format(str(n),str(w)))
  
  X=np.arange(0, 1, 0.0001)
  axis.plot(X, coseno_amplDes(X), color=colores[2])
  if label_derecha == True:
    formato_axis_derecha(axis)
  else:
    formato_axis_izquierda(axis)

#Funciones para el caso 2

def elementos_basicos_caso2(x, w):
  x = np.array(x) #Convirtiendo a 'x' (array) a np.array en caso de ser necesario.
  n=len(x)
  f_nw = vectores_frecuencia(n, w)
  q = np.dot(x, f_nw)
  s = np.dot(x, x)

  return f_nw, q, s


def sigma_caso2(x, w):
  """
  "x" es un np.array, w es un float mayor o igual a cero tal que
  n/2 (donde n= len(x)) no divide a w.
  """

  f_nw, q, s = elementos_basicos_caso2(x, w)

  sigma = abs(q)/math.sqrt(s)
  return sigma

def amplDesfase_caso2(x,w):
  n = len(x)
  f_nw, q, s = elementos_basicos_caso2(x, w)
  A = q/math.sqrt(n)
  return A, 0

def grafica_sigma_amplDesfase_axis_caso2(x, w, axis, label_derecha = True):
  """
  Se dibujan la gráfica de x junto con la de su proyección al espacio de frecuencias P_{w}.

  """
  n = len(x)

  sigma = sigma_caso2(x, w)
  A, phi = amplDesfase_caso2(x, w)

  def coseno_amplDes(t):
    return A * np.cos(2*pi*w*t-2*pi*phi)

  dominio=[k/n for k in range(n)]
  proyeccion_Pnw = [coseno_amplDes(m/n) for m in range(n)]

  #axis.scatter(dominio, x, color=colores[0], s=80)
  axis.scatter(dominio, proyeccion_Pnw, color=colores[2], s=80, label=r'Gráfica de $\Pi_{{P_{{ {0}, {1} }} }}(x)$'.format(str(n),str(w)))
  
  X=np.arange(0, 1, 0.0001)
  axis.plot(X, coseno_amplDes(X), color=colores[2])
  if label_derecha == True:
    formato_axis_derecha(axis)
  else:
    formato_axis_izquierda(axis)


def analisis_espectral_espaciosMonofrecuenciales(x, n, frecuencias, nombre, axis0, axis1, legenda_derecha = True):
  """
  x es un array
  'frecuencias' es un vector de frecuencias.
  En 'axis0' se grafica la señal x, en el 'axis1' el espectro.
  """
  dominio_tiempo=[m/n for m in range(n)] 

  axis0.scatter(dominio_tiempo, x, color= colores[0], s= 50, label = "${0}$".format(nombre), zorder = 3)
  axis0.set_title('Gráfica de '+ "${0}$".format(nombre))

  sigmas = []
  for w in frecuencias: 
    if w % (n/2) == 0 :
      sigma = sigma_caso2(x, w)
      sigmas.append(sigma)
    else:
      sigma = sigma_caso1(x, w)
      print(sigma)
      sigmas.append(sigma)

  
  axis1.scatter(frecuencias, sigmas, color=colores[6], s=5, marker = '*', zorder = 1) #TODO descomenta.
  #axis1.scatter(frecuencias, sigmas, color= 'red', s=50, marker = '*', zorder = 5) #Útil para ver que el último punto de la TDF en efecto coincide con el del espectro de espacios monofrecuenciales.
  sigma_max = max(sigmas)
  frec_max = frecuencias[sigmas.index(sigma_max)]
    
  #axis1.scatter(frec_max, sigma_max, s = 70, color = colores[2], label = '$FP1({0})$'.format(nombre) + '='+str(frec_max), marker = 'v')
  axis1.scatter(frec_max, sigma_max, s = 100, color = colores[2], label = '( ' + str(frec_max) + ', ' + str(round(sigma_max, 2)) + ' )', marker = 'v')
  
  if frec_max % (n/2) == 0:
    grafica_sigma_amplDesfase_axis_caso2(x, frec_max, axis0)
  else: 
    grafica_sigma_amplDesfase_axis_caso1(x, frec_max, axis0)

  axis0.set_xlabel('Tiempo')
  axis0.set_ylabel('Coeficiente respecto a la base canónica')
  axis1.set_xlabel('Frecuencias ' + r'$\omega$')
  axis1.set_ylabel(r'$\sigma_{{{0}}}($'.format(str(n))+"${0}$".format(nombre)+ r'$, \omega)$' + r', $\tau_{{{0}}}($'.format(str(n))+"${0}$".format(nombre)+ r'$, \omega)$') 

  if legenda_derecha == True:
    formato_axis_derecha(axis0)
    formato_axis_derecha(axis1)
  else:
    formato_axis_izquierda(axis0)
    formato_axis_izquierda(axis1)

 
# ---------------------- Análisis espectrales para cualquier señal finita  -----------------------------------

def analisis_espectrales_mostrarGrafica(x, frecuencias, nombre):
  """
  'x' (tipo array) tiene las mediciones.
  'frecuencias' (tipo array) es un array de frecuencias.
  'nombre' (tipo string) es el nombre de la señal.
  """
  n = len(x) 
  fig = plt.figure()
  gs = fig.add_gridspec(2,2)
  
  ax0 = fig.add_subplot(gs[1, :1])
  ax1 = fig.add_subplot(gs[1, 1:])
  ax2 = fig.add_subplot(gs[0, :2])

  ax2.set_title('Espectros')
  ax2.axhline(y = 1, color = 'red', linewidth = 2, linestyle = 'dotted')

  grafica_taus_axis(x, n, nombre, ax0, ax2)
  analisis_espectral_espaciosMonofrecuenciales(x, n, frecuencias, nombre, ax1, ax2)

  fig.suptitle('Análisis espectrales de '+ nombre, fontsize = 18)
  fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
  return plt.show()

def analisis_espectrales_guardarGrafica(x, frecuencias, nombre):
  """
  'x' (tipo array) tiene las mediciones.
  'frecuencias' (tipo array) es un array de frecuencias.
  'nombre' (tipo string) es el nombre de la señal.
  Las proporciones se ven perfectas!
  """
  n = len(x) #TODO pasar como argumento a las otras funciones! 
  fig, axis = plt.subplots(2,2)
  fig.set_size_inches(11.25, 12.34) 
  #fig.tight_layout(pad=0.4)
  
  analisis_espectral_espaciosMonofrecuenciales(x,n, frecuencias, nombre, axis[1,0], axis[1,1])
  grafica_taus_axis(x,n, nombre, axis[0,0], axis[0,1])
  fig.suptitle('Análisis espectrales de ' + "${0}$".format(nombre) + r'$\in \mathbb{{R}}^{{{0}}}$'.format(str(n)), fontsize = 18)
  plt.savefig("/home/ame/GitHub/tesis-licenciatura/imagenes/estudios_espectrales/prueba.png")
  
  
# ------------------------------------- Análisis espectrales de los PDL ------------------------------------
def analisis_espectrales_PDL_mostrarGrafica(n, k, label_derecha = True):
  """
  'x' (tipo array) tiene las mediciones.
  'frecuencias' (tipo array) es un array de frecuencias.
  'nombre' (tipo string) es el nombre de la señal.
  """
  x = legendre.calculo_base(n)[k]
  frecuencias = [a/100 for a in range(int(n*100/2) + 1)]
  nombre = r'\mathcal{{L}}^{{{0}}}'.format(str(n)+','+str(k)) 

  fig = plt.figure()
  gs = fig.add_gridspec(2,2)
  #ax0 = fig.add_subplot(gs[0, :1])
  #ax1 = fig.add_subplot(gs[0, 1:])
  #ax2 = fig.add_subplot(gs[1, :2])
  
  ax0 = fig.add_subplot(gs[1, :1])
  ax1 = fig.add_subplot(gs[1, 1:])
  ax2 = fig.add_subplot(gs[0, :2])

  ax2.set_title('Espectros')
  ax2.axvline(x = k/2, color = 'red', linewidth = 2, label = "x = " + str(k/2))
  ax2.axhline(y = 1, color = 'red', linewidth = 2, linestyle = 'dotted')

  grafica_taus_axis(x, n, nombre, ax0, ax2, label_derecha)
  analisis_espectral_espaciosMonofrecuenciales(x, n, frecuencias, nombre, ax1, ax2, label_derecha)

  fig.suptitle('Análisis espectrales de ' + "${0}$".format(nombre) + r'$\in \mathbb{{R}}^{{{0}}}$'.format(str(n)), fontsize = 18)
  fig.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.5)
  #fig.tight_layout(top=0.914, bottom=0.058, left=0.041, right=0.992, hspace=0.202, wspace=0.09)
  return plt.show()

def analisis_espectrales_PDL_guardarGrafica(n,k):
   #TODO modificar! 
  """
  'x' (tipo array) tiene las mediciones.
  'frecuencias' (tipo array) es un array de frecuencias.
  'nombre' (tipo string) es el nombre de la señal.
  """
  #TODO poner a la ruta como argumento.
  x = legendre.calculo_base(n)[k]
  frecuencias = [a/100 for a in range(int(n*100/2) + 1)]
  nombre = '\mathcal{L}^{ {0} }'.format(str(n)+','+str(k)) 
  fig, axis = plt.subplots(2,2)
  fig.set_size_inches(11.25, 12.34) 

  grafica_taus_axis(x, n, nombre, axis[0,0], axis[0,1])
  analisis_espectral_espaciosMonofrecuenciales(x, n, frecuencias, nombre, axis[1,0], axis[1,1])
  fig.suptitle('Análisis espectrales de '+ "${0}$".format(nombre) +r'$\in \mathbb{{R}}^{{{0}}}$'.format(str(n)), fontsize = 18)
  #fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
  plt.savefig("/home/ame/GitHub/tesis-licenciatura/imagenes/estudios_espectrales/"+str(n)+'_'+str(k))



#---------------------------- Análisis espectral global ------------------------------------


def calculo_sigmaMax(x, n, frecuencias):
    """
    'x' es un array que representa las mediciones de la señal,
    'frecuencias' es un array con frecuencias respecto a la cuales comparar a 'x'.
    Esta función calcula todas las sigmas de 'x' respecto a las frecuencias contenidas en 'frecuencias'
    y regresa la sigma máxima (la que se resalta cuando se grafica el espectro.)
    """
    sigmas = []
    for w in frecuencias: 
      if w % (n/2) == 0 :
        sigma = sigma_caso2(x, w)
        sigmas.append(sigma)
      else:
        sigma = sigma_caso1(x, w)
        sigmas.append(sigma)
    
    sigma_max = max(sigmas) #buscamos la sigma mayor
    frec_max = frecuencias[sigmas.index(sigma_max)] #y la frecuencia asociada a esta.
    #regresamos tal frecuencia
    return frec_max


def calculo_tauMax(x):
    taus = coeficientes_tau(x)
    tau_max = max(taus)
    frec_max = taus.index(tau_max)
    return frec_max

def graficar_analisis_espectralPDL_global(n, label_derecha = True):
    fig, axis = plt.subplots(2,1)
    dominio_grados = [k for k in range(n)]
    X = np.arange(0,n,0.05)
#    
#   #extrayendo los datos de los .txt en los que se ha guardado
    with open('data_AE.txt', 'rb') as f:
        data_AE = pickle.load(f)

    datos_espectrales_dim_n = data_AE[n]
    sigmasMax_n, tausMax_n = datos_espectrales_dim_n[0], datos_espectrales_dim_n[1] 
    b0, m0 = datos_espectrales_dim_n[2], datos_espectrales_dim_n[3]
    b1, m1 = datos_espectrales_dim_n[4], datos_espectrales_dim_n[5]

    axis[1].scatter(dominio_grados, sigmasMax_n, marker = 'D', label = '$FP1$', color = colores[2], s = 150)
    axis[1].plot(X, X/2, label = r'$f(t)=\frac{t}{2}$', color = 'gray', linestyle = 'dotted')
    axis[1].plot(X, b1+X*m1, color = colores[4], label = r'RMC: $l(t) = {{{0}}}t + {{{1}}}$'.format(str(round(m1,2)), str(round(b1,2))))
    axis[0].scatter(dominio_grados, tausMax_n, label = '$FP0$', marker = 'D', color = colores[3] , s = 150)
    axis[0].plot(X, X/2, label = r'$f(t)=\frac{t}{2}$', color = 'gray', linestyle = 'dotted')
    axis[0].plot(X, b0+X*m0, color = colores[3], label = r'RMC: $l(t) = {{{0}}}t + {{{1}}}$'.format(str(round(m0,2)), str(round(b0,2))))
    plt.suptitle("Frecuencias máximas encontradas en los análisis espectrales \n de los PDL de dimensión "+str(n), fontsize = 14)

    if label_derecha == True:
        formato_axis_derecha(axis[0])
        formato_axis_derecha(axis[1])
    else:
        formato_axis_izquierda(axis[0])
        formato_axis_izquierda(axis[1])

    axis[0].set_xlabel('Grado $0 \leq k \leq {0}$'.format(str(n-1)))
    axis[1].set_xlabel('Grado $0 \leq k \leq {0}$'.format(str(n-1)))
    axis[0].set_ylabel(r'$FP0( \mathcal{{L}}^{{ {{{0}}}, k }} )$'.format(str(n)))
    axis[1].set_ylabel(r'$FP1( \mathcal{{L}}^{{ {{{0}}}, k }} )$'.format(str(n)))

    return plt.show()

def analisis_espectralPDL_global(n):
    base_legendre = legendre.calculo_base(n)
    dominio_grados = [k for k in range(n)]
    frecuencias = [a/100 for a in range(int(n*100/2) + 1)]

    sigmasMax_n, tausMax_n = [], [] 
    for k in range(n): #iteramos en los grados de los PDL de dimensión n
        vector_legendre = base_legendre[k]
        sigmasMax_n.append(calculo_sigmaMax(vector_legendre, n, frecuencias))
        tausMax_n.append(calculo_tauMax(vector_legendre))

    b1, m1 = proy.coef_RMC(dominio_grados, sigmasMax_n)
    b0, m0 = proy.coef_RMC(dominio_grados, tausMax_n)

    #Vamos a guardar los valores en el diccionario definido en el script 'diccionario_RMC.py'
    #los separadores no serán comas, sino dobles espacios.

    with open('data_AE.txt', 'rb') as f:
        diccionario = pickle.load(f)
    diccionario[n] = (sigmasMax_n, tausMax_n, b0, m0, b1, m1)

    with open('data_AE.txt', 'wb') as f:
        pickle.dump(diccionario, f)

    return sigmasMax_n, tausMax_n, b0, m0, b1, m1


def grafica_3d_n_k_FP(N):
    """
    N: int mayor o igual a 2 y menor o igual a 69, es la dimensión máxima cuyos datos
    se van a graficar.
    """
    fig = plt.figure()

    axis = fig.add_subplot(1, 1, 1, projection = '3d')
    funciones_figuras3d.dibuja_ejes_labelsPersonalizados(axis, 15, 'n', 'k', 'PR')
    axis.invert_zaxis()
    
    with open('data_AE.txt', 'rb') as f:
        data_AE = pickle.load(f)

    for n in range(2, N): #rango de las dimensiones para las que los datos se guardaron en data_AE.txt: 2-69
        sigmasMax_n = data_AE[n][0]
        tausMax_n = data_AE[n][1]

        for k in range(0, n-1): #iteramos en la variable de grado
            axis.scatter(n, k, tausMax_n[k], color = colores[3], s = 70, marker = '^')
            axis.scatter(n, k, sigmasMax_n[k], color = colores[2], s = 70, marker = 'v')
    
    #graficamos dos puntos de nuevo para poner etiquetas.
    sigmasMax_2 = data_AE[2][0]
    tausMax_2 = data_AE[2][1]
    axis.scatter(2, 0, tausMax_2[0], color = colores[3], s = 70, marker = '^', label = '$FP_{0}$')
    axis.scatter(2, 0, sigmasMax_2[0], color = colores[2], s = 70, marker = 'v', label = '$FP_{1}$')
    

    xx, yy = np.meshgrid(range(0, N), range(0, N))
    zz = 0*xx + 0.5*yy    
    axis.plot_surface(xx, yy, zz, color = 'gray', alpha = 0.6)
    axis.plot_wireframe(xx, yy, zz, color = 'white', alpha = 0.4)

    plt.legend()
    return plt.show()



def grafica_analisisGlobal_k_fijo(k, label_derecha = True):
    #TODO reescribir. Mejor guarda los valores en arrays para que puedas buscar máximos y mínimos y agustar
    #la longitud del y-eje.
    """
    k: int, con 0 \leq k \leq 68. 
    """
    fig, axis = plt.subplots(1,1)
    with open('data_AE.txt', 'rb') as f:
        data_AE = pickle.load(f)
    
    def graficando_puntos(n, graph_label = False):
        datos_dim_n = data_AE[n] #extraemos la información n-dimensional
        sigmasMax_n, tausMax_n = datos_dim_n[0], datos_dim_n[1]
        sigmaMax_n_k = sigmasMax_n[k]
        tauMax_n_k = tausMax_n[k]
        if graph_label == True:
            plt.scatter(n, tauMax_n_k, color = colores[3], marker = '^', label = r'$(n,FP0 (\mathcal{{L}}^{{n, {{{0}}} }}) )$'.format(str(k)))
            plt.scatter(n, sigmaMax_n_k, color = colores[2], marker = 'v', label = r'$(n,FP1 (\mathcal{{L}}^{{n, {{{0}}} }}) )$'.format(str(k)))
        else: 
            plt.scatter(n, tauMax_n_k, color = colores[3], marker = '^')
            plt.scatter(n, sigmaMax_n_k, color = colores[2], marker = 'v')
    
    #axis.set_ylim(k/2 - 1.5, k/2 + 1.5) #TODO no puedes recortar así. No todos se comportan bien.
    if k == 0 or k == 1:
        axis.axhline(y = k/2, color = 'red', label = r'$y = \frac{{ {{{0}}}  }}{{2}}$'.format(str(k)))
        graficando_puntos(3, True)
        for n in range(4, 70):
            graficando_puntos(n, False)
    else:
        for n in range(k+1, 69):
            graficando_puntos(n, False)
        graficando_puntos(69, True)
        axis.axhline(y = k/2, color = 'red', label = r'$y = \frac{{ {{{0}}}  }}{{2}}$'.format(str(k)))

    if label_derecha == True:
        formato_axis_derecha(axis)
    else:
        formato_axis_izquierda(axis)

    plt.suptitle('Gráficas de las frecuencias principales FP para los \n PDL de grado k = '+str(k))
    axis.set_xlim(k+0.5, 69.5)
    axis.set_xlabel('Dimensión $n$')
    axis.set_ylabel('FP del polinomio discreto de Legendre '+ r'$\mathcal{{L}}^{{n, {{{0}}} }} \in \mathbb{{R}}^{{n}}$'.format(str(k)))
    return plt.show()


def grafica_nube_b0m0_b1m1():
    fig, axis = plt.subplots(1,2)
    with open('data_AE.txt', 'rb') as f:
        data_AE = pickle.load(f)

    lista_m0, lista_m1 = [], [] #listas de pendientes
    lista_b0, lista_b1 = [], [] #listas de ordenadas al origen 

    for n in range(3, 70):
        data_AE_dim_n = data_AE[n] #extraemos la información de dimensión n
        b0_n, m0_n = data_AE_dim_n[2], data_AE_dim_n[3]
        #TODO cambiar el nombre de b0 a b0_n, etc
        b1_n, m1_n = data_AE_dim_n[4], data_AE_dim_n[5]
        
        axis[0].scatter(b0_n, m0_n, marker = 'x', color = colores[3])
        axis[1].scatter(b1_n, m1_n, marker = 'x', color = colores[2])
        
        lista_m0.append(m0_n)    
        lista_m1.append(m1_n)    
        lista_b0.append(b0_n)    
        lista_b1.append(b1_n)    
        
    m0_max = max(lista_m0)    
    axis[0].set_xlabel('$b_{0,n}$')
    axis[1].set_xlabel('$b_{1,n}$')

    axis[0].set_ylabel('$m_{0,n}$')
    axis[1].set_ylabel('$m_{1,n}$')
    
    for i in range(2):
        axis[i].axhline(y = 0.5, color = 'red')
        axis[i].set_ylim(0, 0.7)
        axis[i].set_xlim(-1.5, 0.5)

    plt.suptitle('Gráficas de los coeficientes de las rectas \n de mínimos cuadrados para toda $n$')
    
    for i in range(2):
        axis[i].set_xlim(-1.5, 0.5)
        axis[i].set_ylim(0.3, 0.65)
        formato_axis_derecha(axis[i])

    plt.show()


def tabla_informacion():
    #Vamos a extraer el diccionario del archivo binario:

    with open('data_AE.txt', 'rb') as f:
        data_AE = pickle.load(f)
    
    data_AE = pd.DataFrame.from_dict(data_AE).T
    data_AE.columns = ['sigmas', 'taus', 'b0', 'm0', 'b1', 'm1']
    data_AE = data_AE.drop(['sigmas', 'taus'], axis = 1)
    
    data_AE_html = data_AE.to_html() #tipo str, código de html de la tabla con la información
    #Guardamos el código html en el siguiente archivo:
    f = open('EEspectral_tabla_parametrosRectas.html', 'w')
    f.write(data_AE_html)
    f.close()

#------- ejemplos comparación de espectro-0 (TDF) con espectro-1 (espacios monofrecuenciales)--- 

def cos_con_ruido(t, A, w, phi):
    return A * math.cos(2*pi*w*t + 2*pi*phi) + np.random.uniform(0, 0.5)

def cos_sin_ruido(t, A, w, phi):
    return A * math.cos(2*pi*w*t + 2*pi*phi) 

def sinusoide_espectros(n, w, A, phi, nombre, ruido = True):
    """
    n: tamaño de la muestra, w: frecuencia, A: amplitud,
    phi: (entre 0 y 1) desfase, ruido == True sii se usa a la función
    cos_con_ruido para muestrear.
    """
    frecuencias = [a/100 for a in range(int(n*100/2) + 1)]
    #frecuencias = [w for w in range(math.floor(n/2)+1)] #TODO comenta cuando encuentres el error.
    if ruido == True : 
        x = [cos_con_ruido(m/n, A, w, phi) for m in range(n)]
    else:
        x = [cos_sin_ruido(m/n, A, w, phi) for m in range(n)]

    return analisis_espectrales_mostrarGrafica(x, frecuencias, nombre)
        

    

"""
Notas: 
    para n=100, analisis_espectralPDL_global(100) encontró coeficientes que eran tan pequeños
    que se redondeaban a cero, por lo que se tenían errores del tipo 'división por cero'.
"""

if __name__=='__main__':
  n, w, A, phi, nombre = 16, 3, 2.3, 0, 'x' #TODO checa este último tau...
  #sinusoide_espectros(n, w, A, phi, nombre, ruido = False)
  n, w, A, phi, nombre = 17, 3, 1, 0, 'x'
  #sinusoide_espectros(n, w, A, phi, nombre, ruido = False)
  
  n, w, A, phi, nombre = 36, 3.4, -1.5, 0.2, 'x' #este ejemplo está bien!
  #sinusoide_espectros(n, w, A, phi, nombre, ruido = True)

  n, w, A, phi, nombre = 36, 5, -1.5, 0.2, 'x' #este ejemplo está bien!
  #sinusoide_espectros(n, w, A, phi, nombre, ruido = False)


  n, w, A, phi, nombre = 32, 3, math.sqrt(2/32), 0, 'x' #este ejemplo está bien!
  
  #x = [cos_sin_ruido(m/n, A, w, phi) for m in range(n)]
  #vect_cos = [c_w(n, t/n, 16) for t in range(n)]
  #print(abs(np.dot(x, vect_cos))/np.linalg.norm(x))
  #sinusoide_espectros(n, w, A, phi, nombre, ruido = False)

  #TODO recuerda que también querías muestrear al sinusoide morado!

  #analisis_espectrales_guardarGrafica(x, frecuencias, nombre)

  #-------------Gráficas para la pregunta 1
  #analisis_espectrales_PDL_mostrarGrafica(7, 0) 
  #analisis_espectrales_PDL_mostrarGrafica(16, 12) 
  #for k in range(7):
  #  analisis_espectrales_PDL_mostrarGrafica(7,k) 


  #TODO quitamos a n=2 del análisis? De todas formas, no era una dimensión muy interesante.
  #analisis_espectrales_PDL_mostrarGrafica(2,0)
  #analisis_espectrales_PDL_mostrarGrafica(4,0)
  #analisis_espectrales_PDL_mostrarGrafica(3,1)
  #analisis_espectrales_PDL_mostrarGrafica(8,1)
  #analisis_espectrales_PDL_mostrarGrafica(34,1)
  #analisis_espectrales_PDL_mostrarGrafica(39, 15)
  #analisis_espectrales_PDL_mostrarGrafica(39, 8)
  #analisis_espectrales_PDL_mostrarGrafica(39, 5)
  analisis_espectrales_PDL_mostrarGrafica(39, 38)

  #analisis_espectrales_PDL_mostrarGrafica(84, 12)
  #analisis_espectrales_PDL_mostrarGrafica(84, 1)
  #analisis_espectrales_PDL_mostrarGrafica(84, 38)
  #analisis_espectrales_PDL_mostrarGrafica(85, 19) #Muy buen ejemplo, ambos sinusoides se ajustan con desfase.
  #analisis_espectrales_PDL_mostrarGrafica(32,5) 
  #analisis_espectrales_PDL_mostrarGrafica(33,5) 
  #analisis_espectrales_PDL_guardarGrafica(18,15) 


  #-------------Gráficas para la pregunta 2
  #grafica_analisisGlobal_k_fijo(0) 
  #grafica_analisisGlobal_k_fijo(1) 

  #grafica_analisisGlobal_k_fijo(7) 
  #grafica_analisisGlobal_k_fijo(20) 
  grafica_analisisGlobal_k_fijo(44) 
  #grafica_3d_n_k_FP(10)

  #-------------Gráficas para la pregunta 3
  #graficar_analisis_espectralPDL_global(30, False) 
  #graficar_analisis_espectralPDL_global(65, False) 
  #graficar_analisis_espectralPDL_global(58, False) 
  #graficar_analisis_espectralPDL_global(7, False) 
  #grafica_nube_b0m0_b1m1()


  # ----------------------------------------------------------
  #                scripts que ejecuto una sola vez
  # ----------------------------------------------------------

  #tabla_informacion() #Sólo debía ejecutar una vez esta función para obtener el código html de la tabla.
  #Como era código plano, yo escribí un .css sencillo y le agregué algunos elementos más al archivo
  #.html.

  #------------- Guardando datos con pickle
  #Ya ejecuté este loop y guardé los datos en data_AE.txt. No lo ejecutes más!
  #for n in range(2, 70):
  #    analisis_espectralPDL_global(n)

  #------- figuras para ejemplificar el rango de frecuencias usado.
  x = [-0.8, -2.5, -9.23, -0.7, -9.36, -29, -9]
  frecuencias = [a/100 for a in range(int(700/2))]
  nombre = 'x'

  #analisis_espectrales_mostrarGrafica(x, frecuencias, nombre)

  frecuencias = [a/100 for a in range(int(1400/2))]
  #analisis_espectrales_mostrarGrafica(x, frecuencias, nombre)
    
