#script para realizar los dos análisis espectrales finales.

#TODO (14 de mayo, 2023) vamos a limpiar el código para aijar una banda de frecuencias y no checar siempre la 
#pertenecia de la frecuencia a n/2\IZ, pues eso sólo pasará para los extremos del espectro que fijemos. 

"""
#  ---------------------------------------- -- ----------------------------------------

				IMPORTANDO LIBERÍAS Y DEFINIENDO
			      	 FUNCIONES Y OBJETOS AUXILIARES

#  ---------------------------------------- -- ----------------------------------------
"""

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

ruta_carpeta = '/home/ame/GitHub/tesis-licenciatura/imagenes/graficas_analisisEspectrales/'

#colores=['#fe01b1', 'gray', '#00F7CE', '#5f34e7', '#feb308', '#8f99fb', 'gray', '#8e82fe', '#CE00D8']
#colores=['#fe01b1', 'gray', '#0ae408', '#5f34e7', '#feb308', '#8f99fb', 'gray', '#8e82fe','#CE00D8' ]

colores=['#fe01b1', '#3b2747', '#feb308', '#6832e3', '#feb308', '#8f99fb', 'gray', '#8e82fe', '#fd183d']

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
    8: Color para lineas de referencia
"""

#Este formato de axis lo uso para los espectros
def formato_axis(axis):
  """
  Agregando los elementos que me gustan a un axis
  """
  axis.axhline(y=0, color='gray', zorder = 0)
  axis.axvline(x=0, color='gray', zorder = 0)
  axis.grid(True)
  axis.legend(loc = 'best')
  #axis.legend(loc = 'lower right')

#Este formato de axis lo uso para las gráficas de los sinusoides
def formato_axis_leyenda_exterior(axis):
  """
  Agregando los elementos que me gustan a un axis
  """
  axis.axhline(y=0, color='gray', zorder = 0)
  axis.axvline(x=0, color='gray', zorder = 0)
  axis.grid(True)
  axis.legend(loc = 'best', bbox_to_anchor = (0.7, 1.35))
  #axis.legend(loc = 'upper left', frameon = False)
  
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


"""
#  ---------------------------------------- -- ----------------------------------------
                    N U E V O 

                    CÁLCULOS PARA GRAFICAR
                    EL LÍMITE POR CERO DEL ESPECTRO

X0 (primer momento): sum(x)
X1 (segundo momento): \suma{m=0}{n-1}{m*x[m]}

Por las pruebas numéricas que he hecho, parece que sí funciona.
#  ---------------------------------------- -- ----------------------------------------
"""

def limite_cero(x, n):
    """
    Función que da el límite por la derecha del cero del 
    espectro de la señal x de dimensión n.
    """
    pi = np.pi
    def momento_k(x, k):
        suma = 0
        for m in range(n):
            suma += m**k  * x[m]
        return suma

    norma = np.linalg.norm(x)
    X0 = momento_k(x, 0)
    X1 = momento_k(x, 1)
    
    num = 2*X0**2*(2*n-1)*(n-1) +12*X1**2-12*X0*X1*(n-1)
    den = norma**2 * n * (n-1)*(n+1)

    return math.sqrt(num/den) 

#TODO optimiza.
def alternado(x):
    """
    Dado un array 'x', la función regresa el vector alternado de 'x'.
    """
    Ax = []
    n = len(x)
    for m in range(n):
        Ax.append((-1)**m * x[m])
    return Ax

def limite_n_medios(x,n):
    """
    Función que da el límite por la izquierda de n/2 del 
    espectro de la señal x de dimensión n.
    """
    return limite_cero(alternado(x), n)

"""
#  ---------------------------------------- -- ----------------------------------------

				CÁLCULOS PARA REALIZAR EL PRIMER
			      ANÁLISIS ESPETRAL (BASADO EN LA TDF)

#  ---------------------------------------- -- ----------------------------------------
"""

#Funciones coseno y coseno a partir de las que se construye todo lo que sigue.
#En realidad, estos son casos particulares de (???), pero las programamos para evitar
#realizar cálculos que sabemos son iguales a cero.

def calculo_baseFourier(n): 
  """
  función que calcula la BON de Fourier (versión real) de dimensión n.
  """
  #TODO cambiar notacion de f_w a c_nw y de g_w a s_nw.
  dominio=[k/n for k in range(n)]
  M=math.ceil(n/2) #TODO poner como argumento.
  
  base_F=[math.sqrt(1/n)*np.ones([n])] #inicializamos la base, que será un array. Ya incluimos la primera entrada.

  for w in range(1,M): #Nota que, a diferencia del caso complejo, el rango de frecuencia no tiene a 'n-1' como cota superior!
    f_w=[math.sqrt(2/n)*np.cos(2*pi*w*t) for t in dominio]
    g_w=[math.sqrt(2/n)*np.sin(2*pi*w*t) for t in dominio]
    base_F.append(f_w)
    base_F.append(g_w)
    
  if n%2==1: #si n es impar, ya terminamos
    return base_F 
  else: #en caso contrario, falta agregar un vector coseno de frecuencia M
    f_w = [ math.sqrt(1/n)*np.cos(2*pi*M*t) for t in dominio ]
    base_F.append(f_w)
    return base_F
    
def coeficientes_base_Fourier(x): 
  """
  'x' es un array de dimensión (digamos, 'n') mayor a dos
  Se regresan los coeficientes de x respecto a la BON Fn de Rn separados en dos listas: 
  una correspondiente a los vectores cosenos, y otra a los vectores senos.
  """
  n=len(x) #TODO poner como argumento.
  M=math.ceil(n/2) #TODO poner como argumento.
  base_frecuencias = calculo_baseFourier(n)

  coef_cosenos = [np.dot(x, base_frecuencias[0])] #agregamos el primer coeficiente, que corresponde a un coseno.
  coef_senos = [] #inicializamos la lista de coef que son productos punto con un seno

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
  
  Como para calcularlos se usan los coeficietes_espectrales(x) y estos son usados
  también en otras partes justo después de usar coef
  """
  #TODO como también se usan los coef espectrales después, sería muy bueno devolver a estos como
  #output de esta función también, y no sólo los taus.
  
  #TODO: mucho mejor si primero calculas los coef_espectrales y LUEGO los mandas como
  #argumento de esta función. #TODO ya intenté hacer esto pero, no sé por qué, me rompe todo!!!
  
  n = len(x) #TODO poner como argumento.
  M = math.ceil(n/2) #TODO poner como argumento.
  coef_cosenos, coef_senos = coeficientes_base_Fourier(x) 
  norma = np.linalg.norm(x)
  if norma == 0:
      return None
  if norma == 1: #Distinguir este caso nos ahorra tener que realizar algunas divisiones innecesarias.
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


def grafica_taus_axis(x, n, nombre, axis1, axis2):
  """
  'x' es un array de dimensión mayor a dos. 
  Esta función dibuja la gráfica de 'x' con dominio el tiempo, junto
  con la gráfica de los coeficientes tau de x.

  En 'axis1' se grafica la señal x, #TODO ver si uso lo mismo para grafica sigmas.
  en 'axis2' se grafica el espectro.
  """
  M = math.ceil(n/2) #TODO poner como argumento
  coef_cosenos, coef_senos = coeficientes_base_Fourier(x)
  
  dominio=[m/n for m in range(n)]
  taus = coeficientes_tau(x)
  cant_freq=len(taus)

  axis1.scatter(dominio, x, color= colores[0], s=60, label= "${0}$".format(nombre), zorder = 3)
  axis1.set_title('Gráfica de '+ "${0}$".format(nombre))

  for i in range(cant_freq):
    axis2.scatter(i, taus[i], color=colores[7], s=100, marker = '*', zorder = 2)

  axis1.set_ylabel('Coeficiente respecto a la base canónica')
  axis1.set_xlabel('Tiempo')

  axis2.set_xlabel('Frecuencias enteras $\omega$')
  axis2.set_ylabel(r'$\tau_{{{0}}}($'.format(str(n))+"${0}$".format(nombre)+ r'$, \omega)$' ) 

  X=np.arange(0, 1, 0.0001)
  
  max_w = taus.index(max(taus)) #máxima frecuencia
  axis2.scatter(max_w, max(taus), color = colores[3], s = 70, label = '( ' + str(max_w) + ', ' + str(round(max(taus), 4))  + ' )', marker = '^', zorder = 3)
  stemlines = axis2.stem(max_w, max(taus), markerfmt = ' ', linefmt = '--')
  plt.setp(stemlines, color = colores[3]) #artista Python

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
  

  formato_axis_leyenda_exterior(axis1)
  formato_axis(axis2)


"""
#  ---------------------------------------- -- ----------------------------------------

				CÁLCULOS PARA REALIZAR EL SEGUNDO
			ANÁLISIS ESPETRAL (BASADO EN ESPACIOS MONOFRECUENCIALES)

#  ---------------------------------------- -- ----------------------------------------
"""

def vectores_frecuencia(n, w): #TODO quita a n del argumento.
  """
  n y w ambas de tipo int, n mayor o igual a dos, w no negativa.
  n es la dimensión, w la frecuencia 
  
  """
  C_nw = np.array([np.cos(2*pi*w*m/n) for m in range(n) ])

  if w % (n/2) == 0: # caso 2
    c_nw= 1/math.sqrt(n) * C_nw
    return c_nw
  
  else: #caso 1
    a= math.sin(2*pi*w) 
    b= math.cos(2*pi*w*(1-1/n))
    c= math.sin(2*pi*w/n)

    xi_nw=math.sqrt(2) * (n+a*b/c)**(-1/2)
    eta_nw=math.sqrt(2) * (n-a*b/c)**(-1/2)

    S_nw = np.array([np.sin(2*pi*w*m/n) for m in range(n) ])

    c_nw= xi_nw * C_nw
    s_nw= eta_nw * S_nw

  return c_nw, s_nw, xi_nw, eta_nw

#Funciones para el caso 1

def elementos_basicos_caso1(x, w):
  x = np.array(x) #Convirtiendo a 'x' (array) a np.array en caso de ser necesario.
  n = len(x)
  c_nw, s_nw, xi_nw, eta_nw = vectores_frecuencia(n, w)

  p = np.dot(c_nw, s_nw) #TODO no estoy usando la fórmula que calculé, no importa?
  q = np.dot(x, c_nw)
  r = np.dot(x, s_nw) 
  s = np.dot(x, x)

  return c_nw, s_nw, xi_nw, eta_nw, p, q, r, s


def sigma_caso1(x, w):
  """
  "x" es un np.array, w es un float mayor o igual a cero tal que
  n/2 (donde n= len(x)) no divide a w.
  """

  c_nw, s_nw, xi_nw, eta_nw, p, q, r, s = elementos_basicos_caso1(x, w)
  sigma = math.sqrt( (q**2 + r**2 -2*p*q*r)/ (s*(1-p**2)) )
  return sigma

def amplDesfase_caso1(x,w):
  c_nw, s_nw, xi_nw, eta_nw, p, q, r, s = elementos_basicos_caso1(x, w)
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


#TODO quita lo de 'caso 1' del nombre
def grafica_sigma_amplDesfase_axis_caso1(x, w, axis):
  """
  
  Se dibujan en 'axis' (que es el axis de alguna figura)
  la gráfica de x junto con la de su proyección al espacio de frecuencias P_{n,w}.
  donde n = len(x) y w>0 es la frecuencia dada como input.
  
  """
  #TODO cambiar nombre a grafica_amplDesfase_axis_caso1
  n = len(x)

  #sigma = sigma_caso1(x, w) #TODO por qué calculo la sigma? esto no es necesario.
  A, phi = amplDesfase_caso1(x, w)

  def coseno_amplDes(t):
    return A * np.cos(2*pi*w*t-2*pi*phi)

  dominio=[k/n for k in range(n)]
  proyeccion_Pnw = [coseno_amplDes(m/n) for m in range(n)]

  #axis.scatter(dominio, x, color=colores[0], s=80)

  #axis.scatter(dominio, proyeccion_Pnw, color='green', s=80, label = r'Gráfica de $\Pi_{{P_{{ {0}, {1} }} }}(x)$'.format(str(n), str(round(w, 2))) )
  
  X=np.arange(0, 1, 0.0001)
  axis.plot(X, coseno_amplDes(X), color=colores[2])

  formato_axis_leyenda_exterior(axis)

#Funciones para el caso 2
#TODO ver si ya puedo quitar todo lo del caso 2
#def elementos_basicos_caso2(x, w):
#  x = np.array(x) #Convirtiendo a 'x' (array) a np.array en caso de ser necesario.
#  n=len(x)
#  c_nw = vectores_frecuencia(n, w)
#  q = np.dot(x, c_nw)
#  s = np.dot(x, x)
#
#  return c_nw, q, s
#
#
#def sigma_caso2(x, w):
#  """
#  "x" es un np.array, w es un float mayor o igual a cero tal que
#  n/2 (donde n= len(x)) no divide a w.
#  """
#
#  c_nw, q, s = elementos_basicos_caso2(x, w)
#
#  sigma = abs(q)/math.sqrt(s)
#  return sigma
#
##TODO hay muchísimas repeticiones en todas estas funciones !!!
#def amplDesfase_caso2(x,w):
#  n = len(x)
#  c_nw, q, s = elementos_basicos_caso2(x, w)
#  A = q/math.sqrt(n)
#  return A, 0
#
#def grafica_sigma_amplDesfase_axis_caso2(x, w, axis):
#  """
#  Se dibujan la gráfica de x junto con la de su proyección al espacio de frecuencias P_{w}.
#
#  """
#  n = len(x)
#
#  sigma = sigma_caso2(x, w)
#  A, phi = amplDesfase_caso2(x, w)

  def coseno_amplDes(t):
    return A * np.cos(2*pi*w*t-2*pi*phi)

  dominio=[k/n for k in range(n)]
  proyeccion_Pnw = [coseno_amplDes(m/n) for m in range(n)]

  #axis.scatter(dominio, x, color=colores[0], s=80)
  axis.scatter(dominio, proyeccion_Pnw, color=colores[2], s=80, label=r'Gráfica de $\Pi_{{P_{{ {0}, {1} }} }}(x)$'.format(str(n),str( round(w, 2) ) ))
  
  X=np.arange(0, 1, 0.0001)
  axis.plot(X, coseno_amplDes(X), color=colores[2])

  formato_axis_leyenda_exterior(axis)


def analisis_espectral_espaciosMonofrecuenciales(x, n, frecuencias, nombre, axis0, axis1):
  """
  x es un array
  'frecuencias' es un vector de frecuencias.
  En 'axis0' se grafica la señal x, en el 'axis1' el espectro.
  """
  dominio_tiempo=[m/n for m in range(n)] 

  axis0.scatter(dominio_tiempo, x, color= colores[0], s= 50, label = "${0}$".format(nombre), zorder = 3)
  axis0.set_title('Gráfica de '+ "${0}$".format(nombre))

  #Calculando el espectro
  sigmas = []
  for w in frecuencias:
    sigmas.append(sigma_caso1(x,w))
  axis1.scatter(frecuencias, sigmas, color=colores[1], s=5, marker = '*', zorder = 2)

  sigmas.insert(0, limite_cero(x,n))
  axis1.scatter(0, sigmas[0],color=colores[1], s=5, marker = '*', zorder = 2)
  sigmas.append(limite_n_medios(x,n))
  axis1.scatter(n/2, sigmas[-1], color=colores[1], s=5, marker = '*', zorder = 2)
  #limite_0 = limite_cero(x, n) #NUEVO 
  #limite_n_2 = limite_n_medios(x, n) #NUEVO 
  sigma_max = max(sigmas)

  
  #Se tiene que agregar un ciclo if-else porque acorté el vector de frecuencias.
  # posiciones de sigmas =[0, 1, 2, ... , n*100/2]
  # posiciones de frecuencias = [1, 2, ..., n*100/2 -1]
  indice_sigmaMax = sigmas.index(sigma_max) #Mayor frecuencia en la que el espectro tiene su máximo
  if indice_sigmaMax == 0:
      frec_max = 0
  elif indice_sigmaMax == n*100/2:
      frec_max = n/2
  else:
    frec_max = frecuencias[indice_sigmaMax - 1] #TODO esto está bien ??

  axis1.scatter(frec_max, sigma_max, s = 100, color = colores[2], label = '( ' + str(round(frec_max, 2)) + ', ' + str(round(sigma_max, 2)) + ' )', marker = 'v', zorder=4)
  stemlines = axis1.stem(frec_max, sigma_max, markerfmt = ' ', linefmt = '--')
  plt.setp(stemlines, color = colores[2]) 
  
  if frec_max == 0:
    X = np.arange(0,1, 0.0001)
    proyeccion_W_n1 = proy.proyeccion(x,1)
    b0, b1 = proy.coef_RMC(dominio_tiempo, proyeccion_W_n1)
    axis0.scatter(dominio_tiempo, proyeccion_W_n1, color=colores[2], s=80, label = r'Gráfica de $\Pi_{{ W_{{ {0}, 1 }} }}(x)$'.format(str(n)) )
    axis0.plot(X, b0 + b1*X, color = colores[2])
    formato_axis_leyenda_exterior(axis0)
  elif frec_max == n/2:
    X = np.arange(0,1, 0.0001)
    x_alter = alternado(x)
    proyeccion_W_n1 = proy.proyeccion(x_alter,1)
    axis0.scatter(dominio_tiempo, proyeccion_W_n1, color=colores[2], s=80, label = r'Gráfica de $\Pi_{{ W_{{ {0}, 1 }} }}(A_{{ {0}  }} (x))$'.format(str(n)) )
    axis0.plot(dominio_tiempo, x_alter, color = colores[2])
    formato_axis_leyenda_exterior(axis0)
  else:    
    grafica_sigma_amplDesfase_axis_caso1(x, frec_max, axis0)

  axis0.set_xlabel('Tiempo')
  axis0.set_ylabel('Coeficiente respecto a la base canónica')
  axis1.set_xlabel('Frecuencias ' + r'$\omega$')
  axis1.set_ylabel(r'$\sigma_{{{0}}}($'.format(str(n))+"${0}$".format(nombre)+ r'$, \omega)$' + r', $\tau_{{{0}}}($'.format(str(n))+"${0}$".format(nombre)+ r'$, \omega)$') 


  formato_axis_leyenda_exterior(axis0)
  formato_axis(axis1)


"""
#  ---------------------------------------- -- ----------------------------------------

				   FUNCIONES A LLAMAR PARA
				REALIZAR ANÁLISES ESPECTRALES

#  ---------------------------------------- -- ----------------------------------------
"""

#TODO en el análisis de monof quita el if w in n/2 \Z, pues ya sabes que eso pasa solo en los extremos!

def grafica_analisisEspectrales(x, nombre, graficar = True):

  """
  Función para graficar o guardar (dependiendo del valor del booleano 'graficar')
  el análisis espectral de cualquier señal finita (no cero).
  'x' (tipo array) tiene las mediciones.
  'nombre' (tipo string) es el nombre de la señal.

  """
  n = len(x) 
  frecuencias = [a/100 for a in range(int(n*100/2) + 1)] 
  #Quitamos las frecuencias extremas, que son casos especiales.
  frecuencias.pop(0)
  frecuencias.pop(-1)

  fig = plt.figure()
  gs = fig.add_gridspec(2,2)
  
  ax0 = fig.add_subplot(gs[1, :1])
  ax1 = fig.add_subplot(gs[1, 1:])
  ax2 = fig.add_subplot(gs[0, :2])

  ax2.set_title('Espectros')
  ax2.axhline(y = 1, color = colores[8], linewidth = 2, linestyle = 'dotted')

  grafica_taus_axis(x, n, nombre, ax0, ax2)
  analisis_espectral_espaciosMonofrecuenciales(x, n, frecuencias, nombre, ax1, ax2)

  fig.suptitle('Análisis espectrales de '+ '${0}$'.format(nombre), fontsize = 18)
  #fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
  if graficar == True:
    fig.tight_layout()
    return plt.show()
  else:
  #TODO poner a la ruta en la que se guardará la imagen como argumento de la función.
    fig.set_size_inches(11.25, 12.34) 
    fig.tight_layout()
    plt.savefig(ruta_carpeta + str(n) + '_' + nombre)
  
def grafica_analisisEspectrales_PDL(n,k, graficar = True):
  """
  
  Función para graficar el análisis espectral de un polinomio discreto de Legendre (PDL).
  
  'n' (entero mayor o igual a dos) es la dimensión del PDL, 
  'k' \in \{0, 1, ... , n-1 \} es su grado.
  
  """
  #TODO nombre antiguo: analisis_espectrales_PDL_guardarGrafica(n,k)
  #TODO poner a la ruta como argumento. Como es siempre la misma ruta, ponla como una variable global,
  #para no tener que pedirla siempre como argumento.
  x = legendre.calculo_base(n)[k]
  frecuencias = [a/100 for a in range(int(n*100/2) + 1)] 
  #Quitamos las frecuencias extremas, que son casos especiales.
  frecuencias.pop(0)
  frecuencias.pop(-1)
  nombre = r'\mathcal{{L}}^{{{0}}}'.format(str(n)+','+str(k)) 

  fig = plt.figure()
  gs = fig.add_gridspec(2,2)
  ax0 = fig.add_subplot(gs[1, :1])
  ax1 = fig.add_subplot(gs[1, 1:])
  ax2 = fig.add_subplot(gs[0, :2])

  ax2.set_title('Espectros')
  ax2.axvline(x = k/2, color = colores[8], linewidth = 2, label = "x = " + str(k/2))
  ax2.axhline(y = 1, color = colores[8], linewidth = 2, linestyle = 'dotted')

  grafica_taus_axis(x, n, nombre, ax0, ax2)
  analisis_espectral_espaciosMonofrecuenciales(x, n, frecuencias, nombre, ax1, ax2)

  fig.suptitle('Análisis espectrales de ' + "${0}$".format(nombre) + r'$\in \mathbb{{R}}^{{{0}}}$'.format(str(n)), fontsize = 18)

  #fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
  if graficar == True: 
    fig.tight_layout()
    return plt.show()
  else:
    #fig.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.5)
    fig.set_size_inches(11.25, 12.34) 
    fig.tight_layout() #NOTA esto tiene que ir sólo una vez y justo antes de guardar la figura
                       #para no deformar nada!
    plt.savefig(ruta_carpeta + str(n) + '_' + str(k))



"""
#  ---------------------------------------- -- ----------------------------------------

				  EXTRAYENDO Y GUARDANDO INFORMACIÓN

#  ---------------------------------------- -- ----------------------------------------
"""

def calculo_sigmaMax(x, n, frecuencias):
    """
    'x' es un array que representa las mediciones de la señal,
    'frecuencias' es un array con frecuencias respecto a la cuales comparar a 'x'.
    Esta función calcula todas las sigmas de 'x' respecto a las frecuencias contenidas en 'frecuencias'
    y regresa la sigma máxima (la que se resalta cuando se grafica el espectro.)
    #TODO en realidad, la descripción es incorrecta.
    """
    #TODO cambiar el vector de frecuencias del argumento de la función.
    sigmas = []
    for w in frecuencias: 
        sigmas.append(sigma_caso1(x,w))
    sigmas.insert(0, limite_cero(x,n))
    sigmas.append(limite_n_medios(x, n))

    sigma_max = max(sigmas) #Buscamos la sigma mayor...
    #Y ahora buscamos la frecuencia asociada a esta sigma.
    indice_sigmaMax = sigmas.index(sigma_max) #Mayor frecuencia en la que el espectro tiene su máximo
    if indice_sigmaMax == 0:
        frec_max = 0
    elif indice_sigmaMax == n*100/2:
        frec_max = n/2
    else:
      frec_max = frecuencias[indice_sigmaMax - 1] 
    #Regresamos tal frecuencia
    return frec_max


def calculo_tauMax(x):
    taus = coeficientes_tau(x)
    tau_max = max(taus)
    frec_max = taus.index(tau_max)
    return frec_max


def analisis_espectralPDL_global(n):
    """
    'n' es un entero mayor o igual a dos, y es la dimensión de los PDL para los que se
    calcula
    	1.- el vector de las frecuencias principales 1 de los PDL de dimensión n
    	2.- el vector de las frecuencias principales 0 de los PDL de dimensión n
    	3.- la ordenada al origen y pendiente (b0 _n y m0_n respect) de la recta de mínimos cuadrados (?)
    	4.- la la ordenada al origen y pendiente (b1_n y m1_n respect) de la recta de mínimos cuadrados (?)
    
    Tales coeficientes se regresan y guardan en el archivo 'data_AE.txt'
    """
    base_legendre = legendre.calculo_base(n)
    dominio_grados = [k for k in range(n)]

    frecuencias = [a/100 for a in range(int(n*100/2) + 1)] 
    #Quitamos las frecuencias extremas, que son casos especiales.
    frecuencias.pop(0)
    frecuencias.pop(-1)

    sigmasMax_n, tausMax_n = [], [] 
    for k in range(n): #iteramos en los grados de los PDL de dimensión n
        vector_legendre = base_legendre[k]
        sigmasMax_n.append(calculo_sigmaMax(vector_legendre, n, frecuencias))
        tausMax_n.append(calculo_tauMax(vector_legendre))

    b1_n, m1_n = proy.coef_RMC(dominio_grados, sigmasMax_n)
    b0_n, m0_n = proy.coef_RMC(dominio_grados, tausMax_n)

    #Vamos a guardar los valores en el diccionario definido en el script 'diccionario_RMC.py'
    #los separadores no serán comas, sino dobles espacios.

    with open('data_AE.txt', 'rb') as f:
        diccionario = pickle.load(f)
    diccionario[n] = (sigmasMax_n, tausMax_n, b0_n, m0_n, b1_n, m1_n)

    with open('data_AE.txt', 'wb') as f:
        pickle.dump(diccionario, f)

    return sigmasMax_n, tausMax_n, b0_n, m0_n, b1_n, m1_n

def tabla_informacion():
    """
    Función que se ejecuta UNA SOLA VEZ para extraer los
    datos escritos en el archivio 'data_AE.txt' y guardarlos en formato 
    de tabla html en el documento 
    'EEspectral_tabla_parametrosRectas.html'
    """
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



#TODO NUEVO, espectro basado en espacios monofrecuenciales
def Sigma(x, n, w):
    """
    Se calcula $\Sigma_{x}(w)$, con n la dimensión de la señal x.
    """
    if w == 0:
        return limite_cero(x, n)
    elif w == n/2:
        return limite_n_medios(x, n)
    else:
        return sigma_caso1(x, w)

def Tau_legendre(x, n, w):
    """
    Se calcula $\Tau_{x}(w)$, con n la dimensión de la señal x, que es un PDL.
    No tenemos entonces que dividir por la norma de x.
    """
    M = math.ceil(n/2)
    dominio = [m/n for m in range(n)]
    if w == 0:
        c_n=[math.sqrt(1/n)*np.cos(2*np.pi*t*0) for t in dominio]
        return abs(np.dot(c_n,x))
    elif (1 <= w and w <= M-1):
        c_n=[math.sqrt(2/n)*np.cos(2*np.pi*t*w) for t in dominio]
        s_n=[math.sqrt(2/n)*np.sin(2*np.pi*t*w) for t in dominio]
        return math.sqrt(np.dot(c_n, x)**2 + np.dot(s_n, x)**2 )
    elif w == M : #caso que sólo se da cuando n es par.
        c_n=[math.sqrt(1/n)*np.cos(2*np.pi*t*M) for t in dominio]
        return abs(np.dot(c_n,x))


#TODO creo que antes usas el nombre 'sigmaMax' de forma incorrecta. Sería mucho mejor usar 'frecMax0'
def guardando_sigmasYtaus_eval_en_frecPrin(n):
    """
    'n' es un entero mayor o igual a dos, y es la dimensión de los PDL para los que se
    calcula, para toda 0 \leq k \leq n-1
    	1.- la frecuencia principal 1 'w1' del PDL de dim. n y grado k junto con sigma evaluada en tal f.p.
    	2.- la frecuencia principal 0 'w0' del PDL de dim. n y grado k  junto con tau evaluada en tal f.p.
    
    Tales coeficientes se almacenan en el diccionario guardado en 'data_AE_FP.txt'
        * llave: (n,k)
        *valor: (w1, Sigma_{L^{n,k}}(w1)),  (w0, Tau_{L^{n,k}}(w0))
    """
    dominio = [m/n for m in range(n)]

    # TODO En vez de calcular todo de nuevo, voy a abrir en data_AE.txt y extraer los datos
    #	1.- la frecuencia principal 1
    #	2.- la frecuencia principal 0

    with open('data_AE_FP.txt', 'rb') as f:
        diccionario = pickle.load(f)
        
    with open('data_AE.txt', 'rb') as dataAE:
        data_AE = pickle.load(dataAE)

    #Extrayendo los vectores de frecuencias principales de los PDL de dimensión n
    F1_legendre_dim_n = data_AE[n][0]
    F0_legendre_dim_n = data_AE[n][1]
    
    #Calculamos la BLD de dimensión n.
    BLD_dim_n = legendre.calculo_base(n)

    for k in range(n):
        x = BLD_dim_n[k] #extraemos el k-ésimo vector de Legendre de dim m
        FP1_nk = F1_legendre_dim_n[k]
        FP0_nk = F0_legendre_dim_n[k]
        diccionario[(n,k)] = (FP1_nk, Sigma(x, n, FP1_nk)), (FP0_nk, Tau_legendre(x, n, FP0_nk))  

    with open('data_AE_FP.txt', 'wb') as f:
        pickle.dump(diccionario, f)


"""
#  ---------------------------------------- -- ----------------------------------------

				  FUNCIONES PARA RESPONDER A LAS
				   PREGUNTAS 2 y 3 DE MI LISTA
			      (EXLUSIVAS PARA EL ESTUDIO DE LOS PDL)

#  ---------------------------------------- -- ----------------------------------------
"""

#TODO cambiar el grafica por graficar
def grafica_analisisGlobal_k_fijo(k, graficar = True):
    """
    k es un entero entre 0 y 68.
    Se regresa la gráfica de los puntos de la forma $(n, FP0(\mathcal{L}^{n,k}))$
    y $(n, FP1(\mathcal{L}^{n,k}))$, con $k < n \leq 69$. 
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
            plt.scatter(n, tauMax_n_k, color = colores[3], marker = '^', label = r'$(n,FP0 (\mathcal{{L}}^{{n, {{{0}}} }}) )$'.format(str(k)), s = 100)
            plt.scatter(n, sigmaMax_n_k, color = colores[2], marker = 'v', label = r'$(n,FP1 (\mathcal{{L}}^{{n, {{{0}}} }}) )$'.format(str(k)), s = 100)
        else: 
            plt.scatter(n, tauMax_n_k, color = colores[3], marker = '^', s = 100)
            plt.scatter(n, sigmaMax_n_k, color = colores[2], marker = 'v', s = 100)
    
    if k == 0 or k == 1:
        axis.axhline(y = k/2, color = colores[8], label = r'$y = \frac{{ {{{0}}}  }}{{2}}$'.format(str(k)))
        graficando_puntos(3, True)
        for n in range(4, 70):
            graficando_puntos(n, False)
    else:
        for n in range(k+1, 69):
            graficando_puntos(n, False)
        graficando_puntos(69, True)
        axis.axhline(y = k/2, color = colores[8], label = r'$y = \frac{{ {{{0}}}  }}{{2}}$'.format(str(k)))
    
    axis.legend()
    axis.grid()

    plt.suptitle('Gráficas de las frecuencias principales FP para los \n PDL de grado k = '+str(k))
    axis.set_xlim(k+0.5, 69.5)
    axis.set_xlabel('Dimensión $n$')
    axis.set_ylabel('FP del polinomio discreto de Legendre '+ r'$\mathcal{{L}}^{{n, {{{0}}} }} \in \mathbb{{R}}^{{n}}$'.format(str(k)))
    if graficar == True:
        plt.tight_layout()
        return plt.show()
    else:
        fig.set_size_inches(11.25, 12.34) 
        plt.tight_layout()
        return plt.savefig(ruta_carpeta + 'k_' + str(k))


def grafica_3d_n_k_FP(N):
    """
    N: int mayor o igual a 2 y menor o igual a 69, es la dimensión máxima cuyos datos
    se van a graficar.
    
    A esta función no se le puso la opción de guardar en vez de graficar, pues por lo general
    hay que mover un poco la imagen en 3d para que se vea mejor.
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

def grafica_analisisGlobal_n_fijo(n, graficar = True):
    fig, axis = plt.subplots(2,1)
    dominio_grados = [k for k in range(n)]
    X = np.arange(0,n,0.05)
    
    #extrayendo los datos de los .txt en los que se ha guardado
    with open('data_AE.txt', 'rb') as f:
        data_AE = pickle.load(f)

    datos_espectrales_dim_n = data_AE[n]
    sigmasMax_n, tausMax_n = datos_espectrales_dim_n[0], datos_espectrales_dim_n[1] 
    b0_n, m0_n = datos_espectrales_dim_n[2], datos_espectrales_dim_n[3]
    b1_n, m1_n = datos_espectrales_dim_n[4], datos_espectrales_dim_n[5]

    etiqueta_1 = '$(k, $' + r'$FP1( \mathcal{{L}}^{{ {0}, k  }} )$'.format(str(n)) + ')'
    axis[1].scatter(dominio_grados, sigmasMax_n, marker = 'D', label = etiqueta_1, color = colores[2], s = 150)
    axis[1].plot(X, X/2, label = r'$f(t)=\frac{t}{2}$', color = 'gray', linestyle = 'dotted')
    axis[1].plot(X, b1_n+X*m1_n, color = colores[8], label = r'RMC: $l_{{ {2}, 1}}(t) = {{{0}}}t + {{{1}}}$'.format(str(round(m1_n,2)), str(round(b1_n,2)), str(n)))
    
    etiqueta_0 = '$(k, $' + r'$FP0( \mathcal{{L}}^{{ {0}, k  }} )$'.format(str(n)) + ')'
    axis[0].scatter(dominio_grados, tausMax_n, label = etiqueta_0, marker = 'D', color = colores[3] , s = 150)
    axis[0].plot(X, X/2, label = r'$f(t)=\frac{t}{2}$', color = 'gray', linestyle = 'dotted')
    axis[0].plot(X, b0_n+X*m0_n, color = colores[8], label = r'RMC: $l_{{{2}, 0 }}(t) = {{{0}}}t + {{{1}}}$'.format(str(round(m0_n,2)), str(round(b0_n,2)), str(n)))
    plt.suptitle("Frecuencias máximas encontradas en los análisis espectrales \n de los PDL de dimensión "+str(n), fontsize = 14)

    formato_axis(axis[0])
    formato_axis(axis[1])
    
    axis[0].set_xlabel('Grado $0 \leq k \leq {0}$'.format(str(n-1)))
    axis[1].set_xlabel('Grado $0 \leq k \leq {0}$'.format(str(n-1)))
    axis[0].set_ylabel(r'$FP0( \mathcal{{L}}^{{ {{{0}}}, k }} )$'.format(str(n)))
    axis[1].set_ylabel(r'$FP1( \mathcal{{L}}^{{ {{{0}}}, k }} )$'.format(str(n)))
    
    if graficar == True:
        plt.tight_layout()
        return plt.show()
    else:
        fig.set_size_inches(11.25, 12.34) 
        plt.tight_layout()
        return plt.savefig(ruta_carpeta + 'n_' + str(n))
    

def grafica_nube_b0m0_b1m1():
    """
    Gráfica de nube para responder a la pregunta 3
    """
    fig, axis = plt.subplots(1,2)
    with open('data_AE.txt', 'rb') as f:
        data_AE = pickle.load(f)

    lista_m0, lista_m1 = [], [] #listas de pendientes
    lista_b0, lista_b1 = [], [] #listas de ordenadas al origen 

    for n in range(3, 70):
        data_AE_dim_n = data_AE[n] #extraemos la información de dimensión n
        b0_n, m0_n = data_AE_dim_n[2], data_AE_dim_n[3]
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
        axis[i].axhline(y = 0.5, color = colores[8])
        axis[i].set_ylim(0, 0.7)
        axis[i].set_xlim(-1.5, 0.5)

    plt.suptitle('Gráficas de los coeficientes de las rectas \n de mínimos cuadrados para toda $n$')
    
    for i in range(2):
        axis[i].set_xlim(-1.5, 0.5)
        axis[i].set_ylim(0.3, 0.65)
        formato_axis(axis[i])

    plt.show()

def grafica_pendientes_oOrigen_RMC():
    """
    Función que grafica los puntos de la forma
    (n, m0_n) y (n, m1_n) para n entre 2 y 69.
    """
    fig, axis = plt.subplots(2,1)
    with open('data_AE.txt', 'rb') as f:
        data_AE = pickle.load(f)

    data_AE_dim_n = data_AE[3] #extraemos la información de dimensión 3
    b0_n, m0_n = data_AE_dim_n[2], data_AE_dim_n[3]
    b1_n, m1_n = data_AE_dim_n[4], data_AE_dim_n[5]
    axis[0].scatter(3, m0_n, color = colores[3], label = r'$(n, m_{n,0})$')
    axis[0].scatter(3, m1_n, color = colores[2], label = r'$(n, m_{n,2})$')
    axis[1].scatter(3, b0_n, color = colores[3], label = r'$(n, b_{n,0})$')
    axis[1].scatter(3, b1_n, color = colores[2], label = r'$(n, b_{n,1})$')

    for n in range(4, 70):
        data_AE_dim_n = data_AE[n] #extraemos la información de dimensión n
        b0_n, m0_n = data_AE_dim_n[2], data_AE_dim_n[3]
        b1_n, m1_n = data_AE_dim_n[4], data_AE_dim_n[5]
        axis[0].scatter(n, m0_n, color = colores[3])
        axis[0].scatter(n, m1_n, color = colores[2])
        axis[1].scatter(n, b0_n, color = colores[3])
        axis[1].scatter(n, b1_n, color = colores[2])

    axis[0].axhline(y = 0.5, color = 'red')
    axis[1].axhline(y = 0, color = 'red')

    axis[0].grid()
    axis[1].grid()

    axis[0].legend(loc = 'best')
    axis[1].legend(loc = 'best')

    axis[0].set_title('Pendientes de las RMC \n $l_{n,0}$ y $l_{n,1}$ para $3 \leq n \leq 69$')
    axis[1].set_title('Ordenadas al origen de las RMC \n $l_{n,0}$ y $l_{n,1}$ para $3 \leq n \leq 69$')
    
    return plt.show()


def grafica_coefEspectrales_de_frecPrincipales(k, graficar = True):
    """
    k es un entero entre 2 y 68 (inclusivo)
    Se grafican los puntos de la forma (n, tau_{n}(FP0(L^{n,k}))) y 
    (n, sigma_{n}(FP1(L^{n,k})))
    """
    ruta_carpeta = '/home/ame/GitHub/tesis-licenciatura/imagenes/graficas_analisisEspectrales/'
    fig, axis = plt.subplots(1)
    with open('data_AE_FP.txt', 'rb') as f:
        data_AE = pickle.load(f)

    for n in range(k+1, 69):
        sigmaMax_n_k = data_AE[(n,k)][0][1]
        tauMax_n_k = data_AE[(n,k)][1][1] #TODO no funciona k impar??
        axis.scatter(n, sigmaMax_n_k, color = colores[2], marker = 's')
        axis.scatter(n, tauMax_n_k, color = colores[3], marker = 's')

    sigmaMax_n_k = data_AE[(69,k)][0][1]
    tauMax_n_k = data_AE[(69,k)][1][1]
    axis.scatter(n, sigmaMax_n_k, marker = 's', color = colores[2], label = r'$(n, \sigma_{{n}}(\mathcal{{L}}^{{ n, {0}  }}, FP1(\mathcal{{L}}^{{ n, {0}  }})))$'.format(str(k)))
    axis.scatter(n, tauMax_n_k, color = colores[3], label = r'$(n, \tau_{{n}}(\mathcal{{L}}^{{ n, {0}  }}, FP0(\mathcal{{L}}^{{ n, {0}  }})))$'.format(str(k)), marker = 's')

    axis.axhline(y = 1, color = 'red')
    plt.grid()
    plt.legend()
    plt.suptitle('Valores de los coeficientes espectrales evaluados en las \n frecuencias principales de los PDL de grado '+str(k))
    if graficar == True:
      fig.tight_layout()
      return plt.show()
    else:
      #fig.set_size_inches(11.25, 12.34) 
      fig.tight_layout()
      plt.savefig(ruta_carpeta + 'coef_espectralesFP_'+ str(k) )



def espectro_Fourier(x, nombre):
    """
    Función que grafica el espectro de Fourier de una señal finita x (tipo array).
    'nombre' es una cadena que representa el nombre de la señal a analizar.
    """
    fig, axis = plt.subplots(1,2)
    n = len(x)
    grafica_taus_axis(x, n, nombre, axis[0], axis[1])
    
    return plt.show()

def espectro_Fourier_PDL(n,k):
    """
    Función que grafica el espectro de Fourier del polinomio
    discreto de Legendre (PDL) de dimensión n y grado k.
    """
    
    x = legendre.calculo_base(n)[k]
    nombre = r'\mathcal{{L}}^{{{0}}}'.format(str(n)+','+str(k)) 
    
    return espectro_Fourier(x, nombre)




if __name__ == '__main__':
  

  print('No se supone que se ejecute más este script. \n')
  print('Para usar las funciones de graficación definidas en él, hágalo a través del módulo')
  print(' "analisis_espectralesEjecuciones.py" ')

  """
  # ---------------------------------------------------------------------------------
  #                scripts que ejecuté para guardar la información en los .txt
  # ---------------------------------------------------------------------------------

  #tabla_informacion() #Sólo debía ejecutar una vez esta función para obtener el código html de la tabla.
  #Como era código plano, yo escribí un .css sencillo y le agregué algunos elementos más al archivo
  #.html.

  #------------- Guardando datos con pickle
  #Ya ejecuté este loop y guardé los datos en data_AE.txt. No lo ejecutes más!
  #for n in range(2, 70):
  #    analisis_espectralPDL_global(n)
 

  for n in range(2, 75):
    analisis_espectralPDL_global(n)

  for n in range(2, 65):
    guardando_sigmasYtaus_eval_en_frecPrin(n)
  """
