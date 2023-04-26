#TODO optimizar!
#script para realizar los dos análisis espectrales finales.

#librerías de python

import numpy as np
from numpy.linalg import norm
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import pylab
from tqdm import tqdm

#módulos personales
import base_legendreDiscreta as legendre

pi=math.pi
mpl.rcParams.update(mpl.rcParamsDefault)
colores=['#fe01b1', 'gray', '#8f99fb']

def formato_axis(axis):
  """
  Agregando los elementos que me gustan a un axis
  """
  axis.axhline(y=0, color='gray')
  axis.axvline(x=0, color='gray')
  axis.grid(True)
  axis.legend()
  
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
    
def coeficientes_tau(x):
  """
  Estos coeficientes juntan la información de cosenos y senos de la misma
  frecuencia en un solo coeficiente.
  """
  n = len(x)
  M = math.ceil(n/2)

  coef_cosenos, coef_senos = coeficientes_espectrales(x)
  taus = [coef_cosenos[0]**2] #inicializamos la lista de sigmas con la primera entrada
  coef_cosenos.pop(0)
  
  if n%2==1:
    cuadrados_cosenos = np.square(coef_cosenos)
    cuadrados_senos = np.square(coef_senos)
    for i in range(M-1):
      taus.append(cuadrados_cosenos[i]+cuadrados_senos[i])
    return taus
  
  else:
    sigma_final = coef_cosenos[M-1]**2 #guardamos el último sigma
    coef_cosenos.pop(M-1)
    cuadrados_cosenos = np.square(coef_cosenos)
    cuadrados_senos = np.square(coef_senos)
    for i in range(M-1):
      taus.append(cuadrados_cosenos[i]+cuadrados_senos[i])
    taus.append(sigma_final)
    return taus
    

def grafica_taus_axis(x, nombre, axis1, axis2):
  """
  'x' es un array de dimensión mayor a dos. 
  Esta función dibuja la gráfica de 'x' (como se definió en ??), que es una cuyo dominio es el tiempo, junto
  con la gráfica de los coeficientes sigma de x.

  """
  n=len(x)
  M = math.ceil(n/2)

  dominio_tiempo=[m/n for m in range(n)]
  taus = coeficientes_tau(x)
  cant_freq=len(taus)


  axis1.scatter(dominio_tiempo, x, color= colores[0], s=50, label='Gráfica de '+nombre)
  axis1.set_title('Gráfica de '+nombre)

  for i in range(cant_freq):
    axis2.scatter(i, taus[i], color=colores[1], s=70, marker = '*')


  axis1.set_ylabel('Coeficiente respecto a la base canónica')
  axis1.set_xlabel('Tiempo')

  axis2.set_title(r'Gráfica de los coeficientes $\tau$ de '+nombre)
  axis2.set_xlabel('Frecuencias enteras $\omega$')
  axis2.set_ylabel(r'$\tau_{{{0}}}(x, \omega)$'.format(str(n))) 

  X=np.arange(0, 1, 0.0001)

  coef_cosenos, coef_senos = coeficientes_espectrales(x) # :( código con muchas partes repetidas !
  max_w = taus.index(max(taus)) #máxima frecuencia
  axis2.scatter(max_w, max(taus), color = colores[2], s = 70, label = 'FM = '+str(max_w), marker = '*')

  if n %2 ==0:
    if max_w == 0 or max_w == M:
      axis1.plot(X, coef_cosenos[max_w]*np.cos(2*pi*max_w*X), color = colores[2], label = r'${{{0}}} \cdot cos(2 \pi \cdot {{{1}}} t) $'.format(str(round(coef_cosenos[max_w],2)), str(max_w)))
    else:
      axis1.plot(X, coef_cosenos[max_w]*np.cos(2*pi*max_w*X) + coef_senos[max_w]*np.sin(2*pi*max_w*X), color = colores[2], label = r'${{{0}}} \cdot cos(2 \pi \cdot {{{1}}} t) + {{{2}}} \cdot sen(2 \pi \cdot {{{1}}} t) $'.format(str(round(coef_cosenos[max_w],2)), str(max_w), str(round(coef_senos[max_w],2))))
  else:
    if max_w == 0:
      axis1.plot(X, coef_cosenos[max_w]*np.cos(2*pi*max_w*X), color = colores[2], label = r'${{{0}}} \cdot cos(2 \pi \cdot {{{1}}} t) $'.format(str(round(coef_cosenos[max_w],2)), str(max_w))) #redundante
    else:
      axis1.plot(X, coef_cosenos[max_w]*np.cos(2*pi*max_w*X) + coef_senos[max_w]*np.sin(2*pi*max_w*X), color = colores[2], label = r'${{{0}}} \cdot cos(2 \pi \cdot {{{1}}} t) + {{{2}}} \cdot sen(2 \pi \cdot {{{1}}} t) $'.format(str(round(coef_cosenos[max_w],2)), str(max_w), str(round(coef_senos[max_w],2))))

  formato_axis(axis1)
  formato_axis(axis2)


#  ------------------------------------- espacios monofrecuenciales ----------------------------------------

def vectores_frecuencia(n, w): #TODO quita a n del argumento.
  """
  n y w ambas de tipo int, n mayor o igual a dos, w no negativa.
  n es la dimensión, w la frecuencia 
  
  """
  n = len(x)
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
  elif (c<0 and d>0) or (c<0 and d>0):
    phi = (tan + pi)/(2*pi)
  else:
    phi = (tan + 2*pi)/(2*pi)
  return A, phi

def grafica_sigma_amplDesfase_axis_caso1(x, w, axis):
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

  axis.scatter(dominio, x, color=colores[0], s=80)
  axis.scatter(dominio, proyeccion_Pnw, color=colores[2], s=80, label = r'Gráfica de $\Pi_{{P_{{ {0}, {1} }} }}(x)$'.format(str(n),str(w)))
  
  X=np.arange(0, 1, 0.0001)
  axis.plot(X, coseno_amplDes(X), color=colores[2])
  formato_axis(axis)

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

def grafica_sigma_amplDesfase_axis_caso2(x, w, axis):
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

  axis.scatter(dominio, x, color=colores[0], s=80)
  axis.scatter(dominio, proyeccion_Pnw, color=colores[2], s=80, label=r'Gráfica de $\Pi_{{P_{{ {0}, {1} }} }}(x)$'.format(str(n),str(w)))
  
  X=np.arange(0, 1, 0.0001)
  axis.plot(X, coseno_amplDes(X), color=colores[2])
  formato_axis(axis)


def analisis_espectral_espaciosMonofrecuenciales(x, frecuencias, nombre, axis0, axis1):
  """
  x es un array
  'frecuencias' es un vector de frecuencias.
  """
  n=len(x)
  dominio_tiempo=[m/n for m in range(n)] 

  axis0.scatter(dominio_tiempo, x, color= colores[0], s=40, label = nombre)
  axis0.set_title('Gráfica de '+nombre)

  axis1.set_title('Análisis espectral generalizado')

  sigmas = []
  for w in frecuencias: 
    if w % (n/2) == 0 :
      sigma = sigma_caso2(x, w)
      sigmas.append(sigma)
    else:
      sigma = sigma_caso1(x, w)
      sigmas.append(sigma)
  
  axis1.scatter(frecuencias, sigmas, color=colores[1], s=10, marker = '*')
  sigma_max = max(sigmas)
  frec_max = frecuencias[sigmas.index(sigma_max)]
  axis1.scatter(frec_max, sigma_max, color = colores[2], label = 'FM = {0}'.format(str(frec_max)), marker = '*')

  if frec_max % (n/2) == 0:
    grafica_sigma_amplDesfase_axis_caso2(x, frec_max, axis0)
  else: 
    grafica_sigma_amplDesfase_axis_caso1(x, frec_max, axis0)

  axis0.set_xlabel('Tiempo')
  axis0.set_ylabel('Coeficiente respecto a la base canónica')
  axis1.set_xlabel('Frecuencias')
  axis1.set_ylabel(r'$\sigma_{{{0}}}(x, \omega)$'.format(str(n)))

  
  formato_axis(axis0)
  formato_axis(axis1)
 
# ------------------------------------- Análisis espectrales ------------------------------------
  
def analisis_espectrales_mostrarGrafica(x, frecuencias, nombre):
  """
  'x' (tipo array) tiene las mediciones.
  'frecuencias' (tipo array) es un array de frecuencias.
  'nombre' (tipo string) es el nombre de la señal.
  """
  n = len(x) #TODO pasar como argumento a las otras funciones!
  fig, axis = plt.subplots(2,2)
  grafica_taus_axis(x, nombre, axis[0,0], axis[0,1])
  analisis_espectral_espaciosMonofrecuenciales(x, frecuencias, nombre, axis[1,0], axis[1,1])
  fig.suptitle('Análisis espectrales de '+nombre+r'$\in \mathbb{{R}}^{{{0}}}$'.format(str(n)))
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
  
  grafica_taus_axis(x, nombre, axis[0,0], axis[0,1])
  analisis_espectral_espaciosMonofrecuenciales(x, frecuencias, nombre, axis[1,0], axis[1,1])
  fig.suptitle('Análisis espectrales de '+nombre+r'$\in \mathbb{{R}}^{{{0}}}$'.format(str(n)), fontsize = 18)
  plt.savefig("/home/ame/GitHub/tesis-licenciatura/imagenes/estudios_espectrales/prueba.png")
  
  
base_legendre = legendre.calculo_base(12)
x= base_legendre[4]
nombre = '$\mathcal{L}^{12,4}$'
frecuencias = [a/100 for a in range(501)]


#analisis_espectrales_mostrarGrafica(x, frecuencias, nombre)
analisis_espectrales_guardarGrafica(x, frecuencias, nombre)
