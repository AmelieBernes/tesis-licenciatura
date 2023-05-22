import matplotlib.pyplot as plt
import numpy as np
import math
import base_legendreDiscreta as legendre 
import funciones_figuras3d as fig_3d 
import random

#El argumento de las siguientes funciones es la frecuencia 'w' 
# n = len(x), pero lo pasamos como argumento para no calcularlo todo el tiempo.

pi = np.pi
# -------- MOMENTOS -------------
def X0(x):
    suma = 0
    for m in range(n):
        suma += x[m]
    return suma

def X1(x):
    suma = 0
    for m in range(n):
        suma += m * x[m]
    return suma

def X2(x):
    suma = 0
    for m in range(n):
        suma += m**2 * x[m]
    return suma

def X3(x):
    suma = 0
    for m in range(n):
        suma += m**3 * x[m]
    return suma

def X4(x):
    suma = 0
    for m in range(n):
        suma += m**4 * x[m]
    return suma

#------------------------

def xi(w):
    return math.sqrt(2) * (n + np.sin(2*np.pi*w)*np.cos(2*np.pi*w*(n-1)/n)/np.sin(2*np.pi*w/n))**(-1/2)

def xi_asintotico(w):
    num = 2*n * (w - 2*pi**2*(2*n**2-3*n+2)*w**3/(3*n**2))
    den = w - 2*pi**2*w**3/(3*n**2)
    return math.sqrt(2)*(num/den)**(-1/2)

def eta(w):
    return math.sqrt(2) * (n - np.sin(2*np.pi*w)*np.cos(2*np.pi*w*(n-1)/n)/np.sin(2*np.pi*w/n))**(-1/2)

def eta_asintotica(w):
    num = 8*pi**3*(2*n-1)*(n-1)*w**3
    den = 6*pi*n*w-4*pi**3*w**3/n
    return math.sqrt(2)*(num/den)**(-1/2)


def prod_punto(w):
    vect_coseno = [np.cos(2*np.pi*w*m/n) for m in range(n)]
    vect_seno = [np.sin(2*np.pi*w*m/n) for m in range(n)]
    suma = 0
    for m in range(n):
        suma += xi(w)*eta(w)*vect_coseno[m]*vect_seno[m]
    return suma
    
def prod_punto_asintotico(w):
    prim_factor = 2*pi*xi(w)*eta(w)/n
    seg_factor = (n-1)*n*w/2-2*pi**2*(n-1)**2*w**3/3
    return prim_factor*seg_factor


def prod_punto_asintotico_v2(w):
    a = 2*math.sqrt(2)*pi*(n-1)/(n * math.sqrt(n))
    b = 1.5*pi*n**3*w**3/(8*pi**3*(2*n-1)*(n-1)*w**3)
    return a*b**(1/2)

def prod_x_sen(x, w):
    vect_seno = [np.sin(2*np.pi*w*(m)/n) for m in range(n)]
    suma = 0
    for m in range(n):
        suma += eta(w) * vect_seno[m] * x[m]
    return suma
    
def prod_x_sen_asint(x, w):
    fact = 2*pi*X1(x)*w/n - 4*pi**3*X3(x)*w**3/(3*n**3)
    return eta(w) * fact

#TODO no se ve bien....
def prod_x_sen_asint_v2(x, w):
    num = 24*pi**3*X1(x)**2*w**3*math.sqrt(2)
    den = 8*pi**3*(2*n-1)*(n-1)
    return math.sqrt(2)*np.sqrt(num/den)

#TODO no se ve bien....
def prod_x_sen_asint_v2(x, w):
    a = 6*pi*n*w-4*pi**3*w**3/n
    b = 8*pi**3*(2*n-1)*(n-1)*w**3

    c = 4*pi*2*X1(x)**2*w*-16*pi**4*X1(x)*X3(x)*w**4/(3*n**4)
    return np.sqrt(2)*np.sqrt(a/b)*np.sqrt(c)


def prod_x_cos(x, w):
    vect_cos = [np.cos(2*np.pi*w*(m)/n) for m in range(n)]
    suma = 0
    for m in range(n):
        suma += xi(n) * vect_cos[m] * x[m]
    return suma

def prod_x_cos_asint(x, w):
    fact = X0(x) - 2*pi**2*X2(x)*w**2/n**2 + 2*pi**4*X4(x)*w**4/(3*n**4)
    return xi(n) * fact
#-----------------------------------------------------------------------

#Se ve bien !!
def figura_xi_eta():
    fig, axis = plt.subplots(1,2)
    W= np.arange(0.001, 0.25, 0.001) # dominio de w
    
    axis[0].axhline(y = 1/math.sqrt(n), color = 'red', linestyle = 'dotted')
    axis[0].plot(W, xi(W), color = 'green') #en efecto tiende a 1/sqrt(n)
    axis[0].plot(W, xi_asintotico(W), color = 'yellowgreen') #en efecto tiende a 1/sqrt(n)
    axis[0].set_title(r'$\xi_{{ {0}, \omega  }}$'.format(str(n)))
    
    axis[1].plot(W, eta(W), color = 'purple') #en efecto tiende a infty
    axis[1].plot(W, eta_asintotica(W), color = 'mediumpurple') #en efecto tiende a infty
    axis[1].set_title(r'$\eta_{{ {0}, \omega  }}$'.format(str(n)))
    
    axis[0].grid()
    axis[1].grid()
    plt.show()

#Se ve bien!
def figura_prod_punto():
    fig, axis = plt.subplots()
    W= np.arange(0.00001, 0.25, 0.001) # dominio de w
    
    axis.plot(W, prod_punto(W), color = 'green')
    axis.plot(W, prod_punto_asintotico(W), color = 'yellowgreen') 
    axis.plot(W, prod_punto_asintotico_v2(W), color = 'hotpink') 
   
    a = math.sqrt(6*(n-1))
    b = 2*math.sqrt(2*n-1)
    axis.axhline(y = a/b, color = 'red', linestyle = 'dotted')

    axis.set_title(r'$\langle c_{{ {0}, \omega }}, s_{{ {0}, \omega }}  \rangle$'.format(str(n)))
    
    axis.grid()
    plt.show()


#TODO creo que aquí hay un problema.
def figura_seno(x):
    fig, axis = plt.subplots()
    W= np.arange(0.00001, 0.25, 0.0001) # dominio de w
    
    axis.plot(W, prod_x_sen(x, W), color = 'green') 
    axis.plot(W, prod_x_sen_asint(x, W), color = 'yellowgreen')
    #axis.plot(W, prod_x_sen_asint_v2(x, W), color = 'blue', zorder = 3)
    
    discrim = 2*pi*X1(x)*0.01/n - 4*pi**3*X3(x)*0.01**3/(3*n**3)
    if discrim > 0:
        axis.axhline(y = math.sqrt(6*X1(x)**2/(n*(2*n-1)*(n-1))) , color = 'red', linestyle = 'dotted' )
    else: 
        axis.axhline(y = - math.sqrt(6*X1(x)**2/(n*(2*n-1)*(n-1))) , color = 'red', linestyle = 'dotted' )
    axis.set_title(r'$\langle x, s_{{ {0}, \omega }}  \rangle$'.format(str(n)))
    
    axis.grid()
    plt.show()

#Perfecto !
def figura_coseno(x):
    fig, axis = plt.subplots()
    W= np.arange(0.00001, 0.25, 0.0001) # dominio de w
    
    axis.plot(W, prod_x_cos(x, W), color = 'green') 
    axis.plot(W, prod_x_cos_asint(x, W), color = 'yellowgreen')
    axis.axhline(y = X0(x)/math.sqrt(n), color = 'red', linestyle = 'dotted') 
    axis.set_title(r'$\langle x, c_{{ {0}, \omega }}  \rangle$'.format(str(n)))
    
    axis.grid()
    plt.show()
#n = 50
#x = legendre.calculo_base(n)[1]
#figura_seno(x)
#figura_coseno(x)
#
#n = 13
##x = legendre.calculo_base(n)[12]
##figura_seno(x)
##figura_coseno(x)
#n = 15
#x = legendre.calculo_base(15)[5]

x = np.array([-1,2,5,6,6,7.3,2.4])
n = len(x)
figura_seno(x)
figura_coseno(x)
#TODO noto que para los PDL de grado alto las aproximaciones no parecen muy creíbles. Serán errores numéricos?


#x = legendre.calculo_base(n)[1]
#figura_seno(x)
#figura_coseno(x)


n = 150 
#figura_prod_punto()
#figura_xi_eta()
