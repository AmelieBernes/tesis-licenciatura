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
        suma += vect_coseno[m]*vect_seno[m]
    return suma
    
def prod_punto_asintotico(w):
    prim_factor = 2*pi*xi(w)*eta(w)/n
    seg_factor = (n-1)*n*w/2-2*pi**2*(n-1)**2*w**3/3
    return prim_factor*seg_factor


def prod_x_sen(x, w):
    vect_seno = [np.sin(2*np.pi*w*(m)/n) for m in range(n)]
    suma = 0
    for m in range(n):
        suma += vect_seno[m] * x[m]
    return suma
    
def prod_x_sen_asint(x, w):
    fact = 2*pi*X1(x)*w/n - 4*pi**3*X3(x)*w**3/(3*n**3)
    return eta(w) * fact

def prod_x_cos(x, w):
    vect_cos = [np.cos(2*np.pi*w*(m)/n) for m in range(n)]
    suma = 0
    for m in range(n):
        suma += vect_cos[m] * x[m]
    return suma

def prod_x_cos_asint(x, w):
    fact = X0(x) - 2*pi**2*X2(x)*w**2/n**2 + 2*pi**4*X4(x)*w**4/(3*n**4)
    return xi(n) * fact
#-----------------------------------------------------------------------

#Se ve bien !!
def figura_xi_eta():
    fig, axis = plt.subplots(1,2)
    W= np.arange(0.001, 0.25, 0.001) # dominio de w
    
    axis[0].axhline(y = 1/math.sqrt(n), color = 'red')
    axis[0].plot(W, xi(W), color = 'green') #en efecto tiende a 1/sqrt(n)
    axis[0].plot(W, xi_asintotico(W), color = 'yellowgreen') #en efecto tiende a 1/sqrt(n)
    axis[0].set_title(r'$\xi_{{ {0}, \omega  }}$'.format(str(n)))
    
    axis[1].plot(W, eta(W), color = 'purple') #en efecto tiende a infty
    axis[1].plot(W, eta_asintotica(W), color = 'mediumpurple') #en efecto tiende a infty
    axis[1].set_title(r'$\eta_{{ {0}, \omega  }}$'.format(str(n)))
    
    axis[0].grid()
    axis[1].grid()
    plt.show()

# TODO No se ve bien :(
def figura_prod_punto():
    fig, axis = plt.subplots()
    W= np.arange(0.00001, 0.25, 0.001) # dominio de w
    
    axis.plot(W, prod_punto(W), color = 'green')
    axis.plot(W, prod_punto_asintotico(W), color = 'yellowgreen') 
    axis.set_title(r'$\langle c_{{ {0}, \omega }}, s_{{ {0}, \omega }}  \rangle$'.format(str(n)))
    
    axis.grid()
    plt.show()

#se ve bien !!!
def figura_seno(x):
    fig, axis = plt.subplots()
    W= np.arange(0.00001, 0.25, 0.001) # dominio de w
    
    axis.plot(W, prod_x_sen(x, W), color = 'green') 
    axis.plot(W, prod_x_sen_asint(x, W), color = 'yellowgreen')
    axis.set_title(r'$\langle x, s_{{ {0}, \omega }}  \rangle$'.format(str(n)))
    
    axis.grid()
    plt.show()

#se ve super bien!!
def figura_coseno(x):
    fig, axis = plt.subplots()
    W= np.arange(0.00001, 0.25, 0.001) # dominio de w
    
    axis.plot(W, prod_x_cos(x, W), color = 'green') 
    axis.plot(W, prod_x_cos_asint(x, W), color = 'yellowgreen')
    axis.axhline(y = X0(x)/math.sqrt(n), color = 'red', linestyle = 'dotted')
    axis.set_title(r'$\langle x, c_{{ {0}, \omega }}  \rangle$'.format(str(n)))
    
    axis.grid()
    plt.show()
n = 20
x = legendre.calculo_base(n)[8]
#figura_coseno(x)
figura_prod_punto()
#figura_xi_eta()
