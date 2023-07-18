import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import base_legendreDiscreta as legendre
import random


def distancias_espacios_Wnk(x):
    """
    Dado un array x no cero, n:= dim(x) y se regresa un array de n puntos,
    siendo la k-ésima entrada de este (con 0 \leqk \leq n-1) la similitud
    coseno de x al espacio Wnk.
    """
    n = len(x)
    base_legendre = legendre.calculo_base(n)
    coeficientes = [np.dot(x, base_legendre[k]) for k in range(n)]

    sumas_coeficientesCuadrados = []
    suma = 0
    for k in range(n):
        suma += coeficientes[k]**2
        sumas_coeficientesCuadrados.append(suma)

    #Al final, lavariable suma guarda el valor del denominador.
    similitudes_coseno = []
    for k in range(n):
        similitudes_coseno.append(math.sqrt( sumas_coeficientesCuadrados[k]/suma  ))

    return similitudes_coseno



def grado_senal_redondeado(x, distancias_espacios, epsilon = 0.00001):
    """
    Función que, a partir de las distancias coseno del vector x 
    a los espacios de polinomios discretos (calculadas con la función
    "distancias_espacios_Wnk") indica de qué grado es la señal x.
    Se da un margen de error indicado por 'epsilon'.
    """
    cota_superior = 1 - epsilon
    if distancias_espacios[0] >= cota_superior:
        return 0
    n = len(x)
    for k in range(1,n-1):
        if distancias_espacios[k] >= cota_superior:
            return k
    #en el caso en el que ninguna de las distancias sea lo suficientemente cercana a 1
    #como para haber sido redondeada a uno, regresamos el mayor grado posible.
    return n

def formato_axis(axis):
    axis.axhline(y=0, color = 'gray')
    axis.axvline(x=0, color = 'gray')
    axis.grid(True)

def graficar_distancias_espacios_Wnk(x, nombre):
    nombre_latex = "${0}$".format(nombre)
    similitudes_coseno = distancias_espacios_Wnk(x)
    print(similitudes_coseno)
    n = len(x)
    dominio = [k for k in range(n)]

    fig, axis = plt.subplots(1,2)
    axis[0].scatter(dominio, x, color = 'hotpink')
    axis[0].set_title('Gráfica de ' + nombre_latex)

    axis[1].scatter(dominio, similitudes_coseno, color = 'indigo')
    axis[1].set_title("Distancias de " + nombre_latex  + " a los espacios " + r"$W_{{ {0}, k }}$".format(str(n)))
    
    k = grado_senal_redondeado(x, similitudes_coseno)
    fig.suptitle("Se estima que el grado de la señal " + nombre_latex + " es " + str(k))
    axis[1].scatter(k,1, color= 'lightgreen')

    axis[1].axhline(y=1, color = 'orangered')

    
    for i in range(2):
        formato_axis(axis[i])
    return plt.show()

def graficar_distancias_dePDL_espacios_Wnk(n,k):
    """
    n es la dimensión del polinomio discreto de Legendre (PDL)
    y k es su grado.
    """
    x = legendre.calculo_base(n)[k]
    graficar_distancias_espacios_Wnk(x, r"\mathcal{{L}}^{{{0}}}".format(str(n) + ', ' + (str(k))))

#------------------------------------------------

def f(t):
    return (t-2)**2 * (t-5) * (t-7)**2

def f(t):
    return (t-3)*(t-5)*(t-7)*(t-9)


#x = [f(t) for t in range(12)]
x = [f(t) + random.uniform(-1, 1) for t in range(12)]
nombre = "x"
graficar_distancias_espacios_Wnk(x, nombre)
#graficar_distancias_dePDL_espacios_Wnk(30,4)
