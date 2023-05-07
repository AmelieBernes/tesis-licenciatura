#-------------------------------------------
#-- Descripción ----------------------------
#Script para graficar los polinomios discretos de Legendre
#junto con sus versiones continuas, es decir, los polinomios
#que discretizados en la malla Pn dan lugar al correspondiente PDL.
#-------------------------------------------
#-------------------------------------------



import numpy as np
from numpy.linalg import norm
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import pylab
from tqdm import tqdm
import os

#módulos personales
import base_legendreDiscreta as legendre

colores = ['hotpink']

def formato_axis(axis):
  """
  Agregando los elementos que me gustan a un axis
  """
  axis.axhline(y=0, color='gray')
  axis.axvline(x=0, color='gray')
  axis.grid(True)
  axis.legend()

"""
Nota: en realidad, bien pudimos haber graficado el polinomio continuo y despues muestrearlo
uniformemente para obtener el correspondiente PDL, pero no hacemos esto en el código para
comprobar que la discretización del polinomio continuo es, en efecto, el PDL.
"""

#-----------------------------------------------------------------------------
def comprobando_ortonormalidad(n,k1,k2):
    base_legendre = legendre.calculo_base(n)
    vector1 = base_legendre[k1]
    vector2 = base_legendre[k2]

    print(np.dot(vector1, vector2))


def graficaPDL_continuo(n,k, fig, axis, ver_valores = False):
    def Pnk(n,k,t):
        A = math.factorial(n-1)
        coef = A * math.sqrt( (2*k+1)/( math.factorial(n-k-1) * math.factorial(n+k) )  )
        suma = 0
        argumento = 1

        for j in range(k):
            suma += (-1)**j * math.factorial(k+j)*math.factorial(n-1-j)/( (math.factorial(j))**2 * math.factorial(k-j) * A )* argumento
            argumento = argumento * (t - j) 

        suma += (-1)**k * math.factorial(k+k)*math.factorial(n-1-k)/( (math.factorial(k))**2 * math.factorial(k-k) * A )* argumento
        resultado = coef * suma
        if k % 2 == 0:
            return resultado
        else:
            return -1 * resultado
    
    vector_legendre = legendre.calculo_base(n)[k]
    dominio = [m for m in range(n)]
    X = np.arange(0, n-1+0.01, 0.001)
    nombre = r'$\mathcal{{L}}^{{ {0} }}$'.format(str(n) + ', ' + str(k)) 

    if ver_valores == True:
        print(vector_legendre)

    if k == 0:
        axis.scatter(dominio, vector_legendre, color = colores[0], label = nombre)
        axis.plot(X, 0*X + 1/math.sqrt(n), color = colores[0])
        formato_axis(axis)
        axis.set_title('Gráfica de ' + nombre, fontsize = 12)

    else:
        axis.scatter(dominio, vector_legendre, color = colores[0], label = nombre)
        axis.plot(X, Pnk(n, k, X), color = colores[0])
        formato_axis(axis)
        axis.set_title('Gráfica de ' + nombre, fontsize = 12)


def guardar_graficaPDL_continuo(n,k, ruta):
    def Pnk(n,k,t):

        A = math.factorial(n-1)
        coef = A * math.sqrt( (2*k+1)/( math.factorial(n-k-1) * math.factorial(n+k) )  )
        suma = 0
        argumento = 1

        for j in range(k):
            suma += (-1)**j * math.factorial(k+j)*math.factorial(n-1-j)/( (math.factorial(j))**2 * math.factorial(k-j) * A )* argumento
            argumento = argumento * (t - j) 

        suma += (-1)**k * math.factorial(k+k)*math.factorial(n-1-k)/( (math.factorial(k))**2 * math.factorial(k-k) * A )* argumento
        resultado = coef * suma
        if k % 2 == 0:
            return resultado
        else:
            return -1 * resultado
    fig, axis = plt.subplots(1,1)
    X = np.arange(0, n-1+0.01, 0.001)
    nombre = r'$\mathcal{{L}}^{{ {0} }}$'.format(str(n) + ', ' + str(k)) 

    if k == 0:
        vector_legendre = legendre.calculo_base(n)[0]
        dominio = [m for m in range(n)]

        axis.scatter(dominio, vector_legendre, color = colores[0], label = nombre)
        axis.plot(X, 0*X + 1/math.sqrt(n), color = colores[0])
        formato_axis(axis)
        axis.set_title('Gráfica de ' + nombre, fontsize = 12)

        return 
    vector_legendre = legendre.calculo_base(n)[k]
    dominio = [m for m in range(n)]


    axis.scatter(dominio, vector_legendre, color = colores[0], label = nombre)
    axis.plot(X, Pnk(n, k, X), color = colores[0])
    formato_axis(axis)
    axis.set_title('Gráfica de ' + nombre, fontsize = 12)

    final_ruta = '/dim_' + str(n) + '/'
    if os.path.isdir(ruta + final_ruta):
        plt.savefig(ruta + final_ruta +  str(n) + '-' + str(k))
    else:
        os.mkdir(ruta + final_ruta)
        plt.savefig(ruta + final_ruta +  str(n) + '-' + str(k))

def ver_graficaPDL_continuo(n, k, ver_valores = True):
    fig, axis = plt.subplots(1,1)
    graficaPDL_continuo(n,k, fig, axis, ver_valores)
    return plt.show()

def guardar_graficaPDL_continuo(n,k, ruta, ver_valores):
    fig, axis = plt.subplots(1,1)
    graficaPDL_continuo(n,k, fig, axis, ver_valores)
    final_ruta = '/dim_' + str(n) + '/'
    if os.path.isdir(ruta + final_ruta):
        plt.savefig(ruta + final_ruta +  str(n) + '-' + str(k))
    else:
        os.mkdir(ruta + final_ruta)
        plt.savefig(ruta + final_ruta +  str(n) + '-' + str(k))



if __name__=='__main__':
    #ver_graficaPDL_continuo(7, 6)
    ruta_amelie = '/home/ame/GitHub/tesis-licenciatura/imagenes/PDL/'
   
    ver_graficaPDL_continuo(84, 47) #TODO sospechoso, tal vez aquí sí hay errores numéricos...

    #for n in range(2,21):
    #    for k in range(n):
    #        guardar_graficaPDL_continuo(n, k, ruta_amelie, ver_valores = False)
    #for k in range(5):
    #    guardar_graficaPDL_continuo(5, k, ruta_amelie, ver_valores = False)
    #comprobando_ortonormalidad(84, 82, 82)
