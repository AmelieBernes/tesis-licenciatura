import numpy as np
import matplotlib.pyplot as plt
#import pylab #Para usar LaTeX en captions
import math
import legendre #En este archivo de python se guardan las BON de Legendre discretas de dimensión desde 2 hasta 6

#Descripción-------------------------------------------------
#Script para dibujar las graficas de las bases almacenadas en "legendre".
#TODO: lo que quiero es que, en este script, solamente se almacenen todas las funciones
#que se requieran; el proceso de graficación en si lo voy a hacer en varios scripts, uno
#por cada dimensión.

#Nota que las entradas de todos los vectores involucrados tienen todas valores absolutos no mayores a uno.
#------------------------------------------------------------


#Problemas---------------------------------------------------
#No puedo usar el método format cuando en la cadena uso LaTeX
#------------------------------------------------------------

def dibujar_ejes(axis, dominio):
    """
    Función para dibujar los ejes en los axis. Necesito la lista 'dominio'
    para saber qué tan largo debe ser el eje x.
    """
    Y=np.arange(-1.2, 1.2)
    n=len(dominio)
    X=np.arange(-0.2,n+0.2)
    axis.plot(X, 0*X, color='black')
    axis.plot(0*Y,Y, color='black')

def dibujos_axis(axis, dominio):
    axis.legend()
    dibujar_ejes(axis, dominio)










