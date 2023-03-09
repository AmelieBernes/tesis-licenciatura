import numpy as np
import matplotlib.pyplot as plt
import pylab
import math


def coef_RMC(x,y):
    """
    dados 'x' y 'y' (ambos np.arrays del mismo tamaño), se dan la ordenada al origen y
    la pendiente de la recta de mínimos cuadrados (en ese orden) del conjunto de puntos con abscisas 
    las entradas de 'x' y ordenadas las de 'y'. Obtuve las fórmulas de internet.
    """
    dim=np.size(x)
    
    m_x=np.mean(x)
    m_y=np.mean(y)

    SS_xy=np.sum(y*x)-dim*m_y*m_x
    SS_xx=np.sum(x*x)-dim*m_x*m_x

    b_1=SS_xy/SS_xx
    b_0=m_y-b_1*m_x

    return (b_0 , b_1)

def parab3Puntos(A): #A es una lista con las ordenadas en cero, uno y dos
    """se regresa una lista con los coeficientes a, b y c de la parábola 
    y=ax^2+bx+c que pasa por los puntos con primeras coordenadas 1,2,3 (no son 0,1 y 2?)
    y con 
    segundas coordenadas las dadas por la lista A"""
    return [ (A[2]-2*A[1]+A[0])/2, ( -A[2]+4*A[1]-3*A[0])/2, A[0]]
