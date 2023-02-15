import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import legendre #En este archivo de python se guardan las BON de Legendre discretas de dimensión desde 2 hasta 6


#------------------------------------ FUNCIONES --------------------------------------------------------

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


def proyeccion(x, i):
    """
    Dado el vector x (que se supone tiene dimensión dim=len(x) de al menos 2 y a lo más 6)
    y un entero 0 \leq i \leq dim-1, se regresa la proyección de x al espacio Wi, que
    es un np.array
    """
    dim=len(x) #determinamos la dimensión
    baseLegendre=legendre.BON_L[dim] #llamamos a la base de dimensión adecuada
    proyeccion=np.zeros(dim) #inicializamos el vector proyección.
    for j in range(i+1): #sumamos las componentes necesarias.
        #print(proyeccion)
        proyeccion=proyeccion+np.dot(x, baseLegendre[j])*baseLegendre[j]
    return proyeccion


def graficas_senial_parteAfin_parteCuadratica(mediciones):
    """
    Dimensión mínima de 'mediciones':3
    Dimensión máxima de 'mediciones':6
    Esta función dibuja la gráfica de la señal 'mediciones', junto con sus partes afín
    y cuadrática, además de las rectas y parábolas de mínimos cuadrados. Comprobamos
    que las partes afín y cuadráticas son las discretizaciones de las rectas y parábolas
    de mínimos cuadrados.
    """
    #inicializamos algunos objetos
    dim=len(mediciones)
    X=np.arange(0, dim-1+0.2, 0.05)
    dominio=np.zeros(dim)
    for i in range(dim):
        dominio[i]=i
    colores_ame=['goldenrod', 'hotpink', 'rebeccapurple']

    #Graficamos la señal. Los puntos son más grandes que los demás.
    plt.scatter(dominio, mediciones, s=100, color=colores_ame[0], label='Gráfica de x')

    #Calculamos su parte afin y graficamos, junto con la recta de minimos cuadrados.
    proyW1=proyeccion(mediciones,1)
    plt.scatter(dominio, proyW1, s=50, color=colores_ame[1], label='Parte afin de x')
    b=coef_RMC(dominio,mediciones)
    y=b[0]+b[1]*dominio
    plt.plot(dominio, y, color=colores_ame[1], linestyle='dotted')
    
    #Calculamos su parte cuadratica y graficamos, junto con la parabola de minimos cuadrados.
    proyW2=proyeccion(mediciones,2)
    plt.scatter(dominio, proyW2, s=50, color=colores_ame[2], label='Parte cuadrática de x')
    parab=parab3Puntos([proyW2[0], proyW2[1], proyW2[2]])
    plt.plot(X, parab[0]*X**2 + parab[1]*X + parab[2], color=colores_ame[2], linestyle='dotted')

    #Graficamos
    plt.grid()
    plt.legend()
    plt.show()

def graficas_senial_parteCte_Afin_Cuadratica(mediciones):
    """
    Dimensión mínima de 'mediciones':3
    Dimensión máxima de 'mediciones':6
    Esta función dibuja la gráfica de la señal 'mediciones', junto con sus partes
    constante, afín
    y cuadrática, además de las rectas y parábolas de mínimos cuadrados. Comprobamos
    que las partes afín y cuadráticas son las discretizaciones de las rectas y parábolas
    de mínimos cuadrados.
    """
    #inicializamos algunos objetos
    dim=len(mediciones)
    X=np.arange(0, dim-1+0.2, 0.05)
    dominio=np.zeros(dim)
    for i in range(dim):
        dominio[i]=i
    colores_ame=['goldenrod', 'hotpink', 'rebeccapurple', 'gray']

    #Graficamos la señal. Los puntos son más grandes que los demás.
    plt.scatter(dominio, mediciones, s=100, color=colores_ame[0], label='Gráfica de x')

    #Calculamos su parte constante.
    proyW0=proyeccion(mediciones,0)
    plt.scatter(dominio, proyW0, s=50, color=colores_ame[3], label='Parte afin de x')
    plt.plot(dominio, proyW0, color=colores_ame[3], linestyle='dotted')

    #Calculamos su parte afin y graficamos, junto con la recta de minimos cuadrados.
    proyW1=proyeccion(mediciones,1)
    plt.scatter(dominio, proyW1, s=50, color=colores_ame[1], label='Parte afin de x')
    b=coef_RMC(dominio,mediciones)
    y=b[0]+b[1]*dominio
    plt.plot(dominio, y, color=colores_ame[1], linestyle='dotted')
    
    #Calculamos su parte cuadratica y graficamos, junto con la parabola de minimos cuadrados.
    proyW2=proyeccion(mediciones,2)
    plt.scatter(dominio, proyW2, s=50, color=colores_ame[2], label='Parte cuadrática de x')
    parab=parab3Puntos([proyW2[0], proyW2[1], proyW2[2]])
    plt.plot(X, parab[0]*X**2 + parab[1]*X + parab[2], color=colores_ame[2], linestyle='dotted')

    #Graficamos
    plt.grid()
    plt.legend()
    plt.show()

#3 de diciembre ------------------------------------------------------------------------

def graficas_senial_rectaMC(mediciones):
    """
    Dimensión mínima de 'mediciones':3
    Dimensión máxima de 'mediciones':6
    Esta función dibuja la gráfica de la señal 'mediciones', junto con su parte afín
    """
    #inicializamos algunos objetos
    dim=len(mediciones)
    X=np.arange(0, dim-1+0.2, 0.05)
    dominio=np.zeros(dim)
    for i in range(dim):
        dominio[i]=i
    colores_ame=['goldenrod', 'hotpink', 'rebeccapurple']

    #Graficamos la señal. Los puntos son más grandes que los demás.
    plt.scatter(dominio, mediciones, s=100, color=colores_ame[0], label='Conjunto $G_{x}$')

    #Calculamos su parte afin y graficamos, junto con la recta de minimos cuadrados.
    b=coef_RMC(dominio,mediciones)
    y=b[0]+b[1]*dominio
    plt.plot(dominio, y, color=colores_ame[1], linestyle='dotted', label='Recta de mínimos cuadrados')
    print('Ordenada al origen: '+str(b[0]))
    print('Pendiente: '+str(b[1]))

    plt.grid()
    plt.legend()
    plt.show()
    

def graficas_senial_parteAfin(mediciones):
    """
    Dimensión mínima de 'mediciones':3
    Dimensión máxima de 'mediciones':6
    Esta función dibuja la gráfica de la señal 'mediciones', junto con su parte
    afín.
    """
    #inicializamos algunos objetos
    dim=len(mediciones)
    X=np.arange(0, dim-1+0.2, 0.05)
    dominio=np.zeros(dim)
    for i in range(dim):
        dominio[i]=i
    colores_ame=['goldenrod', 'hotpink', 'rebeccapurple', 'gray']

    plt.grid()
    #Graficamos la señal. Los puntos son más grandes que los demás.
    plt.scatter(dominio, mediciones, s=100, color=colores_ame[0], label='Gráfica de $x \in \mathbb{R}^{n}$')

    #Calculamos su parte afin y graficamos, junto con la recta de minimos cuadrados.
    proyW1=proyeccion(mediciones,1)
    print(proyW1)
    plt.scatter(dominio, proyW1, s=50, color=colores_ame[1], label='Gráfica de $\Pi_{W_{1}}(x) \in \mathbb{R}^{n}$')
    b=coef_RMC(dominio,mediciones)
    y=b[0]+b[1]*dominio
    #plt.plot(dominio, y, color=colores_ame[1], linestyle='dotted')
    #plt.legend()
    #plt.show()


#con las modificaciones necesarias para ser usada en mi ejemplo final.
def graficas_senial_parteCuadrV2(mediciones):
    """
    Dimensión mínima de 'mediciones':3
    Esta función dibuja la gráfica de la señal 'mediciones', junto con su parte cuadrática.
    """
    #inicializamos algunos objetos
    dim=len(mediciones)
    X=np.arange(0, dim-1+0.2, 0.05)
    dominio=np.zeros(dim)
    for i in range(dim):
        dominio[i]=i
    colores_ame=['goldenrod', 'hotpink', 'rebeccapurple']

    #Graficamos la señal. Los puntos son más grandes que los demás.
    plt.scatter(dominio, mediciones, s=100, color=colores_ame[0])

    #Calculamos su parte cuadratica y graficamos, junto con la parabola de minimos cuadrados.
    proyW2=proyeccion(mediciones,2)
    plt.scatter(dominio, proyW2, s=50, color=colores_ame[0])
    parab=parab3Puntos([proyW2[0], proyW2[1], proyW2[2]])
    plt.plot(X, parab[0]*X**2 + parab[1]*X + parab[2], color=colores_ame[0], linestyle='dotted')

    plt.grid()
    plt.legend()
    plt.show()
