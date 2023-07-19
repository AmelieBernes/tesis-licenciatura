import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import flechas_2D
import base_legendreDiscreta as legendre
import distancia_espaciosWnk as degree

#TODO no habría sido mejor medir el grado de una señal usando la cosine similarity con
#los de Legendre? Tal vez no, checa qué resultados obtienes así. Calcula también como antes el grado
#de las señales que aquí graficas.

#TODO hacer una animación graficando el vector y y varios vectores cuya proyeccion a
#W sea el vector y y cuyo ángulo con y sea cada vez más cercano a cero.


def ejemplo_Julio():

    def f(t):
        return (t-2)**2 * (t-5) 
    
    dominio = [t for t in range(6)]
    y =np.array([f(t) for t in range(6)] )
    
    base_legendre = legendre.calculo_base(6)
    coeficientes_y = [np.dot(y, base_legendre[k]) for k in range(4)]
    
    suma_squares = 0 
    for i in range(4):
        suma_squares += coeficientes_y[i]**2
    
    print(coeficientes_y)
    
    fig, axis = plt.subplots()
    
    distancias_coseno = degree.distancias_espacios_Wnk(y)
    grado = degree.grado_senal_redondeado(y, distancias_coseno)
    axis.scatter(dominio, y, label = "Proyección, "  + "grado = " + str(grado), color = "hotpink" , s = 100)
    
    
    #Vamos ahora a construir señales x cuya proyección a W6,4 sea "y".
    ultimo_coeficiente = [1, 5, 6, 13]
    ultimo_coeficiente = [0.01, 1, 2.95,  5,  6, 13]
    L_6_5 = np.array( base_legendre[5] )
    for i in range(len(ultimo_coeficiente)):
        x = y + ultimo_coeficiente[i] * L_6_5
        angulo = suma_squares/( np.linalg.norm(x) *  math.sqrt( suma_squares ))
        distancias_espacios = degree.distancias_espacios_Wnk(x)
    
        print("coeficiente " + str(ultimo_coeficiente[i]))
        print(distancias_espacios)
    
        grado = degree.grado_senal_redondeado(x, distancias_espacios, 10**(-2))
        axis.scatter(dominio, x, label = "x, coef. " + str(ultimo_coeficiente[i]) +  ", similitud coseno " + str(round(angulo, 4)) + ", grado " + str(grado), s = 50)
    
    
    
    
    fig.legend()
    plt.grid()
    plt.show()



