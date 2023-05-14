import matplotlib.pyplot as plt
import numpy as np
import math
import base_legendreDiscreta as legendre 
import funciones_figuras3d as fig_3d 
import random

#El argumento de las siguientes funciones es la frecuencia 'w' 
# n = len(x), pero lo pasamos como argumento para no calcularlo todo el tiempo.



def sigma(x, w, n):
    """
    'w' va a ser el argumento de la función.
    se regresan sigma, a, b, c.
    """
    #Calculamos a xi y eta, los factores de normalización
    xi = math.sqrt(2) * (n + np.sin(2*np.pi*w)*np.cos(2*np.pi*w*(n-1)/n)/np.sin(2*np.pi*w/n))**(-1/2)
    eta = math.sqrt(2) * (n - np.sin(2*np.pi*w)*np.cos(2*np.pi*w*(n-1)/n)/np.sin(2*np.pi*w/n))**(-1/2)

    #Calculamos los vectores que resultan de muestrear cosenos y senos de frecuencia w
    vector_coseno = [np.cos(2*np.pi*w*m/n) for m in range(n)]
    vector_seno = [np.sin(2*np.pi*w*m/n) for m in range(n)]

    a = xi * np.dot(x, vector_coseno)
    b = eta * np.dot(x, vector_seno)
    c = xi * eta * np.dot(vector_coseno, vector_seno)
    norm_cuad = np.dot(x,x)

    return [( (a**2 + b**2 - 2*a*b*c)/ (norm_cuad * (1-c**2))  )**(1/2), a, b, c ]

def sigma_extremo_0(x, n):
    norm = np.dot(x,x)**2
    return abs(sum(x))/(norm * math.sqrt(n))

def sigma_extremo_nMedios(x, n):
    vector_coseno = [np.cos(np.pi*m) for m in range(n)]
    norm = np.dot(x,x)**2
    return abs(np.dot(x, vector_coseno))/(norm * math.sqrt(n))
    
#TODO este ejemplo no es nada general. Intenta cambiar esto luego si vale la pena.
def encontrar_frecuencia_real(x, w, Fs, T):
    L = Fs * T
    fig, axis = plt.subplots(2,1)
    tiempos_mediciones = np.arange(0, T, 1/Fs)
    axis[0].scatter(tiempos_mediciones, x, color = 'deeppink')
    axis[1].scatter(tiempos_mediciones, x, color = 'deeppink')
    
    X = np.arange(0, T, 0.001)
    axis[0].plot(X, 0.28*np.cos(2*np.pi*w*X) - 0.49 *np.sin(2*np.pi*w*X), color = 'magenta')

    #Calculemos ahora la frecuencia real
    w_real = Fs*w/6 
    Y = np.arange(0, 5/Fs, 0.001)
    axis[1].plot(X, 0.28*np.cos(2*np.pi*w_real*X) - 0.49 *np.sin(2*np.pi*w_real*X), color = 'magenta')
    #axis[1].set_title()

    for i in range(2):
        axis[i].grid()


    return plt.show()

if __name__ == '__main__':
    
    base = legendre.calculo_base(6)
    puntos_iniciales = [base[4][k] for k in range(6)]
    x = base[4]  #inicializamos la mediciones con el PDL de dimensión 6 y grado 4
    #para mi ejemplo, me falta agregar 34 mediciones:
    for i in range(34):
        x.append(puntos_iniciales[i%4] + np.random.uniform(-0.3,0.3))
    encontrar_frecuencia_real(x, 2, 10, 4)
    #vector_coseno = [np.cos(2*np.pi*w*m/n) for m in range(n)]
    #vector_seno = [np.sin(2*np.pi*w*m/n) for m in range(n)]
    #print(vector_coseno)
    #print(vector_seno)


    #vector_coseno = [np.cos(2*np.pi*wPrima*m/n) for m in range(n)]
    #vector_seno = [np.sin(2*np.pi*wPrima*m/n) for m in range(n)]
    #print(vector_coseno)
    #print(vector_seno)

"""
if __name__=='__main__':
    fig, axis = plt.subplots()
    #x = legendre.calculo_base(15)[1] #aquí sí se tiene discontinuidad a la derecha de 0
    base = legendre.calculo_base(8)
    x = base[1]
    #x = base[4]
    #x = legendre.calculo_base(10)[0] #aquí no se tiene discontinuidad a la derecha de 0
    n = len(x)
    W = np.arange(0.0001,  1, 0.001)
    #W = np.arange(0.0001,  n/2, 0.001)
    Xis = [ math.sqrt(2) * (n + np.sin(2*np.pi*w)*np.cos(2*np.pi*w*(n-1)/n)/np.sin(2*np.pi*w/n))**(-1/2) for w in W ]

    Etas =[math.sqrt(2) * (n - np.sin(2*np.pi*w)*np.cos(2*np.pi*w*(n-1)/n)/np.sin(2*np.pi*w/n))**(-1/2) for w in W]
    #axis.scatter(W, Xis, color = 'hotpink') 
    cos_sen= [np.dot(np.sin(2*np.pi*w), np.cos(2*np.pi*w)) for w in W]
    #axis.scatter(W, cos_sen, color = 'mediumslateblue')

    def cab(n, w):
        return np.sin(2*np.pi*w)*np.cos(2*np.pi*w*(n-1)/n)/np.sin(2*np.pi*w/n)
    #axis.scatter(W, cab(n, W), color = 'lightseagreen')

    #axis.scatter(W, Xis, color = 'green', s = 5, marker = '*', zorder = 1)
    #axis.axhline(y = 1/math.sqrt(n))
    #axis.axhline(y = math.sqrt(2)/math.sqrt(n))
    
    sigmas, A, B, C = [], [], [], []
    for w in W:
        sig, a, b, c = sigma(x, w, n)
        sigmas.append(sig)
        A.append(a)
        B.append(b)
        C.append(c)
    axis.scatter(W, sigmas, color = 'green', s = 5, marker = '*', zorder = 1)
    axis.scatter(W, A, color = 'magenta', s = 5, marker = '*', zorder = 2)
    axis.scatter(W, B, color = 'orange', s = 5, marker = '*', zorder = 2)
    axis.axhline(y = sum(x)/math.sqrt(n), color = 'red') #TODO esto es lo que creo que es el límite por 0.

   
    dominio = [m/n for m in range(n)]
    axis.scatter(dominio, x)
    amelie = [np.cos(np.pi*m) for m in range(n)]
    #axis.axhline(y = np.dot(x, amelie)/math.sqrt(n), color = 'blue') #TODO este es lo que creo el límite por n/2
    

    plt.grid()
    plt.show()
"""
