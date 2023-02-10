#TODO: parece que sí funciona correctamente.
#TODO: Sí use las simetrías en las entradas de los vectores para hacer menos cálculos.


#TODO: recuerda que planeas usar programación dinámica y diagramas de mariposa para optimizar el código.
#TODO: ve si es util usar 'assert'
#TODO: Deberías usar numpy.

import math


#--- Función factorial, escrita por mí, usando un diccionario como memoria.
# sacrificamos memoria a cambio de obtener rapidez.

memo={0: 1, 1:1} #inicializamos la memoria. Guardamos los primeros valores de factorial.
#La llave es un int>=0 n, el valor es el factorial de n (la llave).


def factorial(n, memoria=memo):
    """
    n es un int>=0. 'memoria' es un diccionario en el que se almacenan los factoriales
    calculados por esta función cuando ha sido llamada con distintos valores de n.
    """
    try: 
        return memoria[n]
    except KeyError: #si la llave 'n' no está en 'memoria', es porque no se ha calculado el valor n!
        #calculamos pues n! y después guardamos el valor en 'memoria'
        resultado=n*factorial(n-1)
        memoria[n]=resultado
        return resultado
#--------------------------------

def calcular_sumando_Legendre(n,k,m,limite):
        suma=0 #inicializamos la suma
        factor=factorial(m)/(factorial(n-1))
        for j in range(limite):
            suma+=factor*(-1)**j*factorial(k+j)*factorial(n-(j+1))/((factorial(j))**2*factorial(k-j)*factorial(m-j))
        return suma



def base_Legendre_dimensionImpar(n):
    """
    Dada una dimensión n (int>1) que sea impar, esta función regresa un array con n arrays, siendo estos
    los vectores de la base de Legendre discreta de la dimensión n especificada.
    """
    N=n//2
    matriz_limite_sumacion=[]
    for m in range(N+1):
        row=[] #inicializamos la m-ésima fila.
        #ahora, lo llenamos.
        for k in range(n):
            row.append(min(m,k))
        #por último, lo agregamos a la matriz_limite_sumacion
        matriz_limite_sumacion.append(row)
    
    BaseLegendre=[]
    for k in range(n): #iteramos primero en la variable de grado.
        A_nk=((-1)**k)*factorial(n-1)*math.sqrt((2*k+1)/(factorial(n-k-1)*factorial(n+k)))
        vector_Legendre=[0]*n #inicializamos el vector de Legendre de grado k
        #vamos a calcular la m-ésima entrada del vector.
        for m in range(N):
            #se calcula la m-esima entrada del vector de Legendre (salvo el factor A_nk).
            suma=calcular_sumando_Legendre(n,k,m,matriz_limite_sumacion[m][k]+1)
            #La agregamos...
            vector_Legendre[m]=A_nk*suma
            #... y la reflejamos en la posición opuesta con el signo correcto.
            vector_Legendre[2*N-m]=(-1)**k*A_nk*suma
        #por último, calculamos la entrada central.
        if k%2==1: #si k es impar
            vector_Legendre[N]=0
        else:
            suma=calcular_sumando_Legendre(n,k,N,matriz_limite_sumacion[N][k]+1)
            vector_Legendre[N]=(-1)**k*A_nk*suma
        BaseLegendre.append(vector_Legendre)
        
    return BaseLegendre

def base_Legendre_dimensionPar(n):
    """
    Dada una dimensión n (int>1) que sea par, esta función regresa un array con n arrays, siendo estos
    los vectores de la base de Legendre discreta de la dimensión n especificada.
    """
    N=n//2
    matriz_limite_sumacion=[]
    for m in range(N+1):
        row=[] #inicializamos la m-ésima fila.
        #ahora, lo llenamos.
        for k in range(n):
            row.append(min(m,k))
        #por último, lo agregamos a la matriz_limite_sumacion
        matriz_limite_sumacion.append(row)
    
    BaseLegendre=[]
    for k in range(n): #iteramos primero en la variable de grado.
        A_nk=((-1)**k)*factorial(n-1)*math.sqrt((2*k+1)/(factorial(n-k-1)*factorial(n+k)))
        vector_Legendre=[0]*n #inicializamos el vector de Legendre de grado k
        #vamos a calcular la m-ésima entrada del vector.
        for m in range(N):
            #se calcula la m-esima entrada del vector de Legendre (salvo el factor A_nk).
            suma=calcular_sumando_Legendre(n,k,m,matriz_limite_sumacion[m][k]+1)
            #La agregamos...
            vector_Legendre[m]=A_nk*suma
            #... y la reflejamos en la posición opuesta con el signo correcto.
            vector_Legendre[2*N-m-1]=(-1)**k*A_nk*suma
        BaseLegendre.append(vector_Legendre)
        
    return BaseLegendre


def base_Legendre(n):
    """
    Dada una dimensión n (int>1), esta función regresa un array con n arrays, siendo estos
    los vectores de la base de Legendre discreta de la dimensión n especificada.
    """
    if n%2==0:
        return base_Legendre_dimensionPar(n)
    return base_Legendre_dimensionImpar(n)


#print(base_Legendre(2))
#print(base_Legendre(3))
print(base_Legendre(10))
#print(base_Legendre(20))
#print(base_Legendre(50))
#print(base_Legendre(70))
