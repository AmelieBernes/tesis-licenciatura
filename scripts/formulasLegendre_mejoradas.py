#2023


#TODO: parece que sí funciona correctamente.
#TODO: Sí use las simetrías en las entradas de los vectores para hacer menos cálculos.


#TODO: recuerda que planeas usar programación dinámica y diagramas de mariposa para optimizar el código.
#TODO: ve si es util usar 'assert'
#TODO: Deberías usar numpy.
#TODO: Cambia 'Base_legendre' por 'base_legendre' (evita mayúsculas en nombres de variables)

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

def sumatoria(n,k,m):
        limite=min(m,k) #calculamos el límite superior de la sumatoria
        factor=factorial(m)/(factorial(n-1))
        suma=0 #iniciamos la sumatoria
        for j in range(limite+1):
            suma+=(-1)**j*factorial(k+j)*factorial(n-(j+1))/((factorial(j))**2*factorial(k-j)*factorial(m-j))
        return factor*suma


def base_Legendre_dimensionImpar(n):
    """
    Dada una dimensión n (int>1) que sea impar, esta función regresa un array con n arrays, siendo estos
    los vectores de la base de Legendre discreta de la dimensión n especificada.
    """
    N=n//2 #Damos el entero N tal que n=2*N+1

    #primero, calculamos y guardamos en una matriz los límites de sumación. Recuerde que tales límites son simplemente min(k,m) #TODO: si es tan simple, ve si puedes acortar esta parte. 

    BaseLegendre=[]
    for k in range(n): #iteramos primero en la variable de grado.
        #como a estas alturas del algoritmo ya se han fijado 'n' y 'k', calculamos A_nk
        A_nk=((-1)**k)*factorial(n-1)*math.sqrt((2*k+1)/(factorial(n-k-1)*factorial(n+k)))

        vector_Legendre=[0]*n #inicializamos el vector de Legendre de grado k
        #vamos a calcular la m-ésima entrada del vector.
        for m in range(N):
            #Salvo el factor A_nk, se calcula la m-esima entrada del vector de Legendre.
            suma=sumatoria(n,k,m)
            entrada_Legendre=A_nk*suma
            #La agregamos...
            vector_Legendre[m]=entrada_Legendre
            #... y la reflejamos en la posición opuesta con el signo correcto.
            vector_Legendre[2*N-m]=(-1)**k*entrada_Legendre
        #por último, calculamos la entrada central.
        if k%2==1: #si k es impar
            vector_Legendre[N]=0
        else:
            suma=sumatoria(n,k,N)
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
            suma=sumatoria(n,k,m)
            entrada_Legendre=A_nk*suma
            #La agregamos...
            vector_Legendre[m]=entrada_Legendre
            #... y la reflejamos en la posición opuesta con el signo correcto.
            vector_Legendre[2*N-m-1]=(-1)**k*entrada_Legendre
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


print(base_Legendre(4))
#print(base_Legendre(3))
#print(base_Legendre(5))
#print(base_Legendre(10))
#print(base_Legendre(20))
#print(base_Legendre(50))
#print(base_Legendre(70))
