import math 
import numpy as np
from tqdm import tqdm


import sys
import cProfile



########## Funciones sumatoria (versión 1) y sumando (versiones 1.0 y 1.1)--------------------

#--- Función factorial, escrita por mí, usando un diccionario como memoria.-
# sacrificamos memoria a cambio de obtener rapidez.

memo={0: 1, 1:1} #inicializamos la memoria. Guardamos los primeros valores de factorial.
# La llave es un int>=0 n, el valor es el factorial de n (la llave).


def ame_factorial(n, memoria=memo):
    """
    n es un int>=0. 'memoria' es un diccionario en el que se almacenarán los factoriales
    calculados por esta función cuando ha sido llamada con distintos valores de n
    (la llave será el entero, el valor el factorial de dicho entero).
    Ya tiene guardados los factoriales de 0 y 1.
    """
    try: 
        return memoria[n]
    except KeyError: #si la llave 'n' no está en 'memoria', es porque no se ha calculado el valor n!
        #calculamos pues n! y después guardamos el valor en 'memoria'
        resultado=n*factorial(n-1)
        memoria[n]=resultado
        return resultado

#SUMANDO: VERSION 1.0 
#Que usa directamente la expresión de los sumandos como un cociente de factoriales, y para el cálculo de factoriales
#usa la built-in function de Python "factorial", que está en el módulo "math".

from math import factorial

def sumando_V1_0(n,k,m,j):
  return math.factorial(k+j)*math.factorial(n-(j+1))/((math.factorial(j))**2*math.factorial(k-j)*math.factorial(m-j)) 

#SUMANDO: VERSION 1.1
#Que usa directamente la expresión de los sumandos como un cociente de factoriales, y para el cálculo de factoriales
#usa la función "ame_factorial" que definí más arriba.
def sumando_V1_1(n,k,m,j):
  return ame_factorial(k+j)*ame_factorial(n-(j+1))/((ame_factorial(j))**2*ame_factorial(k-j)*ame_factorial(m-j)) 


#SUMATORIA: VERSION 1. Además de los enteros n,k y m, esta función tiene como argumento a una de las dos funciones
#"sumando_V1_0" y "sumando_V1_1", que usará para calcular (salvo cierto factor) a los números que va a sumar.
#Se ha puesto de default a sumando=sumando_V1_0.

def sumatoria_V1(n,k,m, sumando=sumando_V1_1):
  limite=min(m,k) #Calculamos el límite superior de la sumatoria
  factor=math.factorial(m)/(math.factorial(n-1))

  #Calculamos y guardamos los sumandos dependiendo de la paridad del índice j correspondiente.
  sumandos_pares=[sumando(n,k,m,j) for j in range(0,limite+1,2)]
  sumandos_impares=[sumando(n,k,m,j) for j in range(1,limite+1,2)] 

  #reordenamos los valores almacenados en los dos anteriores arrays de menor a mayor con el método "sort".
  sumandos_pares.sort() 
  sumandos_impares.sort()

  #TODO: Al momento de ejecutar a continuación la built-in function "sum", ¿estoy segura de que se suman de menor a mayor?

  return (math.factorial(m)/(math.factorial(n-1)))*(sum(sumandos_pares)-sum(sumandos_impares))



########## Funciones sumatoria (versión 2) y sumando (versiones 2.0)-------------------- 

def eliminar_entradas_comunes(A, B):
  """
  A y B son ambas listas.
  El output es una lista cuyas entradas son todas las de A menos las que aparecen también en B.
  Nota entonces que el orden de los argumentos es importante.
  """
  copia_A=A.copy()
  copia_B=B.copy()
  for i in A:
    if i in copia_B:
      copia_A.remove(i)
      copia_B.remove(i)
  return copia_A

#SUMATORIA: VERSION 2. En la que se usa la expresión (??) de la tesis para el sumando B_{n,k,m,j}

def sumando_V2_0(n,k,m,j):
  B1=[t for t in range(m-j+1, m+1)]
  B2=[t for t in range(k-j+1, k+j+1)]
  B3=[t for t in range(n-j, n)]
  B4=[t for t in range(1, j+1)]

  #concatenamos cadenas
  numerador=B1+B2 
  denominador=B3+2*B4

  #simplificamos #TODO: deberíamos ordenar de menor a mayor?
  numerador_simplificado=eliminar_entradas_comunes(numerador, denominador)
  denominador_simplificado=eliminar_entradas_comunes(denominador, numerador)

  #Si alguna de esas listas es vacía, np. da como producto de sus entradas a 1 
  num=np.prod(numerador_simplificado)
  den=np.prod(denominador_simplificado)

  return num/den

def sumatoria_V2(n,k,m, sumando=sumando_V2_0):
  limite=min(m,k) #Calculamos el límite superior de la sumatoria

  #Calculamos y guardamos los sumandos dependiendo de la paridad del índice j correspondiente.
  #Nota: No te confundas; no estás usando Python generator expressions.
  sumandos_pares=[sumando(n,k,m,j) for j in range(0,limite+1,2)]
  sumandos_impares=[sumando(n,k,m,j) for j in range(1,limite+1,2)] 

  #reordenamos los valores almacenados en los dos anteriores arrays de menor a mayor con el método "sort".
  sumandos_pares.sort() 
  sumandos_impares.sort()

  #TODO: Al momento de ejecutar a continuación la built-in function "sum", ¿estoy segura de que se suman de menor a mayor?
  return sum(sumandos_pares)-sum(sumandos_impares)


################--------------------------------------------------------------------------------
"""

Cálculo de las bases de legendre discretas.

A continuación se dan los algoritmos

base_Legendre_dimensionImpar 
base_Legendre_dimensionPar 
calculo_base

"""
def base_Legendre_dimensionImpar(n, sumatoria=sumatoria_V1):
    """
    Input: n natural impar (variable de tipo int), que se interpreta como la dimensión del espacio,
           función "sumatoria", que se usará para calcular los sumandos (??). El valor de default de esta variable es "sumatoria_V1"
    Output: array con n arrays, siendo estos los vectores de la base de Legendre discreta de la dimensión n especificada.
    """ 
    N=n//2 #Damos el entero N tal que n=2*N+1

    BaseLegendre=[]
    for k in tqdm(range(n)): #iteramos primero en la variable de grado.
        #como a estas alturas del algoritmo ya se han fijado 'n' y 'k', calculamos A_nk
        A_nk=((-1)**k)*math.factorial(n-1)*math.sqrt((2*k+1)/(math.factorial(n-k-1)*math.factorial(n+k)))

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



def base_Legendre_dimensionPar(n, sumatoria=sumatoria_V1):
    """
    Input: n natural par (variable de tipo int), que se interpreta como la dimensión del espacio,
           función "sumatoria", que se usará para calcular los sumandos (??). El valor de default de esta variable es "sumatoria_V1"
    Output: array con n arrays, siendo estos los vectores de la base de Legendre discreta de la dimensión n especificada.
    """ 
    N=n//2
    
    BaseLegendre=[]
    for k in tqdm(range(n)): #iteramos primero en la variable de grado.
        A_nk=((-1)**k)*math.factorial(n-1)*math.sqrt((2*k+1)/(math.factorial(n-k-1)*math.factorial(n+k)))
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


def calculo_base(n, sumatoria=sumatoria_V1):
    """
    Dada una dimensión n (int>1), esta función regresa un array con n arrays, siendo estos
    los vectores de la base de Legendre discreta de la dimensión n especificada.
    """
    if n%2==0:
        return base_Legendre_dimensionPar(n)
    return base_Legendre_dimensionImpar(n)



############### Comparando los dos enfoque entre sí y el desempeño de cada uno.


def comparando_resultados(dim):
  """
  dim es un entero mayor a 2. La función calcula las bases de Legendre de dimensión 2 hasta dim
  con las funciones sumatoria_V2 y sumatoria_V3, para comparar las correspondientes entradas de las
  bases coinciden entre sí. Si, al ir iterando, encuentra una dimensión en la que no coinciden, imprime
  la dimensión y el grado en el que encontró una discrepancia. Si nunca encuentra discrepancias, termina
  de ejecutarse imprimiendo un mensaje en el que se indica esto.
  """
  print("Comparamos las bases de Legendre discretas hasta dimensión "+str(dim)+" calculadas con los dos métodos de sumación")
  for n in range(2, dim+1):
    legendreV2=base_Legendre(n, sumatoria=sumatoria_V1)
    legendreV3=base_Legendre(n, sumatoria=sumatoria_V2)
    for k in range(n):
      if legendreV2[k]!=legendreV3[k]:
        print("Primera discrepancia de resultados en la dimensión "+str(n)+" y grado"+str(k))
        break
  print("Los resultados coincidieron todos!")

#print(comparando_resultados(30)) #Los resultados coincidieron todos!


#print(sys.getsizeof(base_Legendre(100, sumatoria=sumatoria_V2)))
#cProfile.run('base_Legendre(50, sumatoria=sumatoria_V2)')

if __name__ == "__main__":
	#print(np.dot(calculo_base(20)[10], calculo_base(20)[10]))
	base=calculo_base(4)
	print(base)
