import math #para usar la función sqrt
import numpy as np



"""
Versión del programa que calcula los polinomios discretos de Legendre con MÍ propia
función factorial, implementada usando dynamical programming.
"""


memo={0: 1, 1:1} #A dictionary in which I will save all the factorials computed during execution.
#The key is an int>=0, and the value it stores is the factorial of such integer.

def factorial(n, memoria=memo):
    try: 
        return memoria[n]
    except KeyError: #Error raised if the factorial of n was not yet computed. 
		#In this case we need to compute it and then store it in "memo"
        resultado=n*factorial(n-1)
        memoria[n]=resultado
        return resultado


"""**Funciones sumatoria y sumando: versión 1**"""

#SUMATORIA: VERSION 1. En la que implemento algunos de los detalles que discutimos juntos.

def sumando_V1(n,k,m,j):
  return factorial(k+j)*factorial(n-(j+1))/((factorial(j))**2*factorial(k-j)*factorial(m-j)) 


def sumatoria_V1(n,k,m):
  limite=min(m,k) #Calculamos el límite superior de la sumatoria
  factor=factorial(m)/(factorial(n-1))

  #Calculamos y guardamos los sumandos dependiendo de la paridad del índice j correspondiente.
  sumandos_pares=[sumando_V1(n,k,m,j) for j in range(0,limite+1,2)]
  sumandos_impares=[sumando_V1(n,k,m,j) for j in range(1,limite+1,2)] 

  #reordenamos los valores almacenados en los dos anteriores arrays de menor a mayor con el método "sort".
  sumandos_pares.sort() 
  sumandos_impares.sort()
  #TODO: Al momento de ejecutar a continuación la built-in function "sum", ¿estoy segura de que se suman de menor a mayor?

  return (factorial(m)/(factorial(n-1)))*(sum(sumandos_pares)-sum(sumandos_impares))

"""**Funciones sumatoria y sumando: versión 2**"""

#SUMATORIA: VERSION 2. 

def eliminar_entradas_comunes(A, B):
  """
  A y B son ambas listas.
  MEJORAR DESCRIPCION
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

#-----------------------------------------------------------------------------

def sumando_V2(n,k,m,j):
  B1=[t for t in range(m-j+1, m+1)]
  B2=[t for t in range(k-j+1, k+j+1)]
  B3=[t for t in range(n-j, n)]
  B4=[t for t in range(1, j+1)]

  #concatenamos cadenas
  numerador=B1+B2 
  denominador=B3+2*B4

  #simplificamos #TODO: ordenamos de menor a mayor?
  numerador_simplificado=eliminar_entradas_comunes(numerador, denominador)
  denominador_simplificado=eliminar_entradas_comunes(denominador, numerador)

  #Si alguna de esas listas es vacía, np. da como producto de sus entradas a 1 
  num=np.prod(numerador_simplificado)
  den=np.prod(denominador_simplificado)

  return num/den

def sumatoria_V2(n,k,m):
  limite=min(m,k) #Calculamos el límite superior de la sumatoria

  #Calculamos y guardamos los sumandos dependiendo de la paridad del índice j correspondiente.
  #Nota: Note confundas; no estás usando Python generator expressions.
  sumandos_pares=[sumando_V2(n,k,m,j) for j in range(0,limite+1,2)]
  sumandos_impares=[sumando_V2(n,k,m,j) for j in range(1,limite+1,2)] 

  #reordenamos los valores almacenados en los dos anteriores arrays de menor a mayor con el método "sort".
  sumandos_pares.sort() 
  sumandos_impares.sort()
  #TODO: Al momento de ejecutar a continuación la built-in function "sum", ¿estoy segura de que se suman de menor a mayor?

  return sum(sumandos_pares)-sum(sumandos_impares)

#----------------------------------------------------------------------------------------------

"""**Aquí inicia el cálculo de las bases de legendre discretas. 
Los algoritmos usan todas las funciones que definimos antes.**
"""

def base_Legendre_dimensionImpar(n, sumatoria=sumatoria_V2):
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

def base_Legendre_dimensionPar(n, sumatoria=sumatoria_V2):
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

#PROBANDO EL ALGORITMO CON LA FUNCIÓN SUMATORIA_V1
def base_Legendre(n, sumatoria=sumatoria_V1):
    """
    Dada una dimensión n (int>1), esta función regresa un array con n arrays, siendo estos
    los vectores de la base de Legendre discreta de la dimensión n especificada.
    """
    if n%2==0:
        return base_Legendre_dimensionPar(n)
    return base_Legendre_dimensionImpar(n)

base_Legendre(4)

#PROBANDO EL ALGORITMO CON LA FUNCIÓN SUMATORIA_V2
def base_Legendre(n, sumatoria=sumatoria_V2):
    """
    Dada una dimensión n (int>1), esta función regresa un array con n arrays, siendo estos
    los vectores de la base de Legendre discreta de la dimensión n especificada.
    """
    if n%2==0:
        return base_Legendre_dimensionPar(n)
    return base_Legendre_dimensionImpar(n)

base_Legendre(4)

"""**Comparando los dos enfoques**"""

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
    legendreV1=base_Legendre(n, sumatoria=sumatoria_V1)
    legendreV2=base_Legendre(n, sumatoria=sumatoria_V2)
    for k in range(n):
      if legendreV1[k]!=legendreV1[k]:
        print("Primera discrepancia de resultados en la dimensión "+str(n)+" y grado"+str(k))
        break
    print("Dimensión "+str(n)+" examinada.")
  print("Los resultados coincidieron todos!")

#print(comparando_resultados(30)) #Los resultados coincidieron todos!
#print(comparando_resultados(60)) #Los resultados coincidieron todos!
#print(comparando_resultados(100)) #Los resultados coincidieron todos!

#--------- Probando los tiempos de ejecución del programa
import sys 

#print(sys.getsizeof(base_Legendre(100, sumatoria=sumatoria_V1))) #904 en colab. Es mucho? pregunta a Javier
#print(sys.getsizeof(base_Legendre(300, sumatoria=sumatoria_V1))) #OverflowError: int too large to convert to float

import cProfile

cProfile.run('base_Legendre(50, sumatoria=sumatoria_V2)')


