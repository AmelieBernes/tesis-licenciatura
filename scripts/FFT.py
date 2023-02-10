import math
import cmath
import numpy as np


pi=math.pi

#TODO: cambiar tuplas por arrays!
#TODO: por el momento, supongo que la última entrada de a no es cero. Esto en realidad no importa, verdad?
#TODO: no podría hacer esto con recursión? No me sale jaja
def eval_PolCoeficiente(a,x):
    """
    a es una tupla (cuyas entradas son todas de tipo float) 
    que contiene los coeficientes de un polinomio A con coeficientes reales.
    Usando la regla de Horner vamos a evaluar al polinomio A en el float x.
    Si la longitud de a es uno, entonces debe terminar en coma. Cómo puedo evitar esto=
    """
    n=len(a)
    if n==1:
        return a[0]
    #inicializamos la suma
    suma=a[n-2]+x*a[n-1]
    for i in range(1,n-1):
        suma=x*suma+a[n-2-i]
    return suma


def DFT(b):
    """
    Función que, dada la tupla b, regresa la transformada discreta de Fourier de b.
    """
    n=len(b)
    w_n=complex(math.cos(2*pi/n),math.sin(2*pi/n))
    w=1

    #Inicializamos la transformada discreta
    y=[0]*n
    for k in range(n):
        y[k]=(eval_PolCoeficiente(b, w))
        w=w*w_n
    return y


#------------------------------------------

def recursive_FFT(a):
    """
    a es un array cuya longitud es una potencia de dos.
    output: la transformada rápida de Fourier de a.
    Se requiere que len(a) sea una potencia de 2 (esto para que n/2 siempre sea nuevamente divisible por 2).
    """
    n=len(a)
    
    #Base de la recursión
    if n==1:
        return a
    #Inicializamos una raíz n-ésima primitiva de la unidad. 
    w_n=complex(math.cos(2*pi/n),math.sin(2*pi/n))
    w=1

    #Creamos dos tuplas nuevas con las entradas pares e impares (resp.) de a.
    a_par=a[1::2]
    a_impar=a[::2]
    
    #Calculamos las FFT de los dos subproblemas
    y_par=recursive_FFT(a_par)
    y_impar=recursive_FFT(a_impar)

    #Tenemos todo listo para formar la FFT de a.
    #Inicializamos la FFT de a.
    y=[0]*n

    for k in range(int(n/2)):
        y[k]=y_par[k]+w*y_impar[k]
        y[k+int(n/2)]=y_par[k]-w*y_impar[k]
        w=w*w_n #actualizamos el valor de la raíz; 'twiddle factor'

    return y


#TODO: no está bien.
def inversa_FFT(y):
    """
    y es una tupla que representa la FFT de un polinomio A(x) cuya representación
    como coeficiente es una tupla a; el output es tal tupla a.
    Se requiere que len(y) sea una potencia de 2. 
    """
    n=len(y)

    #Para reusar el código de FFT, debemos antes dividir todas las entradas de la tupla y por n.
    #TODO: deberías usar numpy arrays.
    for k in range(n):
        y[k]=y[k]/n

    #Base de la recursión
    if n==1:
        return y
    #Inicializamos una raíz n-ésima primitiva de la unidad. 
    w_n=complex(math.cos(2*pi/n),math.sin(2*pi/n))
    w=1

    #Creamos dos tuplas nuevas con las entradas pares e impares (resp.) de y.
    y_par=y[1::2]
    y_impar=y[::2]
    
    #Calculamos las FFT de los dos subproblemas
    a_par=inversa_FFT(y_par)
    a_impar=inversa_FFT(y_impar)

    #Tenemos todo listo para formar la FFT de a.
    #Inicializamos la FFT de a.
    a=[0]*n

    for k in range(int(n/2)):
        a[k]=a_par[k]+w*a_impar[k]
        a[k+int(n/2)]=a_par[k]-w*a_impar[k]
        w=w*w_n #actualizamos el valor de la raíz; 'twiddle factor'

    return a

y=[1,2,3,4,5,6,7,1]
#print(FFT(inversa_FFT(y)))
#print('------------')

#print(inversa_FFT(y))
#print('------------')
#
#print(FFT(y))



    
#---------------------------------------------------
#Efficient FFT implementations

def bit_reversed(x,n):
    """
    x y n son ambos de tipo int.
    output: un entero cuya representación en binario que, consta
    de n bits, es la opuesta a la de x. Se revisa que n bits sean suficientes para
    representar a x.
    """
    #inicializamos un array de potencias de dos y longitud n.
    potencias_dos=np.array([])
    for i in range(n):
        potencias_dos=np.append(potencias_dos, 2**i)

    binario_cadena="{0:b}".format(x) #str cuyas entradas son los dígitos en binario de x.
    m=len(binario_cadena)
    if n<m:
        return None

    binario=np.zeros(n-m) #inicializamos la representación con los ceros de la derecha necesarios.
    #terminamos de llenar a 'binario' con los dígitos en binario de x.
    for i in binario_cadena:
        binario=np.append(binario,int(i))

    return int(np.dot(potencias_dos, binario))

#print(bit_reversed(1,3))

#def array_bitReverse(array):
#    """
#    'array' es un np.array, n es un int
#    output: un array de longitud m=len(array) cuya i-ésima entrada sea bit_reversed(array[i],n)
#    """
#    n=int(math.log(len(array),2))
#    A=np.array([])
#    m=len(array)
#    for i in range(m):
#        A=np.append(A,bit_reversed(array[i],n))
#    return A

#a=[0,1,2,3,4,5,6,7]
#A=array_bitReverse(a)
#print(A)



def indices_para_bitReverse(n):
    """
    n es un int, potencia de dos.
    output: un array de longitud n cuya i-ésima entrada sea bit_reversed(i,k),
    donde k=int(math.log(n,2))
    """
    A=[]
    k=int(math.log(n,2))
    for i in range(n):
        A.append(bit_reversed(i,k))
    return A


a=[0,1,2,3,4,5,6,7]


#TODO: no me da el mismo resultado que recursive_FFT :(
def iterative_FFT(a):
    """
    a es un array cuya longitud es una potencia de dos.
    output: la transformada rápida de Fourier de a.
    Se requiere que len(a) sea una potencia de 2 (esto para que n/2 siempre sea nuevamente divisible por 2).
    """
    n=len(a)
    k=int(math.log(n,2)) #cantidad de subniveles en el árbol.
    orden=indices_para_bitReverse(n)
    A=[]
    
    for i in range(n):
        A.append(a[orden[i]])
    
    #el array 'A' contiene ahora los datos en el orden requerido para el buen funcionamiento del algoritmo.
    #ya NO trabajamos con el input inicial 'a', sino con 'A'.

    for s in range(1, k+1): #'s' es la variable con la que pasamos de una rama a la otra.
        m=2**s
        
        #inicializamos el twiddle factor.
        w_m=complex(math.cos(2*pi/m),math.sin(2*pi/m))
        for k in range(0,n,m):
            w=1 #primera raíz m-ésima de la unidad
            for j in range(int(m/2)):
                t=w*A[int(k+j+m/2)]
                u=A[k+j]
                A[k+j]=u+t #reemplazamos esa entrada de A por la DFT
                A[int(k+j+m/2)]=u-t
                w= w*w_m #update de la raíz m-ésima de la unidad

    return A


a=[0,1,2,3,4,5,6,7]        
print(recursive_FFT(a))
print(iterative_FFT(a))
