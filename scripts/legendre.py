import numpy as np
import math

#Bases de Legendre discretas (para dimensiones de 2 hasta 6). Cada una es un array, y los vectores de las bases son np.array (para poder hacer operaciones disponibles en numpy con ellas.)

legendre2=[
np.array([1/math.sqrt(2),1/math.sqrt(2)]),
np.array([-1/math.sqrt(2),1/math.sqrt(2)]),
]

legendre3=[
np.array([1/math.sqrt(3),1/math.sqrt(3),1/math.sqrt(3)]),
np.array([-1/math.sqrt(2),0,1/math.sqrt(2)]),
np.array([1/math.sqrt(6), -math.sqrt(2/3), 1/math.sqrt(6)])
]

legendre4=[
np.array([1/2, 1/2, 1/2, 1/2]),
np.array([-3/(2*math.sqrt(5)),-1/(2*math.sqrt(5)),1/(2*math.sqrt(5)),3/(2*math.sqrt(5))]),
np.array([1/2, -1/2, -1/2, 1/2]),
[-1/(2*math.sqrt(5)), 3/(2*math.sqrt(5)),-3/(2*math.sqrt(5)),1/(2*math.sqrt(5))]
        ]

legendre5=[
np.array([1/math.sqrt(5), 1/math.sqrt(5), 1/math.sqrt(5), 1/math.sqrt(5), 1/math.sqrt(5) ]),
np.array([-math.sqrt(2/5), -1/math.sqrt(10), 0, 1/math.sqrt(10), math.sqrt(2/5)]),
np.array([math.sqrt(2/7), -1/math.sqrt(14), -math.sqrt(2/7), -1/math.sqrt(14), math.sqrt(2/7)]),
np.array([-1/math.sqrt(10), math.sqrt(2/5), 0, -math.sqrt(2/5), 1/math.sqrt(10)]),
np.array([1/math.sqrt(70), -2*math.sqrt(2/35), 3*math.sqrt(2/35), -2*math.sqrt(2/35),1/math.sqrt(70)])
]

legendre6=[
np.array([1/math.sqrt(6),1/math.sqrt(6),1/math.sqrt(6),1/math.sqrt(6),1/math.sqrt(6),1/math.sqrt(6)]),
np.array([-math.sqrt(5/14), -3/math.sqrt(70), -1/math.sqrt(70), 1/math.sqrt(70), 3/math.sqrt(70),math.sqrt(5/14)]),
np.array([5/(2*math.sqrt(21)), -1/(2*math.sqrt(21)), -2/(math.sqrt(21)),-2/(math.sqrt(21)), -1/(2*math.sqrt(21)),5/(2*math.sqrt(21))]),
np.array([-math.sqrt(5)/6, 7/(6*math.sqrt(5)), 2/(3*math.sqrt(5)), -2/(3*math.sqrt(5)),- 7/(6*math.sqrt(5)),math.sqrt(5)/6]),
np.array([1/(2*math.sqrt(7)), -3/(2*math.sqrt(7)),1/(math.sqrt(7)),1/(math.sqrt(7)),-3/(2*math.sqrt(7)),1/(2*math.sqrt(7))]),
           [-1/(6*math.sqrt(7)),5/(6*math.sqrt(7)),-5/(3*math.sqrt(7)),5/(3*math.sqrt(7)),-5/(6*math.sqrt(7)),1/(6*math.sqrt(7))] 
]

#-----------------------------------------------------------------------------------------------------------
#Almacenamos todas las bases en un diccionario. La clave es la dimensión.
BON_L={
        2: legendre2,
        3: legendre3,
        4: legendre4,
        5: legendre5,
        6: legendre6,
}




