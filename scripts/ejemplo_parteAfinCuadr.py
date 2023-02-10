import numpy as np
import matplotlib.pyplot as plt
import math
import legendre
import proyecciones as proy

"""
Script que, dado un vector de mediciones (de longitud al menos 2 y a lo más 6),
imprime una imagen con la gráfica del vector, las de su parte afín y cuadrática, 
así como las gráficas (con trazo punteado) de la recta y parábola de mínimos cuadrados
calculadas a partir de la gráfica del vector.
"""

#mediciones=[1,1,1,1,1] 
#mediciones=[-1, -0.5,0,0.5,1] 
#mediciones=[1,0.25,0,0.25,1] 
#mediciones=[-1,-0.125,0,0.125,1] 
mediciones=[1,0.0625,0,0.0625,1] 


#mediciones=[-0.5, 2.4, 1.6, 1.7, 2.3, 0]
if len(mediciones)<2 or len(mediciones)>6:
    print('Error: el vector de mediciones introducido tiene demasiados (o muy pocos) datos.')
    exit()
proy.graficas_senial_parteCte_Afin_Cuadratica(mediciones)
proy.graficas_senial_rectaMC(mediciones)
proy.graficas_senial_parteAfin(mediciones)



#print(proy.proyeccion(mediciones,3)) 
