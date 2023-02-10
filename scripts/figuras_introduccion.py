import numpy as np
import matplotlib.pyplot as plt
import random
from random import uniform
import math
import legendre

#fig, axis= plt.subplots(1,2) #uno para la recta, otro para la parábola

X=np.arange(0, 30.5, 0.05)
dominio=[]
for i in range(30):
    dominio.append(i)
#------------------

mediciones_parabola=[]
for i in range(30):
    x=dominio[i]
    mediciones_parabola.append((0.23*x-2)**2-10+random.uniform(-0.5, 0.5))

mediciones_recta=[]
for i in range(30):
    x=dominio[i]
    mediciones_recta.append(-0.7*x+5+random.uniform(-0.5, 0.5))

plt.scatter(dominio, mediciones_parabola, color='hotpink')
plt.scatter(dominio, mediciones_recta, color='hotpink', marker='x')
plt.grid()
plt.show()
