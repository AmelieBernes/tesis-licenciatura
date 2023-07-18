import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import flechas_2D
import base_legendreDiscreta as legendre

def f(t):
    return (t-2)**2 * (t-5) 

dominio = [t for t in range(6)]
y =np.array([f(t) for t in range(6)] )
L_12_5 = np.array(legendre.calculo_base(6)[4] )


fig, axis = plt.subplots()
axis.scatter(dominio, y, color = 'hotpink', label = "Proyección")
axis.scatter(dominio, y + 2* L_12_5, color = 'green', label = "alpha 2")
axis.scatter(dominio, y + 3* L_12_5, label = "alpha 3")
axis.scatter(dominio, y + 8* L_12_5, label = "alpha 8")
axis.scatter(dominio, y + 20* L_12_5, color = 'black', label = "alpha 20")

fig.legend()
plt.grid()
plt.show()
