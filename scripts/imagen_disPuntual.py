import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import flechas_2D

colores_ame=['hotpink', 'rebeccapurple']

fig, axis= plt.subplots(1,2)

X=np.arange(-1.2,2.2,0.05)

P=[-1,0,1,2]
P_int=[-1,-0.25,0.5,1.25,2]

valores_funcion=[t**6-2*t**5+0.01*t**3+2 for t in P]
print(valores_funcion)
valores_funcion_int=[t**6-2*t**5+0.01*t**3+2 for t in P_int]


#gráficas 

axis[0].plot(X,X**6-2*X**5+0.01*X**3+2 , color=colores_ame[0], label='$f(t)=t^{6}-2t^{5}+0.01X^{3}+2$')
axis[0].scatter(P, valores_funcion, color=colores_ame[0], s=100)

axis[1].plot(X,X**6-2*X**5+0.01*X**3+2 , color=colores_ame[0], label='$f(t)=t^{6}-2t^{5}+0.01X^{3}+2$')
axis[1].scatter(P_int, valores_funcion_int, color=colores_ame[0], s=100)

for i in range(2):
    axis[i].legend()
    axis[i].grid(True)
    flechas_2D.dibujar_flechas_2d(fig, axis[i])



plt.show()



#-----------------------------------------
#X=np.arange(-1.2,2.2,0.05)
#P=[-1,0,1,2]
#valores_funcion=[t**6-2*t**5+0.01*t**3+2 for t in P]
#
#plt.plot(X,X**6-2*X**5+0.01*X**3+2 , color=colores_ame[0], label='$f(t)=t^{6}-2t^{5}+0.01X^{3}+2$')
#plt.scatter(P, valores_funcion, color=colores_ame[0])
#plt.grid(True)
##flechas_2D.dibujar_flechas_2d(fig, axis[0])
#plt.legend()
#
#plt.show()
