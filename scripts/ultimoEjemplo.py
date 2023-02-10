from mpl_toolkits import mplot3d
import numpy as np
import pylab
import math
import matplotlib.pyplot as plt
import ame_flecha
import legendre
import proyecciones as proy
import figuras3d as ame3d
pi=3.1416
#---------parte en R2

#fig, axis=plt.subplots()
#axis.set_xlim(-0.2,2.2)
#axis.set_ylim(-8.2,8.2)
#
#parte_afin=math.sqrt(3)*legendre.BON_L[3][0]+math.sqrt(2)*legendre.BON_L[3][1] #es [0,1,2]
##print(parte_afin)
#
##Lista con todos los ángulos que voy a considerar.
#alpha=[pi/4, pi/3, 0, pi/3, 6*pi/17]
#signo=[1,1,0,-1,-1] #lista de signos (para determinar la región.)
#
#a2_abs=[math.sqrt(5*(math.tan(t))**2) for t in alpha] #lista de coeficientes a_2, que dependen de los ángulos.
#x=parte_afin+signo[4]*(a2_abs[4])*legendre.BON_L[3][2]
#proy.graficas_senial_parteCuadrV2(x)


#INCOMPLETA: checa esta parte en el script 'figuras3d'
#------- parte en R3
fig=plt.figure()
axis=plt.axes(projection='3d')
ame3d.dibuja_W2_de_R3(axis, 3)
ame3d.dibuja_ejes(5)
axis.quiver(0,0,0,0,1,2) #esta es la proyección común a w1 de todos los vectores aquí considerados.
#axis.scatter(0,1,2) #esta es la proyección común a w1 de todos los vectores aquí considerados.

X=np.arange(-5,5,0.05)
axis.plot((0,1,2), (1,-1,3))


axis.invert_zaxis()
plt.show()


