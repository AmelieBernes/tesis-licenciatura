import matplotlib.pyplot as plt
import numpy as np
import pylab


"""

#pruebas str python
#16 Febrero
X=np.linspace(-1,1,100)
Y=2*X

plt.plot(X,Y)
#plt.title("$x^{2}+5$") #Sí funciona

N=20
k=13
potencia=0.3
subindice=-12
plt.title("$x_{{ {1} }}^{{ {0} }}+5$".format(str(subindice), str(potencia)))  #Sí funciona

plt.title(r"$x_{{ \sigma }}^{{ {0} }}+5$".format(str(subindice), str(potencia)))  #Sí funciona
plt.title(r"$\mathbb{L}^{3}$")  #Sí funciona
plt.title(r"$\mathcal{{L}}^{{ {0} }}_{{ 80 }}$".format(str(20) )  ) #Sí funciona.
plt.title(r"$\mathcal{{L}}^{{ {0}, {1} }}_{{ 80 }}$".format(str(20), str(k) )  ) #Sí funciona


plt.title(r"$\sigma_{{\omega}}^{{0}}( \mathcal{{ L }}^{{ {0} , {1} }} )$".format(str(N), str(k) ) ) #Sí funciona



plt.show()
"""

#----------------------------------------------------------------------------------
#17 Febrero


#fig, axis= plt.subplots(1,2)

#X=np.linspace(-1,1,100)
#Y=2*X
#Z=0.5*X 

#axis[0].plot(X, Y)
#axis[1].plot(X, Z)

#plt.show()
#plt.clf()
#plt.show()

#----------------------------------------------------------------------------------
#18 Febrero
#r=0
#for i in range(8):
#	for j in range(3):
#		r+=1
#		print(i+j)
#		if r==4:
#			break


print('-----------------------')

#Y si mejor pongo dos break? Para salir de los dos for loops!
r=0
for i in range(6):
	if r==5:
		break
	for j in range(2):
		print("i= "+str(i)+", j="+str(j))
		r+=1
		if r==5:
			break 
	#r+=1

