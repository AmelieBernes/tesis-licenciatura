import numpy as np
import matplotlib.pyplot as plt
import pylab #Para usar LaTeX en captions
import math
import flechas_2D

plt.style.use('seaborn-v0_8-poster') 
params = {"ytick.color" : "black",
          "xtick.color" : "black",
          "axes.labelcolor" : "black",
          "axes.edgecolor" : "black",
          "text.usetex" : True,
          "font.family" : "serif",
          "font.serif" : ["Computer Modern Serif"]}
plt.rcParams.update(params)

#---------------------------------------------------

fig=plt.figure()
ax=fig.add_subplot(1,1,1)

X=np.linspace(0,3.15,100)
plt.plot(X, np.cos(X), color="gray", linestyle="dotted", label="$y=cos(\\theta)$")

plt.scatter(0, 0, marker="o", color="gray", s=200, label="$v$ y $w$ son paralelos")
plt.scatter(3.14, 0, marker="o", color="gray", s=200)
plt.scatter(3.14/2, 0, marker="s", color="gray", s=200, label="$v$ y $w$ son perpendiculares")
plt.scatter(0.7, np.cos(0.7), marker="*", color="gray", s=300, label="$\\frac{\langle v, w \\rangle}{||v|| \cdot ||w||}$")


plt.text(0,0.1,"0", fontsize=17)
plt.text(3.14,0.1,"$\pi$", fontsize=17)
plt.text(3.14/2,0.1,"$\\frac{\pi}{2}$", fontsize=17)

flechas_2D.dibujar_flechas_2d(fig, ax)

plt.grid()
plt.legend()
plt.show()
