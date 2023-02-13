#TODO averiguar por qué esto funciona.
#plt.style.use('seaborn-v0_8-poster') 

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import pylab

def formato_amelie():
    plt.style.use('seaborn-v0_8-poster') 
    params = {"ytick.color" : "black",
              "xtick.color" : "black",
              "axes.labelcolor" : "black",
              "axes.edgecolor" : "black",
              "text.usetex" : True,
              "font.family" : "serif",
              "font.serif" : ["Computer Modern Serif"]}
    plt.rcParams.update(params)


