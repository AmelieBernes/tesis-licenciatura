import matplotlib.pyplot as plt
import numpy as np


## Construcción de la señal. Vamos a tomar tiempo = 1seg (una sola unidad de tiempo).
#
#sr = 2000 # 1.- Fijamos una frecuencia de muestreo
#ts = 1/sr #Unidad en la que se debe dividir un segundo para realizar el muestreo
#t = np.arange(0,1,ts)
#
##2.- Muestreamos de una suma de sinusoides de frecuencias 1, 4 y 7.
#freq =1
#x = 3*np.sin(2*np.pi*freq*t)
#
#freq =4
#x += np.sin(2*np.pi*freq*t)
#
#freq = 7
##freq = 7.3
#x += 0.5* np.sin(2*np.pi*freq*t)
#
#plt.figure(figsize = (8,6))
#plt.plot(t, x, 'r')
#plt.ylabel('Amplitud')
#plt.grid()
#
#plt.show()
#
## FFT
#from numpy.fft import fft, ifft
#X = fft(x)
#N = len(X) # N = 2000, TODO es siempre igual a dim(x)? Creo que sí! Tiene sentido que así sea. Creo que esta linea es redundante!
#n = np.arange(N) # n = [0, 1, ..., 1999]
#T = N/sr # T = 1, es la cantidad de unidades de tiempo de la medición, verdad?
#freq = n/T # Array [0, 1, 2, ... , 1999] #TODO estas son las frecuencias consideradas por la DFT? No, porque deberían ser enteras...
#print(freq)
#
#
#plt.figure(figsize = (12, 6))
#plt.subplot(121)
#
#plt.stem(freq, np.abs(X), 'b', markerfmt = ' ', basefmt = '-b')
#plt.xlabel('Frecuencia (Hz)')
#plt.ylabel('FFT Amplitude |X(freq)|') # Recuerda que el dominio de la TDF son frecuencias!
##plt.xlim(0, 10)
#plt.xlim(30, 40) #aquí no hay nada:) 
#
##Recuperando y graficando a la señal a partir de la inversa de la transformada discreta de Fourier.
#plt.subplot(122)
#plt.plot(t, ifft(X), 'r')
#plt.xlabel('Tiempo (segundos)')
#plt.ylabel('Amplitud')
#
#plt.tight_layout()
#plt.show()



#----------------------------------------------------------------------
#             Vamos a repetir el ejemplo con un muestreo no tan grande. 
#----------------------------------------------------------------------

# Construcción de la señal. Vamos a tomar tiempo = 1seg (una sola unidad de tiempo).

sr = 25 # 1.- Fijamos una frecuencia de muestreo
ts = 1/sr #Unidad en la que se debe dividir un segundo para realizar el muestreo
t = np.arange(0,1,ts)

#2.- Muestreamos de una suma de sinusoides de frecuencias 1, 4 y 7.
freq =1
x = 3*np.sin(2*np.pi*freq*t)

freq =4
x += np.sin(2*np.pi*freq*t)

freq = 7
#freq = 7.3
x += 0.5* np.sin(2*np.pi*freq*t)

plt.figure(figsize = (8,6))
print(type(x))
print(type(t))
plt.scatter(t, x, color = 'r')
plt.ylabel('Amplitud')
plt.grid()

plt.show()

# FFT
from numpy.fft import fft, ifft
X = fft(x)
N = len(X) 
n = np.arange(N) # n = [0, 1, ..., 1999]
T = N/sr # T = 1, es la cantidad de unidades de tiempo de la medición, verdad?
freq = n/T 


plt.figure(figsize = (12, 6))
plt.subplot(121)

plt.stem(freq, np.abs(X), 'b', markerfmt = ' ', basefmt = '-b')
plt.xlabel('Frecuencia (Hz)')
plt.ylabel('FFT Amplitude |X(freq)|') # Recuerda que el dominio de la TDF son frecuencias!
plt.xlim(0, 10)
#plt.xlim(30, 40) #aquí no hay nada:) 

#Recuperando y graficando a la señal a partir de la inversa de la transformada discreta de Fourier.
plt.subplot(122)
plt.scatter(t, ifft(X), color = 'r')
plt.xlabel('Tiempo (segundos)')
plt.ylabel('Amplitud')

plt.tight_layout()
plt.show()










#----------------------------------------------------------------------
#                       A N T I G U O
#----------------------------------------------------------------------


#Let's create a simple sine wave, with sampling rate=100, amplitude=1 and frequency=3.
#amplitude values are calculated every 1/100th second (sampling rate) and stored into a list called y1.

#import matplotlib.pyplot as plt
#import numpy as np
#from scipy import fft
#
#samples=100
#f=3
#x=np.arange(samples) #sample number
#y1=np.sin(2*np.pi*f*(x/samples)) #amplitude values
#plt.figure()
#plt.stem(x,y1,'r', )
#plt.plot(x,y1)
#plt.xlabel("Time -->")
#plt.ylabel("<-- Amplitude -->")
#plt.show()
#
#"""
#Now we have a sequence of amplitudes stored in list y1. We will pass this sequence to the FFT algorithm implemented by scipy. This algorithm returns a list yf of complex-valued amplitudes of the frequencies found in the signal. The first half of this list returns positive-frequency-terms, and the other half returns negative-frequency-terms which are similar to the positive ones. You can pick out any one half and calculate absolute values to represent the frequencies present in the signal. Following function takes samples as input and plots the frequency graph —
#"""
#
#def fft_plot(audio, sampling_rate):
#    n=len(audio)
#    T=1/sampling_rate
#    yf=fft(audio)
#    xf=np.linspace(0.0, 1.0/(2.0*T), n/2)
#    fig, ax=plt.subplots()
#    ax.plot(xf, 2.0/n* np.abs(yf[:n//2]))
#    plt.grid()
#    plt.xlabel("Frequency -->")
#    plt.ylabel("Magnitude")
#    return plt.show()
#
#
#fft_plot(y1,100)
#
#
#
#
#
#
