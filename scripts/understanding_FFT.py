#Let's create a simple sine wave, with sampling rate=100, amplitude=1 and frequency=3.
#amplitude values are calculated every 1/100th second (sampling rate) and stored into a list called y1.

import matplotlib.pyplot as plt
import numpy as np
from scipy import fft

samples=100
f=3
x=np.arange(samples) #sample number
y1=np.sin(2*np.pi*f*(x/samples)) #amplitude values
plt.figure()
plt.stem(x,y1,'r', )
plt.plot(x,y1)
plt.xlabel("Time -->")
plt.ylabel("<-- Amplitude -->")
plt.show()

"""
Now we have a sequence of amplitudes stored in list y1. We will pass this sequence to the FFT algorithm implemented by scipy. This algorithm returns a list yf of complex-valued amplitudes of the frequencies found in the signal. The first half of this list returns positive-frequency-terms, and the other half returns negative-frequency-terms which are similar to the positive ones. You can pick out any one half and calculate absolute values to represent the frequencies present in the signal. Following function takes samples as input and plots the frequency graph —
"""

def fft_plot(audio, sampling_rate):
    n=len(audio)
    T=1/sampling_rate
    yf=fft(audio)
    xf=np.linspace(0.0, 1.0/(2.0*T), n/2)
    fig, ax=plt.subplots()
    ax.plot(xf, 2.0/n* np.abs(yf[:n//2]))
    plt.grid()
    plt.xlabel("Frequency -->")
    plt.ylabel("Magnitude")
    return plt.show()


fft_plot(y1,100)






