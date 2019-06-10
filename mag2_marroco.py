import math
import numpy as np
import matplotlib.pyplot as plt

# Calculo da permeabilidade
# funcao Marroco

f = 60
t_i = 0
t_f = 5/f # 3 periodos
deltat=0.0001
t = np.linspace(t_i,t_f,(t_f-t_i)/deltat)

alpha = 8.0
c = 2.0
ep = 1e-3
tau = 400e3
mi0 = 4e-7*np.pi

Bs = 1.6
B = Bs*np.sin(2*np.pi*f*t)

H = (B/mi0) * ( ((B**(2*alpha))/((B**(2*alpha))+tau) * (c - ep)) + ep )

fig1 = plt.figure(1)
plt.subplot(211)
plt.plot(t,H)
plt.subplot(212)
plt.plot(t,B)

fig2 = plt.figure(2)
plt.plot(H,B)

plt.show()