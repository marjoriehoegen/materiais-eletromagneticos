import math
import numpy as np
import matplotlib.pyplot as plt

# Calculo da permeabilidade

f = 60
t_i = 0
t_f = 3/f # 3 periodos
deltat=0.0001
t = np.linspace(t_i,t_f,(t_f-t_i)/deltat)

Bs = 1.2
Hm = 8000
Hs = 5000
H = Hm*np.sin(2*np.pi*f*t)
B = np.zeros(round((t_f-t_i)/deltat))

mi0 = 4e-7*np.pi
mi = Bs/Hs
# B = (Bs - mi0*Hs) + mi0*H

for j in range(0,len(H)):
	if H[j] < -Hs:
		B[j] = -Bs-mi0*(H[j]+Hs)
	elif H[j] > Hs:
		B[j] = Bs+mi0*(H[j]-Hs)
	else:
		B[j] = mi*H[j]


fig1 = plt.figure(1)
plt.subplot(211)
plt.plot(t,H)
plt.subplot(212)
plt.plot(t,B)

fig2 = plt.figure(2)
plt.plot(H,B)

plt.show()