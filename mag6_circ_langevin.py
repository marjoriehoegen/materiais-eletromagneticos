import math
import numpy as np
import matplotlib.pyplot as plt

# Circuito magnético
# -V + Ri + Ndphidt = 0
# phi = int Bds
# nesse caso, S constante
# -V + Ri + NSdBdt = 0

# usando Langevin pra permeabilidade

lm = 0.94
S = 105e-6
N = 700
R = 1.5
mi0 = 4e-7*np.pi
mi = 1000*mi0
Vp = 30
f = 60
t_i = 0
t_f = 12/f #5 periodos
deltat=0.00001 #passo
t = np.linspace(t_i,t_f,(t_f-t_i)/deltat)
n = round(t_f/deltat)

# tensao senoidal
V = Vp*np.cos(2*np.pi*f*t)
# condições iniciais
I = np.zeros([n])
H = np.zeros([n])
B = np.zeros([n])

# Langevin
He = np.zeros([n])
M = np.zeros([n])
dManHe = 0
Ms = 1.5e6
alpha = 1e-4
a = 50


for j in range(0, n-1):
	dBdt = (1/(N*S))*(V[j] - R*I[j])
	B[j+1] = B[j] + dBdt*deltat

	if abs(He[j]) <= 1e-3:
		dManHe = Ms/(3*a)
	else:
		dManHe = (Ms/a)*(1-(1/np.tanh(He[j]/a))**2+(a/He[j])**2)

	M[j+1] = (B[j+1]/mi0) - H[j+1]
	deltM = M[j+1] - M[j]
	H[j+1] = deltM*((1/dManHe)-alpha) + H[j]
	I[j+1] = H[j+1] * (lm/N)
	He[j+1] = H[j+1] + alpha*M[j+1]

	I[j+1] = H[j+1] * (lm/N)


fig1 = plt.figure(1)
fig1.suptitle('V e I')
plt.subplot(211)
plt.plot(t,V)
plt.grid(True)
plt.subplot(212)
plt.plot(t,I)
plt.grid(True)

fig2 = plt.figure(2)
fig2.suptitle('B e H')
plt.subplot(211)
plt.plot(t,B)
plt.grid(True)
plt.subplot(212)
plt.plot(t,H)
plt.grid(True)

fig3 = plt.figure(3)
plt.title('BH')
plt.plot(H,B)


plt.show()