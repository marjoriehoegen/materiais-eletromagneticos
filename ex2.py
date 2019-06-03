import math
import numpy as np
import matplotlib.pyplot as plt

# Circuito RI
# 311V, 60Hz, 9.68 Ohms

R = 9.68
f = 60
t_i = 0
t_f = 3/f # 3 periodos
passo = 0.000001

t = np.arange(0,t_f,passo)

# tensao
def V(t):
	return 311*np.sin(2*np.pi*f*t)

# corrente
I = V(t)/R

# potencia instantanea
P = V(t) * I

# tensao eficaz

n = round((1/f)/passo)
soma = 0

V_int = 0
Vsq = V(t) ** 2	

for i in range(1,n+1):
	soma = passo/2 * (Vsq[i] + Vsq[i+1])
	V_int = V_int + soma

Vef = math.sqrt(f * V_int)
print("Vef: ", Vef)

# corrente eficaz

I_int = 0
Isq = I ** 2	


for i in range(1,n+1):
	soma = passo/2 * (Isq[i] + Isq[i+1])
	I_int = I_int + soma

Ief = math.sqrt(f * I_int)
print("Ief: ", Ief)

# potencia media

P_int = 0
Isq = I ** 2	

for i in range(1,n+1):
	soma = passo/2 * (P[i] + P[i+1])
	P_int = P_int + soma

Pmed = f * P_int
print("Pmed: ", Pmed)

# plots
fig, axs = plt.subplots(3, 1, constrained_layout=True)

axs[0].plot(t,V(t))
axs[0].axhline(Vef, linestyle ='--', linewidth = 0.8)
axs[0].set_title('Tensão')
axs[0].set_ylabel('V [V]')
axs[0].set_xlabel('t [s]')

axs[1].plot(t,I,'g')
axs[1].axhline(Ief, linestyle ='--', linewidth = 0.8, color='g')
axs[1].set_title('Corrente')
axs[1].set_ylabel('I [A]')
axs[1].set_xlabel('t [s]')

axs[2].plot(t,P,'r')
axs[2].axhline(Pmed, linestyle ='--', linewidth = 0.8, color='r')
axs[2].set_title('Potência')
axs[2].set_ylabel('P [W]')
axs[2].set_xlabel('t [s]')

plt.show()

