import math
import numpy as np
import matplotlib.pyplot as plt

# circuito RLC
# 110V, 60Hz, R = 20 Ohms, L = 3.18e-2, C = 120.6e-6
# Vc_i = 36 V
# Z = R + jwL - j/wC

# corrente em regime permanente
# triangulo de potencias
# FP

# If = 4.92 A
# Pat = 483 W
# Preat L = 290
# Preat C = -530
# S = 541.2

R = 20
L = 3.18e-2
C = 120.6e-6
f = 60

t_i = 0
t_f = 8/f # 5 periodos

deltat=0.00001;
n = round((t_f)/deltat)

t = np.linspace(t_i,t_f,(t_f-t_i)/deltat)

V = 110*np.sin(2*np.pi*f*t)
i = np.zeros([n])
di = np.zeros([n])
Vc = np.zeros([n])
Vc[0] = 36

for j in range(0,n-1):
	di[j]= (V[j] - R*i[j] - Vc[j])/L
	i[j+1] = i[j] + di[j]*deltat
	Vc[j+1]=Vc[j]+(i[j+1]*deltat)/C

# peúltimo período
nf = round((t_f - (1/f))/deltat)
ni = round((t_f - (2/f))/deltat)

# corrente eficaz
I_int = 0
Isq = i ** 2
soma = 0
for j in range(ni,nf-1):
	soma = deltat/2 * (Isq[j] + Isq[j+1])
	I_int = I_int + soma
Ief = math.sqrt(f * I_int)
print("Ief: ", Ief)

# potencia ativa
P = V * i
Pat = 0
soma = 0
for j in range(ni,nf):
	soma = deltat/2 * (P[j] + P[j+1])
	Pat = Pat + soma
Pat = f * Pat

print("Potência ativa: ", Pat)


# plots
fig, axs = plt.subplots(2, 1, constrained_layout=True)

axs[0].plot(t,Vc)
axs[0].set_title('Tensão no capacitor')
axs[0].set_ylabel('V [V]')
axs[0].set_xlabel('t [s]')

axs[1].plot(t,i,'g')
axs[1].set_title('Corrente')
axs[1].set_ylabel('I [A]')
axs[1].set_xlabel('t [s]')

# axs[2].plot(t,P,'r')
# axs[2].set_title('Potência')
# axs[2].set_ylabel('P [W]')
# axs[2].set_xlabel('t [s]')

plt.show()

# gráfico da tensão e corrente juntos
# fig, ax1 = plt.subplots()

# ax1.plot(t,V)
# ax1.set_xlabel('t [s]')
# ax1.set_ylabel('V [V]', color='b')
# ax1.tick_params('y', colors='b')

# ax2 = ax1.twinx()
# ax2.plot(t, i, 'g')
# ax2.set_ylabel('I [A]', color='g')
# ax2.tick_params('y', colors='g')

# fig.tight_layout()
# plt.show()