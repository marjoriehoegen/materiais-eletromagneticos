import math
import numpy as np
import matplotlib.pyplot as plt

# Circuito RL
# 311V, 60Hz, R = 10 Ohms, L = 26.52e-3
# tensao, corrente, potencia
# Pativa, reativa, aparente, fator de potencia

# di/dt = 1/L (V - Ri)
# euler
# i(t+delta_t) = i(t) + di/dt(delta_t)

R = 10
L = 26.52e-3
f = 60
t_i = 0
t_f = 3/f # 3 periodos
n = 1001
passo = 0.0001
deltat = (t_f - t_i)/(n-1)

t = np.linspace(t_i,t_f,n)

# tensao
def V(t):
	return 311*np.sin(2*np.pi*f*t)
V = V(t)

# corrente
i = np.zeros([n])
di = np.zeros([n])

for j in range(1,n):
	i[j] = i[j-1] + ((di[j-1])*deltat)
	di[j] = (V[j] - (R * i[j]))/L

P = V * i

# Potencia ativa
Pat = 0
soma = 0

#tc = round(t_f - (1/f))

N = round((1/f)/passo)

for j in range(0,N):
	soma = passo/2 * (P[j] + P[j+1])
	Pat = Pat + soma
Pat = f * Pat

print("Pat: ", Pat)




# plots

# 3 gráficos: V, I e P
fig, axs = plt.subplots(3, 1, constrained_layout=True)

axs[0].plot(t,V)
axs[0].set_title('Tensão')
axs[0].set_ylabel('V [V]')
axs[0].set_xlabel('t [s]')

axs[1].plot(t,i,'g')
axs[1].set_title('Corrente')
axs[1].set_ylabel('I [A]')
axs[1].set_xlabel('t [s]')

axs[2].plot(t,P,'r')
axs[2].set_title('Potência')
axs[2].set_ylabel('P [W]')
axs[2].set_xlabel('t [s]')

#plt.show()

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