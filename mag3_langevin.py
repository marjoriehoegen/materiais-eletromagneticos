import math
import numpy as np
import matplotlib.pyplot as plt

# Função de Langevin
# magnetização anisterética

# parâmetros: Ms, a e alpha
# magnetizacao de saturacao
# a  - relacionado à temperatura
# alpha - domínio magnéticos
# Man = Ms (coth((H+aM)/a)-a(H+alphaM))


f = 60
t_i = 0
t_f = 5/f # 5 periodos
deltat=0.00001
t = np.linspace(t_i,t_f,(t_f-t_i)/deltat)
n = round(t_f/deltat)

Ms = 1.5e6
alpha = 1e-4
a = 50
mi0 = 4e-7*np.pi

Bm = 1.2
B = Bm*np.sin(2*np.pi*f*t)

H = np.zeros([n])
He = np.zeros([n])
M = np.zeros([n])
dManHe = 0

for j in range(0,n-1):

	if abs(He[j]) <= 1e-3:
		dManHe = Ms/(3*a)
	else:
		dManHe = (Ms/a)*(1-(1/np.tanh(He[j]/a))**2+(a/He[j])**2)

	M[j+1] = (B[j+1]/mi0) - H[j+1]
	deltM = M[j+1] - M[j]
	H[j+1] = deltM*((1/dManHe)-alpha) + H[j]
	He[j+1] = H[j+1] + alpha*M[j+1]
	

fig1 = plt.figure(1)
plt.subplot(211)
plt.plot(t,B)
plt.grid(True)
plt.subplot(212)
plt.plot(t,H)
plt.grid(True)


fig2 = plt.figure(2)
plt.title('BH')
plt.plot(H,B)

plt.grid(True)
plt.show()