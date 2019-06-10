import math
import numpy as np
import matplotlib.pyplot as plt

# Circuito magnético
# -V + Ri + Ndphidt = 0
# phi = int Bds
# nesse caso, S constante
# -V + Ri + NSdBdt = 0

# usando JA pra permeabilidade

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

Ms = 1.5779e6
k = 5.7311e1
c = 2.6972e-1
alpha = 1.9975e-4
a = 1.0498e2

# tensao senoidal
V = Vp*np.cos(2*np.pi*f*t)
# condições iniciais
I = np.zeros([n])
H = np.zeros([n])
B = np.zeros([n])

# JA
He = np.zeros([n])
M = np.zeros([n])
Man = np.zeros([n])
Mirr = np.zeros([n])
dMandHe = 0
deltB = 0
dMirrdBe = 0
dMdB = 0


for j in range(0, n-1):
	dBdt = (1/(N*S))*(V[j] - R*I[j])
	B[j+1] = B[j] + dBdt*deltat

	deltB = B[j+1] - B[j]

	if deltB > 0:
		delta = 1
	else:
		delta = -1

	M[j] = (B[j]/mi0) - H[j]
	He[j] = H[j] + alpha*M[j]

	if abs(He[j])/a <= 0.1:
		Man[j] = (Ms/(3*a))*He[j]
		dMandHe = Ms/(3*a)
	else:
		Man[j] = Ms*((1/np.tanh(He[j]/a))-(a/He[j]))
		dMandHe = (Ms/a)*(1-(1/np.tanh(He[j]/a))**2+(a/He[j])**2)

	Mirr[j] = (M[j] - c*Man[j])/(1 - c)
	dMirrdBe = (1/mi0)*((Man[j]-Mirr[j])/(k*delta))

	if dMirrdBe < 0:
		dMirrdBe = 0

	dMdB = ((1-c)*dMirrdBe+((c/mi0)*dMandHe))/(1+(mi0*(1-alpha)*(1-c)*dMirrdBe)+(c*(1-alpha)*dMandHe))
	
	M[j+1] = M[j]+ (dMdB*deltB)
	H[j+1] = (B[j+1]/mi0) - M[j+1]

	I[j+1] = H[j+1] * (lm/N)


fig1 = plt.figure(1)

plt.subplot(221)
plt.plot(t,V)
plt.title('V')
plt.grid(True)
plt.subplot(222)
plt.plot(t,I, color='xkcd:fuchsia')
plt.title('I')
plt.grid(True)

plt.subplot(223)
plt.plot(t,B)
plt.title('B')
plt.grid(True)
plt.subplot(224)
plt.plot(t,H, color='red')
plt.title('H')
plt.grid(True)
plt.tight_layout()

fig2 = plt.figure(2)
plt.title('BH')
plt.plot(H,B, color='green')

plt.show()