import math
import numpy as np
import matplotlib.pyplot as plt

# Jiles-Atherton

# parâmetros: Ms, k, c, a e alpha
# magnetização anisterética
# magnetização reversível
# taxa de variação da mag reversível
# magnetização total

f = 60
t_i = 0
t_f = 5/f # 5 periodos
deltat=0.00001
t = np.linspace(t_i,t_f,(t_f-t_i)/deltat)
n = round(t_f/deltat)

Ms = 1.5779e6
k = 5.7311e1
c = 2.6972e-1
alpha = 1.9975e-4
a = 1.0498e2
mi0 = 4e-7*np.pi

Bm = 1.2
B = Bm*np.sin(2*np.pi*f*t)

H = np.zeros([n])
He = np.zeros([n])
M = np.zeros([n])
Man = np.zeros([n])
Mirr = np.zeros([n])
dMandHe = 0
deltB = 0
dMirrdBe = 0
dMdB = 0

for j in range(0,n-1):

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


# calcular perda por histerese pra Bm = 1.2T
area = 0
nf = round(t_f/deltat)
ni = round((t_f-(1/f))/deltat)

for k in range(ni,nf-1):
	area = area + (((H[k+1]+H[k])/2)*(B[k+1]-B[k]))

print('Perdas =', area, 'J')
perdaHist = area * f
perdaHist_k = (area * f)/1000
print("Perdas = ", perdaHist_k, "kW")


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