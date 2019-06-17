import math
import numpy as np
import matplotlib.pyplot as plt

# Circuito magnético
# -V + Ri + Ndphidt = 0
# phi = int Bds
# nesse caso, S constante
# -V + Ri + NSdBdt = 0

# usando JA pra permeabilidade
# acrescentando componentes de perdas por
# correntes induzidas e excedentes

# parcela de campo correntes induzidas
# Hind = (sigma*d**2/12)*dBdt
# parcela de campo correntes induzidas
# Hex = (sqrt(kex)*1/(sqrt(deltat)*abs()dB))

# perda no resistor, perda no núcleo (com separação)
# fp da estrutura

# testar com V sin -> não tá funcionando
# fluxo inicial
# B0 = 0,5
# M0 = B0/mu0

lm = 0.94
S = 105e-6
N = 700
R = 1.5
mi0 = 4e-7*np.pi
mi = 1000*mi0
Vp = 30
f = 60
t_i = 0
t_f = 12/f
deltat=1e-6 #passo
t = np.linspace(t_i,t_f,(t_f-t_i)/deltat)
n = round(t_f/deltat)

sigma = 2.5e6
d = 1e-3

kex = 0.39

Ms = 1.5779e6
k = 5.7311e1
c = 2.6972e-1
alpha = 1.9975e-4
a = 1.0498e2

vol = 9.87e-5

# tensao senoidal
V = Vp*np.sin(2*np.pi*f*t)
# condições iniciais
I = np.zeros([n])
B = np.zeros([n])

Htotal = np.zeros([n])
Hhist = np.zeros([n])
Hind = np.zeros([n])
Hex = np.zeros([n])

# JA
He = np.zeros([n])
M = np.zeros([n])
Man = np.zeros([n])
Mirr = np.zeros([n])
dMandHe = 0
deltB = 0
dMirrdBe = 0
dMdB = 0

# fluxo inicial
B[0] = -0.5
#B[0] = 0.5
M[0] = B[0]/mi0


for j in range(0, n-1):
	dBdt = (1/(N*S))*(V[j] - R*I[j])
	B[j+1] = B[j] + dBdt*deltat
	deltB = B[j+1] - B[j]
	
	# pra V sin
	if (deltB == 0):
		deltB = 1e-8

	if deltB > 0:
		delta = 1
	else:
		delta = -1

	M[j] = (B[j]/mi0) - Hhist[j]
	He[j] = Hhist[j] + alpha*M[j]

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
	Hhist[j+1] = (B[j+1]/mi0) - M[j+1]


	Hind[j+1] = ((sigma*(d**2))/12) * dBdt
	Hex[j+1] = kex * (1/(math.sqrt(deltat)*math.sqrt(abs(deltB)))) * deltB
		
	Htotal[j+1] = Hhist[j+1] + Hind[j+1] + Hex[j+1]
	I[j+1] = Htotal[j+1] * (lm/N)

# cálculo das perdas
nf = round(t_f/deltat)
ni = round((t_f-(1/f))/deltat)

# perdas totais
area_total = 0
for k in range(ni,nf-1):
	area_total = area_total + (((Htotal[k+1]+Htotal[k])/2)*(B[k+1]-B[k]))
print('Perdas Totais =', area_total, 'J')

# perdas por correntes induzidas
area_ind = 0
for k in range(ni,nf-1):
	area_ind = area_ind + (((Hind[k+1]+Hind[k])/2)*(B[k+1]-B[k]))
print('Perdas por correntes induzidas =', area_ind, 'J')

# perdas excedentes
area_ex = 0
for k in range(ni,nf-1):
	area_ex = area_ex + (((Hex[k+1]+Hex[k])/2)*(B[k+1]-B[k]))
print('Perdas Excedentes =', area_ex, 'J')

# perdas por histerese
area_hist = 0
for k in range(ni,nf-1):
	area_hist = area_hist + (((Hhist[k+1]+Hhist[k])/2)*(B[k+1]-B[k]))
print('Perdas por histerese =', area_hist, 'J')


# perdas no ferro em W (considerando o volume)
perdas_totais = area_total * f * vol
print('Perdas totais =', perdas_totais, 'W')
perdas_hist = area_hist * f * vol
print('Perdas por histerese =', perdas_hist, 'W')
perdas_ind = area_ind * f * vol
print('Perdas por correntes induzidas =', perdas_ind, 'W')
perdas_ex = area_ex * f * vol
print('Perdas Excedentes =', perdas_ex, 'W')

# corrente eficaz
I_int = 0
Isq = I ** 2
soma = 0
for j in range(ni,nf-1):
	soma = deltat/2 * (Isq[j] + Isq[j+1])
	I_int = I_int + soma
Ief = math.sqrt(f * I_int)
print("Ief: ", Ief)
# tensão eficaz
V_int = 0
Vsq = V ** 2	
for i in range(ni,nf-1):
	V_int = V_int + (deltat/2 * (Vsq[i] + Vsq[i+1]))
Vef = math.sqrt(f * V_int)
print("Vef: ", Vef)

# potencia aparente
S = 0
for j in range(ni,nf-1):
	S = Vef*Ief
print("Potência aparente: ", S)

# perdas no resistor
potencia_resistor = R * Ief**2
print('Potencia resistor', potencia_resistor)

# perdas totais
potencia_total = 0
P = V * I
for j in range(ni,nf-1):
	potencia_total = potencia_total + ((deltat/2) * (P[j] + P[j+1]))
potencia_total_circ = potencia_total * f
print("Potencia total no circuito: ", potencia_total_circ)

# fp
fp = potencia_total_circ/S
print("fp: ", fp)

# I máx
# print("Imáx: ", max(I), 'A')

fig1 = plt.figure(1)

plt.subplot(411)
plt.plot(t,V)
plt.title('V')
plt.grid(True)
plt.subplot(412)
plt.plot(t,I, color='xkcd:fuchsia')
plt.title('I')
plt.grid(True)

plt.subplot(413)
plt.plot(t,B)
plt.title('B')
plt.grid(True)
plt.subplot(414)
plt.plot(t,Htotal, color='red')
plt.title('H')
plt.grid(True)
plt.tight_layout()

# fig2 = plt.figure(2)
# plt.subplot(221)
# plt.title('BHotal')
# plt.plot(Htotal,B, color='green')
# plt.subplot(222)
# plt.title('BHhist')
# plt.plot(Hhist,B, color='green')
# plt.subplot(223)
# plt.title('BHind')
# plt.plot(Hind,B, color='green')
# plt.subplot(224)
# plt.title('BHex')
# plt.plot(Hex,B, color='green')

# plt.tight_layout()

fig3 = plt.figure(3)
plt.title('BH')
plt.plot(Htotal,B)

# fig3 = plt.figure(3)
# plt.title('BH')
# plt.plot(Hhist,B, color='green')

# fig4 = plt.figure(4)
# plt.title('BH')
# plt.plot(Hind,B, color='green')

# fig5 = plt.figure(5)
# plt.title('BH')
# plt.plot(Hex,B, color='green')

plt.show()