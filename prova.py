import math
import numpy as np
import matplotlib.pyplot as plt

# frequencia e tempo
f = 60
t_i = 0
t_f = 12/f
deltat=1e-6 #passo
t = np.linspace(t_i,t_f,(t_f-t_i)/deltat)
n = round(t_f/deltat)

# parâmetros
N1 = 212
N2 = 110
R1 = 0.5
R2 = 0.3

r_ext = 0.15
r_int = 0.10
rmedio = (r_ext-r_int)/(math.log(r_ext/r_int))
lm = 2*np.pi*rmedio
nchapas = 20
d = 3e-3 #espessura
S = (nchapas*d)*(r_ext-r_int)

densidade = 7900 #kg/m3
resistividade = 30e-8 #ohms.m
sigma = 1/resistividade
mi0 = 4e-7*np.pi

vol = S * lm

kex = 0.53
Ms = 1.25e6
k = 720
c = 350e-3
alpha = 220.05e-5
a = 900

# tensao senoidal
Vp = 25
V = Vp*np.sin(2*np.pi*f*t)

# condições iniciais
I = np.zeros([n])
B = np.zeros([n])

L = np.zeros([n])
phi = np.zeros([n])

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

# fluxo inicial, remanencia
B[0] = 0.4
M[0] = B[0]/mi0

for j in range(0, n-1):

	dBdt = (1/(N1*S))*(V[j] - R1*I[j])
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
	I[j+1] = Htotal[j+1] * (lm/N1)

	# phi[j+1] = B[j+1] * S
	# L[j+1] = ((N1+N2)*phi[j+1])/I[j+1]

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

perdas_totais = area_total * 300 * vol
print('Perdas totais =', perdas_totais, 'W')
perdas_hist = area_hist * 300 * vol
print('Perdas por histerese =', perdas_hist, 'W')
perdas_ind = area_ind * 300 * vol
print('Perdas por correntes induzidas =', perdas_ind, 'W')
perdas_ex = area_ex * 300 * vol
print('Perdas Excedentes =', perdas_ex, 'W')

# 1) Corrente de pico e eficaz

# I máx
print("I pico: ", max(I), 'A')

I_int = 0
Isq = I ** 2
soma = 0
for j in range(ni,nf-1):
	soma = deltat/2 * (Isq[j] + Isq[j+1])
	I_int = I_int + soma
Ief = math.sqrt(f * I_int)
print("Ief: ", Ief)

fig1 = plt.figure(1)
plt.plot(t,I)
plt.xlabel('tempo (s)')
plt.ylabel('corrente (A)')
plt.title('Corrente no primário')

# 2) Tensão induzida no secundário

Vind = V*(N2/N1)
# Vs máx
print("Vs pico: ", max(Vind), 'V')

fig2 = plt.figure(2)
plt.plot(t,Vind)
plt.xlabel('tempo (s)')
plt.ylabel('tensão induzida (V)')
plt.title('Tensão induzida no secundário')


# 3) Perdas na resistência e no núcleo

# perda na resistência
potencia_resistor = R1 * Ief**2
print('Potencia resistor', potencia_resistor)

# curvas BH

fig4 = plt.figure(4)
plt.xlabel('H (A/m)')
plt.ylabel('B (T)')
plt.title('Curva BH')
plt.plot(Htotal,B, color='red')


fig3 = plt.figure(3)
plt.xlabel('H (A/m)')
plt.ylabel('B (T)')
plt.title('Curvas BH')
plt.plot(Htotal,B, label='Perdas totais')
plt.plot(Hhist,B, color='green', label='Histerese')
plt.legend(loc='upper left')

fig5 = plt.figure(5)
plt.subplot(221)
plt.title('Total')
plt.plot(Htotal,B, color='green')
plt.subplot(222)
plt.title('Histerese')
plt.plot(Hhist,B, color='green')
plt.subplot(223)
plt.title('Correntes induzidas')
plt.plot(Hind,B, color='green')
plt.subplot(224)
plt.title('Excedentes')
plt.plot(Hex,B, color='green')

plt.tight_layout()

# 4) triangulo de potenciais e fp

# tensão eficaz
V_int = 0
Vsq = V ** 2	
for i in range(ni,nf-1):
	V_int = V_int + (deltat/2 * (Vsq[i] + Vsq[i+1]))
Vef = math.sqrt(f * V_int)
print("Vef: ", Vef)


# potencia aparente
P_aparente = 0
for j in range(ni,nf-1):
	P_aparente = Vef*Ief
print("Potência aparente: ", P_aparente)

#perdas totais
potencia_total = 0
P = V * I
for j in range(ni,nf-1):
	potencia_total = potencia_total + ((deltat/2) * (P[j] + P[j+1]))
potencia_total_circ = potencia_total * f
print("Potencia total no circuito: ", potencia_total_circ)

#fp
fp = potencia_total_circ/P_aparente
print("fp: ", fp)


plt.show()