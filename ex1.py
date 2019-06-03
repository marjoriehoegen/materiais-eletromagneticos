import numpy as np
import matplotlib.pyplot as plt

# Integral de uma funcao

x = np.arange(0,11)

def y(x):
	return -x**2 + 2*x + 100

plt.plot(x,y(x))
# plt.show()

# calculo da integral
h = 0.0001
area_total = 0
soma = 0

# numpy.arange(start,stop,step)
x_passo = np.arange(0,11,h)

for i in x_passo:
	soma = (h/2) * (y(i) + y(i+1))
	area_total = area_total + soma

print (area_total)