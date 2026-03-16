import numpy as np
import matplotlib.pyplot as plt

x = [50,100,200,400,800]
y = [20,35,60,110,210]

def divided_differences(x, y):
    n = len(y)
    coef = np.zeros([n, n])
    coef[:,0] = y

    for j in range(1,n):
        for i in range(n-j):
            coef[i][j] = (coef[i+1][j-1] - coef[i][j-1])/(x[i+j] - x[i])

    return coef

#ф-цію обч знач інтерполяційного многочлена н
def newton_polynomial(x_data, coef, x):
    n = len(x_data) - 1
    p = coef[0][n]
#будує многочлен н
    for k in range(1,n+1):
        p = coef[0][n-k] + (x-x_data[n-k])*p

    return p

#виклик ф-ції і будується табл р р
coef = divided_differences(x,y)

cpu = newton_polynomial(x,coef,600)

print("CPU при 600 RPS =", cpu)

x_plot = np.linspace(50,800,100)
y_plot = [newton_polynomial(x,coef,i) for i in x_plot]
#експериментальні т і крива н
plt.scatter(x,y)
plt.plot(x_plot,y_plot)

plt.xlabel("RPS")
plt.ylabel("CPU")

plt.grid()

plt.show()
coef = divided_differences(x,y)
print("Таблиця розділених різниць:")
print(coef)

# прогноз для 600
cpu = newton_polynomial(x,coef,600)
print("CPU при 600 RPS =", cpu)

# дослідження вузлів
x_small = x[:3]
y_small = y[:3]

coef_small = divided_differences(x_small,y_small)
cpu_small = newton_polynomial(x_small,coef_small,600)

print("CPU (3 вузли) =", cpu_small)

x_medium = x[:4]
y_medium = y[:4]

coef_medium = divided_differences(x_medium,y_medium)
cpu_medium = newton_polynomial(x_medium,coef_medium,600)

print("CPU (4 вузли) =", cpu_medium)