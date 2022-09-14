import math
import matplotlib.pyplot as plt
import random
import numpy as np
from scipy.optimize import least_squares
from scipy.optimize import minimize


def f1(x):
    return x ** 3


def f2(x):
    return math.fabs(x - 0.2)


def f3(x):
    return x * math.sin(1 / x)


e = 0.001
a1 = 0
a2 = 0.1
b = 1

x = a1
i = 0
Fmin = f1(b)
xMin = 0

n1 = int((b - a1) / e - 1)
d1 = a1 + i * (b - a1) / (n1 + 1)

while i < n1:
    d1 = a1 + i * (b - a1) / (n1 + 1)
    x = a1 + d1
    if Fmin < f1(x):
        xMin = x
    i += 1
print("f = x ** 3")
print("number of 𝑓-calculations is -", n1)
print("Minimal value of function on X =", xMin, " is ", f1(xMin))

x = a1
i = 0
Fmin = f2(b)
xMin = 0

n1 = int((b - a1) / e - 1)
d1 = a1 + i * (b - a1) / (n1 + 1)

while i < n1:
    d1 = a1 + i * (b - a1) / (n1 + 1)
    x = a1 + d1
    if Fmin < f2(x):
        xMin = x
    i += 1
print("f = fabs(x - 0.2)")
print("number of 𝑓-calculations is -", n1)
print("Minimal value of function on X =", xMin, " is ", f2(xMin))

E = 0.001
a = 0.01
x = a
xMin = x
i = 0
n = int((b - a) / E - 1)

while i < n:
    d = a + i * (b - a) / (n + 1)
    x = a + d
    F1 = x * math.sin(1 / x)
    F2 = xMin * math.sin(1 / xMin)
    if F2 > F1:
        xMin = x
    i += 1
print("f = x * sin(1 / x)")
print("number of 𝑓-calculations is -", n)
print("Minimal value of function on X =", xMin, " is ", F2)


# dihotomiya
def MPD(f, a, b, eps=0.001):
    n = 0
    while abs(b - a) > eps:
        n += 1
        x = (a + b) / 2
        if (f(x + eps) > f(x - eps)):
            b = x
        else:
            a = x
    return x, n


x, n = MPD(f1, 0, 1)
print("f = x ** 3")
print("number of 𝑓-calculations is -", n)
print("X =", round(x, 5))
print("Minimal value of function", round(f1(x), 5))

x, n = MPD(f2, 0, 1)
print("f = fabs(x - 0.2)")
print("number of 𝑓-calculations is -", n)
print("X =", round(x, 5))
print("Minimal value of function", round(f2(x), 5))

x, n = MPD(f3, 0.1, 1)
print("f = x * sin(1 / x)")
print("number of 𝑓-calculations is -", n)
print("X =", round(x, 5))
print("Minimal value of function", round(f3(x), 5))


def goldenSectionSearch(f, l, r, eps):
    x1 = l + 0.382 * (r - l)
    x2 = r - 0.382 * (r - l)
    f1 = f(x1)
    f2 = f(x2)
    n = 0
    while (abs(r - l) > eps):
        n += 1
        if f1 < f2:
            r = x2
            x2 = x1
            f2 = f1
            x1 = l + 0.382 * (r - l)
            f1 = f(x1)
        else:
            l = x1
            x1 = x2
            f1 = f2
            x2 = r - 0.382 * (r - l)
            f2 = f(x2)
    return (x1 + x2) / 2


x = goldenSectionSearch(f1, 0, 1, 0.001)
print('X =', x)
print('f =', f1(x))
print('number of 𝑓-calculations', n)

x = goldenSectionSearch(f2, 0, 1, 0.001)
print('X =', x)
print('f =', f2(x))
print('number of 𝑓-calculations', n)

x = goldenSectionSearch(f3, 0.1, 1, 0.001)
print('X =', x)
print('f =', f3(x))
print('number of 𝑓-calculations', n)

massiveY = []
massiveK = list(range(1, 101))

for i in range(1, 101):
    x = i / 100
    a = random.random()
    b = random.random()
    y = a * x + b
    massiveY.append(y)
    y1 = i


def lin(a, x, y):
    return a[0] + a[1] * x - y


def fun(a, x, y):
    for i in range(0, len(x)):
        D = ((a[0] + a[1] * x[i] - y[i]) ** 2)
    return D


def fun_f(a, b, x, y):
    for i in range(0, len(x)):
        D = ((a + b * x[i] - y[i]) ** 2)
    return D


def exhaustive_search(f, x, y):
    a = np.array(np.linspace(0, 1, 100))
    b = np.array(np.linspace(0, 1, 100))
    min = f([0.5, 0], x, y)
    for i in range(len(a)):
        for j in range(len(b)):
            D = f([a[i], b[j]], x, y)
            if min > D:
                min = D
                mass_min = [round(a[i], 2), round(b[j], 2)]
    return mass_min


def MPD1(f, constant, a, b, x, y, eps=0.001):
    n = 0
    while abs(b - a) > eps:
        n += 1
        mid = (a + b) / 2
        if (f(mid + eps, constant, x, y) > f(mid - eps, constant, x, y)):
            b = mid
        else:
            a = mid
    return mid


def MPD2(f, constant, a, b, x, y, eps=0.001):
    n = 0
    while abs(b - a) > eps:
        n += 1
        mid = (a + b) / 2
        if (f(constant, mid + eps, x, y) > f(constant, mid - eps, x, y)):
            b = mid
        else:
            a = mid
    return mid


def gauss(f, a0, x, y):
    a = np.array([0.32, 0])
    n = 0
    while abs(f(a0[0], a0[1], x, y) - f(a[0], a[1], x, y)) >= 0.001:
        n += 1
        a = a0
        a0 = [a0[0] + 0.05, a0[1] + 0.05]
        for i in range(len(a0)):
            if i == 0:
                constant = a0[1]
                a0[0] = MPD1(f, constant, 0, 1, x, y)
            if i == 1:
                constant = a0[0]
                a0[1] = MPD2(f, constant, 0, 1, x, y)
    return a0


x = np.array(massiveK)
y = np.array(massiveY)
a0 = np.array([0.4, 0])

res_lsq = least_squares(lin, x0=a0, args=(x, y))
res_nm = minimize(fun, a0, args=(x, y), method='nelder-mead')
res_es = exhaustive_search(fun, x, y)
res_gauss = gauss(fun_f, a0, x, y)

print("a = %.2f, b = %.2f" % tuple(res_lsq.x))
print("a = %.2f, b = %.2f" % tuple(res_nm.x))
print("a = ", res_es[0], "b = ", res_es[1])
print("a = ", round(res_gauss[0], 2), "b = ", round(res_gauss[1], 2))
print()
f = lambda x: sum([u * v for u, v in zip(res_lsq.x, [1, x])])
f_nm = lambda x: sum([u * v for u, v in zip(res_nm.x, [1, x])])
f_es = lambda x: sum([u * v for u, v in zip(res_es, [1, x])])
f_gauss = lambda x: sum([u * v for u, v in zip(res_gauss, [1, x])])
x_p = np.linspace(min(x), max(x), 20)
y_p = f(x_p)
y_nm = f_nm(x_p)
y_es = f_es(x_p)
y_gauss = f_gauss(x_p)

fig, ax = plt.subplots()
plt.plot(x, y, '.r', label='Generated datas')
plt.plot(x_p, y_p, '-g', label='Least squares method')
plt.plot(x_p, y_nm, '-b', label='nelder-mead')
plt.plot(x_p, y_es, '-y', label='exhaustive_search')
plt.plot(x_p, y_gauss, color='black', label='gauss')
ax.legend()
plt.show()


def rational(a, x, y):
    return a[0] / (1 + a[1] * x) - y


def fun1(a, x, y):
    for i in range(0, len(x)):
        D = ((a[0] / (1 + a[1] * x[i]) - y[i]) ** 2)
    return D


def fun_f1(a, b, x, y):
    for i in range(0, len(x)):
        D = ((a / (1 + b * x[i]) - y[i]) ** 2)
    return D


x1 = np.array(massiveK)
y1 = np.array(massiveY)
a1 = np.array([0.45, 0])

res_lsq1 = least_squares(rational, x0=a1, args=(x1, y1))
res_nm1 = minimize(fun1, a1, args=(x1, y1), method='nelder-mead')
res_es1 = exhaustive_search(fun1, x1, y1)
res_gauss1 = gauss(fun_f1, a1, x1, y1)

print("a = %.2f, b = %.2f" % tuple(res_lsq1.x))
print("a = %.2f, b = %.2f" % tuple(res_nm1.x))
print("a = ", res_es1[0], "b = ", res_es1[1])
print("a = ", round(res_gauss1[0], 2), "b = ", round(res_gauss1[1], 2))
print()
f1 = lambda x1: sum([u * v for u, v in zip(res_lsq1.x, [1, x1])])
f_nm1 = lambda x1: sum([u * v for u, v in zip(res_nm1.x, [1, x1])])
f_es1 = lambda x: sum([u * v for u, v in zip(res_es1, [1, x])])
f_gauss1 = lambda x: sum([u * v for u, v in zip(res_gauss1, [1, x])])
x_p1 = np.linspace(min(x1), max(x1), 20)
y_p1 = f(x_p1)
y_nm1 = f_nm1(x_p1)
y_es1 = f_es1(x_p1)
y_gauss1 = f_gauss1(x_p1)

fig, ax = plt.subplots()
plt.plot(x, y, '.r', label='Generated datas')
plt.plot(x_p1, y_p1, '-g', label='Least squares method')
plt.plot(x_p1, y_nm1, '-b', label='nelder-mead')
plt.plot(x_p1, y_es1, '-y', label='exhaustive_search')
plt.plot(x_p1, y_gauss1, color='black', label='gauss')
ax.legend()
plt.show()
