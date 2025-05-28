import numpy as np
from math import cos, sin, log10

def f1(x1, x2):
    return x1 - cos(x2) - 1

def f2(x1, x2):
    return x2 - log10(x1 + 1) - 3

def jacobian(x1, x2):
    ln10 = np.log(10)
    df1dx1 = 1
    df1dx2 = sin(x2)
    df2dx1 = -1 / ((x1 + 1) * ln10)
    df2dx2 = 1
    return np.array([
        [df1dx1, df1dx2],
        [df2dx1, df2dx2]
    ])

x = np.array([1.0, 1.0])
tolerance = 1e-6
max_iter = 50

print("iter\t   x1\t\t   x2\t\t   f1\t\t   f2")
print("-" * 60)

for i in range(max_iter):
    F = np.array([f1(x[0], x[1]), f2(x[0], x[1])])
    norm = np.linalg.norm(F)
    if norm < tolerance:
        print(f"\nРешение найдено за {i} итераций:")
        print(f"x1 = {x[0]:.7f}, x2 = {x[1]:.7f}")
        break

    print(f"{i:2d}\t{x[0]:12.7f}\t{x[1]:12.7f}\t{F[0]:12.7f}\t{F[1]:12.7f}")
    J = jacobian(x[0], x[1])
    delta = np.linalg.solve(J, -F)
    x += delta

else:
    print("\nМаксимальное число итераций достигнуто")
