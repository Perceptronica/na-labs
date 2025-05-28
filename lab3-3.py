import numpy as np
import matplotlib.pyplot as plt

x = np.array([-3.0, -2.0, -1.0, 0.0, 1.0, 2.0])
y = np.array([0.04979, 0.13534, 0.36788, 1.0, 2.7183, 7.3891])
#x = np.array([0.0, 1.7, 3.4, 5.1, 6.8, 8.5])
#y = np.array([0.0,1.3038, 1.8439, 2.2583, 2.6077, 2.9155])
n = len(x)

sum_x = np.sum(x)
sum_x2 = np.sum(x**2)
sum_x3 = np.sum(x**3)
sum_x4 = np.sum(x**4)
sum_y = np.sum(y)
sum_xy = np.sum(x * y)
sum_x2y = np.sum(x**2 * y)

A = np.array([[n, sum_x],
              [sum_x, sum_x2]])
b = np.array([sum_y, sum_xy])
a0, a1 = np.linalg.solve(A, b)

A = np.array([[n, sum_x, sum_x2],
              [sum_x, sum_x2, sum_x3],
              [sum_x2, sum_x3, sum_x4]])
b = np.array([sum_y, sum_xy, sum_x2y])
a0_q, a1_q, a2_q = np.linalg.solve(A, b)

def linear_poly(x):
    return a0 + a1 * x

def quadratic_poly(x):
    return a0_q + a1_q * x + a2_q * x**2

y_linear = linear_poly(x)
y_quadratic = quadratic_poly(x)

sse_linear = np.sum((y - y_linear)**2)
sse_quadratic = np.sum((y - y_quadratic)**2)

print("Линейная аппроксимация (1-я степень):")
print(f"y = {a0:.4f} + {a1:.4f} * x")
print(f"SSE = {sse_linear:.6f}\n")

print("Квадратичная аппроксимация (2-я степень):")
print(f"y = {a0_q:.4f} + {a1_q:.4f} * x + {a2_q:.4f} * x^2")
print(f"SSE = {sse_quadratic:.6f}")

plt.figure(figsize=(10, 6))
plt.scatter(x, y, color='red', label='Исходные данные')
x_plot = np.linspace(-3, 2, 100)
plt.plot(x_plot, linear_poly(x_plot), label='Линейная аппроксимация')
plt.plot(x_plot, quadratic_poly(x_plot), label='Квадратичная аппроксимация')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid()
plt.title('Аппроксимация методом наименьших квадратов')
plt.show()
