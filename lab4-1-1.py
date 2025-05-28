import numpy as np
import matplotlib.pyplot as plt

def exact_solution(x):
    return (1 + x) * np.exp(-x**2)

def f(x, y1, y2):
    return -4 * x * y2 - (4 * x**2 + 2) * y1

def euler_method(h, x_end):
    x = np.arange(0, x_end + h, h)
    y1 = np.zeros_like(x)
    y2 = np.zeros_like(x)
    y1[0] = 1
    y2[0] = 1
    for i in range(1, len(x)):
        y1[i] = y1[i-1] + h * y2[i-1]
        y2[i] = y2[i-1] + h * f(x[i-1], y1[i-1], y2[i-1])
    return x, y1

def implicit_euler_method(h, x_end):
    x = np.arange(0, x_end + h, h)
    y1 = np.zeros_like(x)
    y2 = np.zeros_like(x)
    y1[0] = 1
    y2[0] = 1

    for i in range(1, len(x)):
        x_next = x[i]
        # y1_next = y1[i-1] + h * y2_next
        # y2_next = y2[i-1] + h * (-4 * x_next * y2_next - (4 * x_next**2 + 2) * y1_next)
        # y2_next = y2[i-1] + h * (-4 * x_next * y2_next - (4 * x_next**2 + 2) * (y1[i-1] + h * y2_next))
        # y2_next * (1 + 4 * h * x_next + (4 * x_next**2 + 2) * h^2) = y2[i-1] - h * (4 * x_next**2 + 2) * y1[i-1]

        denominator = 1 + 4 * h * x_next + (4 * x_next**2 + 2) * h**2
        numerator = y2[i-1] - h * (4 * x_next**2 + 2) * y1[i-1]

        y2_next = numerator / denominator
        y1_next = y1[i-1] + h * y2_next

        y1[i] = y1_next
        y2[i] = y2_next

    return x, y1

def runge_kutta_method(h, x_end):
    x = np.arange(0, x_end + h, h)
    y1 = np.zeros_like(x)
    y2 = np.zeros_like(x)
    y1[0] = 1
    y2[0] = 1
    for i in range(1, len(x)):
        k1_y1 = h * y2[i-1]
        k1_y2 = h * f(x[i-1], y1[i-1], y2[i-1])

        k2_y1 = h * (y2[i-1] + 0.5 * k1_y2)
        k2_y2 = h * f(x[i-1] + 0.5 * h, y1[i-1] + 0.5 * k1_y1, y2[i-1] + 0.5 * k1_y2)

        k3_y1 = h * (y2[i-1] + 0.5 * k2_y2)
        k3_y2 = h * f(x[i-1] + 0.5 * h, y1[i-1] + 0.5 * k2_y1, y2[i-1] + 0.5 * k2_y2)

        k4_y1 = h * (y2[i-1] + k3_y2)
        k4_y2 = h * f(x[i-1] + h, y1[i-1] + k3_y1, y2[i-1] + k3_y2)

        y1[i] = y1[i-1] + (k1_y1 + 2 * k2_y1 + 2 * k3_y1 + k4_y1) / 6
        y2[i] = y2[i-1] + (k1_y2 + 2 * k2_y2 + 2 * k3_y2 + k4_y2) / 6

    return x, y1

def adams_method(h, x_end):
    x = np.arange(0, x_end + h, h)
    if len(x) < 4:
        raise ValueError("Шаг слишком большой для метода Адамса 4-го порядка")
    x_rk, y1_rk = runge_kutta_method(h, 3 * h)
    y1 = np.zeros_like(x)
    y2 = np.zeros_like(x)
    y1[:4] = y1_rk[:4]
    for i in range(4):
        y2[i] = (y1[i+1] - y1[i]) / h if i < 3 else 0
    for i in range(4):
        y2[i] = (y1[i+1] - y1[i]) / h if i < 3 else y2[i]
    for i in range(4, len(x)):
        y1_pred = y1[i-1] + h * (55 * y2[i-1] - 59 * y2[i-2] + 37 * y2[i-3] - 9 * y2[i-4]) / 24
        y2_pred = y2[i-1] + h * (55 * f(x[i-1], y1[i-1], y2[i-1]) -
                                  59 * f(x[i-2], y1[i-2], y2[i-2]) +
                                  37 * f(x[i-3], y1[i-3], y2[i-3]) -
                                  9 * f(x[i-4], y1[i-4], y2[i-4])) / 24

        y1[i] = y1[i-1] + h * (9 * y2_pred + 19 * y2[i-1] - 5 * y2[i-2] + y2[i-3]) / 24
        y2[i] = y2[i-1] + h * (9 * f(x[i], y1_pred, y2_pred) +
                                19 * f(x[i-1], y1[i-1], y2[i-1]) -
                                5 * f(x[i-2], y1[i-2], y2[i-2]) +
                                f(x[i-3], y1[i-3], y2[i-3])) / 24

    return x, y1

h = 0.1
x_end = 1
x_euler, y_euler = euler_method(h, x_end)
x_implicit_euler, y_implicit_euler = implicit_euler_method(h, x_end)
x_rk, y_rk = runge_kutta_method(h, x_end)
x_adams, y_adams = adams_method(h, x_end)

x_exact = np.linspace(0, x_end, 100)
y_exact = exact_solution(x_exact)

plt.figure(figsize=(10, 6))
plt.plot(x_exact, y_exact, 'k-', label='Точное решение')
plt.plot(x_euler, y_euler, 'bo--', label='Явный метод Эйлера')
plt.plot(x_implicit_euler, y_implicit_euler, 'co--', label='Неявный метод Эйлера')
plt.plot(x_rk, y_rk, 'go--', label='Метод Рунге-Кутты 4-го порядка')
plt.plot(x_adams, y_adams, 'ro--', label='Метод Адамса 4-го порядка')
plt.xlabel('x')
plt.ylabel('y(x)')
plt.legend()
plt.grid()
plt.title('Сравнение численных методов решения ОДУ')
plt.show()

def runge_romberg(h1, h2, y_h1, y_h2, p=4):
    error = np.abs(y_h1[-1] - y_h2[-1 * int(h2/h1)]) / (2**p - 1)
    return error

h2 = h / 2
x_rk_h, y_rk_h = runge_kutta_method(h, x_end)
x_rk_h2, y_rk_h2 = runge_kutta_method(h2, x_end)

rr_error = runge_romberg(h, h2, y_rk_h, y_rk_h2)
print(f"Оценка погрешности методом Рунге-Ромберга: {rr_error}")

y_exact_at_points = exact_solution(x_rk)
error_exact = np.abs(y_rk - y_exact_at_points)
print("Погрешность по сравнению с точным решением:")
for x_val, err in zip(x_rk, error_exact):
    print(f"x = {x_val:.1f}, погрешность = {err:.6f}")
