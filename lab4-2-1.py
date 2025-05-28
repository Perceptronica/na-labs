import numpy as np
import matplotlib.pyplot as plt

# Точное решение
def exact_solution(x):
    return (1 + x) * np.exp(-x**2)

# Правая часть системы ОДУ: dy/dx = z, dz/dx = -4xz - (4x^2 + 2)y
def f(x, y, z):
    return z, -4 * x * z - (4 * x**2 + 2) * y

# Метод Эйлера
def euler_method(h, x_end):
    x = np.arange(0, x_end + h, h)
    y = np.zeros_like(x)
    z = np.zeros_like(x)
    y[0], z[0] = 1, 1  # начальные условия

    for i in range(len(x) - 1):
        dy, dz = f(x[i], y[i], z[i])
        y[i + 1] = y[i] + h * dy
        z[i + 1] = z[i] + h * dz

    return x, y

# Метод Рунге-Кутты 4-го порядка
def runge_kutta_4th_order(h, x_end):
    x = np.arange(0, x_end + h, h)
    y = np.zeros_like(x)
    z = np.zeros_like(x)
    y[0], z[0] = 1, 1  # начальные условия

    for i in range(len(x) - 1):
        k1_y, k1_z = f(x[i], y[i], z[i])
        k2_y, k2_z = f(x[i] + h/2, y[i] + h/2 * k1_y, z[i] + h/2 * k1_z)
        k3_y, k3_z = f(x[i] + h/2, y[i] + h/2 * k2_y, z[i] + h/2 * k2_z)
        k4_y, k4_z = f(x[i] + h, y[i] + h * k3_y, z[i] + h * k3_z)

        y[i + 1] = y[i] + h/6 * (k1_y + 2*k2_y + 2*k3_y + k4_y)
        z[i + 1] = z[i] + h/6 * (k1_z + 2*k2_z + 2*k3_z + k4_z)

    return x, y

# Метод Адамса 4-го порядка (использует Рунге-Кутты для начальных точек)
def adams_4th_order(h, x_end):
    x = np.arange(0, x_end + h, h)
    y = np.zeros_like(x)
    z = np.zeros_like(x)
    y[0], z[0] = 1, 1  # начальные условия

    # Используем Рунге-Кутты для первых 4 точек
    for i in range(3):
        k1_y, k1_z = f(x[i], y[i], z[i])
        k2_y, k2_z = f(x[i] + h/2, y[i] + h/2 * k1_y, z[i] + h/2 * k1_z)
        k3_y, k3_z = f(x[i] + h/2, y[i] + h/2 * k2_y, z[i] + h/2 * k2_z)
        k4_y, k4_z = f(x[i] + h, y[i] + h * k3_y, z[i] + h * k3_z)

        y[i + 1] = y[i] + h/6 * (k1_y + 2*k2_y + 2*k3_y + k4_y)
        z[i + 1] = z[i] + h/6 * (k1_z + 2*k2_z + 2*k3_z + k4_z)

    # Метод Адамса для остальных точек
    for i in range(3, len(x) - 1):
        # Предсказание (по Адамсу-Бэшфорту)
        f0_y, f0_z = f(x[i], y[i], z[i])
        f1_y, f1_z = f(x[i - 1], y[i - 1], z[i - 1])
        f2_y, f2_z = f(x[i - 2], y[i - 2], z[i - 2])
        f3_y, f3_z = f(x[i - 3], y[i - 3], z[i - 3])

        y_pred = y[i] + h/24 * (55*f0_y - 59*f1_y + 37*f2_y - 9*f3_y)
        z_pred = z[i] + h/24 * (55*f0_z - 59*f1_z + 37*f2_z - 9*f3_z)

        # Коррекция (по Адамсу-Моултону)
        f_pred_y, f_pred_z = f(x[i + 1], y_pred, z_pred)

        y[i + 1] = y[i] + h/24 * (9*f_pred_y + 19*f0_y - 5*f1_y + f2_y)
        z[i + 1] = z[i] + h/24 * (9*f_pred_z + 19*f0_z - 5*f1_z + f2_z)

    return x, y

# Оценка погрешности методом Рунге-Ромберга
def runge_romberg_error(y_h, y_h2, p=4):
    error = np.abs((y_h - y_h2[::2]) / (2**p - 1))
    return error

# Параметры
h = 0.1
x_end = 1.0

# Решение разными методами
x_euler, y_euler = euler_method(h, x_end)
x_rk, y_rk = runge_kutta_4th_order(h, x_end)
x_adams, y_adams = adams_4th_order(h, x_end)

# Точное решение
x_exact = np.linspace(0, x_end, 100)
y_exact = exact_solution(x_exact)

# Погрешность по сравнению с точным решением
y_exact_at_points = exact_solution(x_rk)
error_euler = np.abs(y_euler - y_exact_at_points)
error_rk = np.abs(y_rk - y_exact_at_points)
error_adams = np.abs(y_adams - y_exact_at_points)

# Оценка погрешности методом Рунге-Ромберга (для Рунге-Кутты)
h2 = h / 2
x_rk_h2, y_rk_h2 = runge_kutta_4th_order(h2, x_end)
y_rk_h2_at_h_points = y_rk_h2[::2]  # Берем точки с шагом h
rr_error_rk = runge_romberg_error(y_rk, y_rk_h2)

# Вывод результатов
print("x\t Euler\t\t RK4\t\t Adams\t\t Exact\t\t Euler Error\t RK4 Error\t Adams Error\t RR Error (RK4)")
for i in range(len(x_rk)):
    print(f"{x_rk[i]:.1f}\t {y_euler[i]:.6f}\t {y_rk[i]:.6f}\t {y_adams[i]:.6f}\t {y_exact_at_points[i]:.6f}\t"
          f" {error_euler[i]:.6f}\t {error_rk[i]:.6f}\t {error_adams[i]:.6f}\t {rr_error_rk[i]:.6f}")

# Графики
plt.figure(figsize=(12, 6))
plt.plot(x_exact, y_exact, label='Точное решение', color='black')
plt.plot(x_euler, y_euler, 'o-', label='Метод Эйлера')
plt.plot(x_rk, y_rk, 's-', label='Метод Рунге-Кутты 4-го порядка')
plt.plot(x_adams, y_adams, '^-', label='Метод Адамса 4-го порядка')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid()
plt.title('Сравнение численных методов решения ОДУ')
plt.show()

# Графики погрешностей
plt.figure(figsize=(12, 6))
plt.plot(x_rk, error_euler, 'o-', label='Погрешность метода Эйлера')
plt.plot(x_rk, error_rk, 's-', label='Погрешность метода Рунге-Кутты')
plt.plot(x_rk, error_adams, '^-', label='Погрешность метода Адамса')
plt.plot(x_rk, rr_error_rk, 'x-', label='Оценка Рунге-Ромберга (для РК4)')
plt.xlabel('x')
plt.ylabel('Абсолютная погрешность')
plt.legend()
plt.grid()
plt.title('Погрешности численных методов')
plt.show()
