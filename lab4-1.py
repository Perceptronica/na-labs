import numpy as np
import matplotlib.pyplot as plt

def f(x, y):
    return np.array([y[1], -4 * x * y[1] - (4 * x ** 2 + 2) * y[0]])

def exact_solution(x):
    return (1 + x) * np.exp(-x**2)

def euler_method(f, y0, x0, h, n):
    y = np.zeros((n+1, len(y0)))
    y[0] = y0
    x = x0
    for i in range(n):
        y[i + 1] = y[i] + h * f(x, y[i])
        x += h
    return y

def runge_kutta_method(f, y0, x0, h, n):
    y = np.zeros((n + 1, len(y0)))
    y[0] = y0
    x = x0
    for i in range(n):
        k1 = f(x, y[i])
        k2 = f(x + h / 2, y[i] + h / 2 * k1)
        k3 = f(x + h / 2, y[i] + h / 2 * k2)
        k4 = f(x + h, y[i] + h * k3)
        y[i + 1] = y[i] + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        x += h
    return y

def adams_method(f, y0, x0, h, n):
    y = np.zeros((n + 1, len(y0)))
    y[0] = y0
    for i in range(3):
        y[i + 1] = runge_kutta_method(f, y0, x0, h, i + 1)[i + 1]
    x = x0 + h * 3
    for i in range(3, n):
        y[i + 1] = (y[i] + (h / 24) * (9 * f(x, y[i]) - 19 * f(x - h, y[i - 1]) + \
                                         5 * f(x - 2 * h, y[i - 2]) - f(x - 3 * h, y[i - 3])))
        x += h
    return y

def main():
    x0 = 0
    y0 = np.array([1, 1])  # y(0) = 1, y'(0) = 1
    h = 0.1
    n = int(1 / h)

    y_euler = euler_method(f, y0, x0, h, n)
    y_rk = runge_kutta_method(f, y0, x0, h, n)
    y_adams = adams_method(f, y0, x0, h, n)

    x_values = np.arange(x0, x0 + (n + 1) * h, h)
    y_exact = exact_solution(x_values)

    error_euler = np.abs(y_euler[:, 0] - y_exact)
    error_rk = np.abs(y_rk[:, 0] - y_exact)
    error_adams = np.abs(y_adams[:, 0] - y_exact)

    print("x\tExact y\tEuler y\tRK y\tAdams y")
    for i in range(n + 1):
        print(f"{x_values[i]:.2f}\t{y_exact[i]:.4f}\t{y_euler[i][0]:.4f}\t{y_rk[i][0]:.4f}\t{y_adams[i][0]:.4f}")

    plt.plot(x_values, y_exact, label='Exact Solution', color='black', linestyle='--')
    plt.plot(x_values, y_euler[:, 0], label='Euler', marker='o', markersize=4)
    plt.plot(x_values, y_rk[:, 0], label='Runge-Kutta', marker='x', markersize=4)
    plt.plot(x_values, y_adams[:, 0], label='Adams', marker='s', markersize=4)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.title('Результаты численного решения')
    plt.grid()
    plt.show()

if __name__ == "__main__":
    main()
