import numpy as np
import matplotlib.pyplot as plt

def exact_solution(x):
    return -np.tan(x)

def shooting_method(f, y0, yp0_guess, x0, xn, n):
    h = (xn - x0) / n
    x_values = np.linspace(x0, xn, n + 1)

    def rk4_system(y0, yp0, x_values):
        y = np.zeros((n + 1, 2))
        y[0] = [y0, yp0]
        for i in range(n):
            k1 = f(x_values[i], y[i])
            k2 = f(x_values[i] + h / 2, y[i] + h / 2 * k1)
            k3 = f(x_values[i] + h / 2, y[i] + h / 2 * k2)
            k4 = f(x_values[i] + h, y[i] + h * k3)
            y[i + 1] = y[i] + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        return y

    yp0 = yp0_guess
    while True:
        y = rk4_system(y0, yp0, x_values)
        if np.isclose(y[-1, 0], -np.sqrt(3) / 3, atol=1e-6):
            break
        yp0 += (y[-1, 0] - (-np.sqrt(3) / 3)) * 0.1
    return y[:, 0], x_values

def finite_difference_method(y0, yn, x0, xn, n):
    h = (xn - x0) / n
    x_values = np.linspace(x0, xn, n + 1)
    y = np.zeros(n + 1)
    y[0] = y0
    y[-1] = yn

    A = np.zeros((n - 1, n - 1))
    b = np.zeros(n - 1)

    for i in range(1, n):
        A[i - 1][i - 1] = -2 / h**2  # Главная диагональ
        if i > 1:
            A[i - 1][i - 2] = 1 / h**2  # Под диагональ
        if i < n - 1:
            A[i - 1][i] = 1 / h**2  # Над диагональю
        b[i - 1] = 0  # Правая часть

    b[0] -= y[0] / h**2  # Условия на границе 1
    b[-1] -= y[-1] / h**2  # Условия на границе 2

    y_inner = np.linalg.solve(A, b)
    y[1:n] = y_inner
    return y, x_values

def main():
    x0 = 0
    xn = np.pi / 6
    y0 = 0
    yn = -np.sqrt(3) / 3
    n = 10

    yp0_guess = -1
    y_shooting, x_values = shooting_method(lambda x, y: np.array([y[1], 2 * (1 + np.tan(x)**2) * y[0]]), y0, yp0_guess, x0, xn, n)
    y_fd, x_values_fd = finite_difference_method(y0, yn, x0, xn, n)
    y_exact = exact_solution(x_values)

    error_shooting = np.abs(y_shooting - y_exact)
    error_fd = np.abs(y_fd - y_exact)

    print("x\tExact y\t\tShooting y\tFDM y")
    for i in range(n + 1):
        print(f"{x_values[i]:.4f}\t{y_exact[i]:.6f}\t{y_shooting[i]:.6f}\t{y_fd[i]:.6f}")

    plt.plot(x_values, y_exact, label='Exact Solution', color='black', linestyle='--')
    plt.plot(x_values, y_shooting, label='Shooting Method', marker='o', markersize=4)
    plt.plot(x_values_fd, y_fd, label='Finite Difference Method', marker='x', markersize=4)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.title('Результаты численного решения')
    plt.grid()
    plt.show()

if __name__ == "__main__":
    main()
