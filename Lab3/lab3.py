import csv
import matplotlib.pyplot as plt


def read_csv(filename):
    x = []
    y = []

    with open(filename, "r", encoding="utf-8") as file:
        reader = csv.reader(file)
        next(reader)

        for row in reader:
            x.append(float(row[0]))
            y.append(float(row[1]))

    return x, y


def form_matrix(x, m):
    #к-сть коеф полінома
    n = m + 1
    A = []

    for i in range(n):
        row = []
        for j in range(n):
            s = 0.0
            #по всіх вузлах табл
            for k in range(len(x)):
                s += x[k] ** (i + j)
            row.append(s)
        A.append(row)

    return A

    #права част сист
def form_vector(x, y, m):
    n = m + 1
    b = []

    for i in range(n):
        s = 0.0
        for k in range(len(x)):
            s += y[k] * (x[k] ** i)
        b.append(s)

    return b

#розв сист рівн.к-сть рівн у системі розмір матриці
def gauss_solve(A, b):
    n = len(A)

    a = [row[:] for row in A]
    rhs = b[:]

    for k in range(n):
        max_row = k
        max_value = abs(a[k][k])

        for i in range(k + 1, n):
            if abs(a[i][k]) > max_value:
                max_value = abs(a[i][k])
                max_row = i

        if abs(a[max_row][k]) < 1e-12:
            raise ValueError("Матриця вироджена або система не має єдиного розв'язку.")

        if max_row != k:
            a[k], a[max_row] = a[max_row], a[k]
            rhs[k], rhs[max_row] = rhs[max_row], rhs[k]

        for i in range(k + 1, n):
            factor = a[i][k] / a[k][k]
            for j in range(k, n):
                a[i][j] = a[i][j] - factor * a[k][j]
            rhs[i] = rhs[i] - factor * rhs[k]

    x_sol = [0.0 for _ in range(n)]

    for i in range(n - 1, -1, -1):
        s = 0.0
        for j in range(i + 1, n):
            s += a[i][j] * x_sol[j]
        x_sol[i] = (rhs[i] - s) / a[i][i]

    return x_sol


def polynomial(x_values, coef):
    y_poly = []

    for x in x_values:
        s = 0.0
        for i in range(len(coef)):
            s += coef[i] * (x ** i)
        y_poly.append(s)

    return y_poly


def variance(y_true, y_approx):
    s = 0.0
    n = len(y_true)

    for i in range(n):
        s += (y_true[i] - y_approx[i]) ** 2

    return s / n


def calculate_error(y_true, y_approx):
    error = []

    for i in range(len(y_true)):
        error.append(y_true[i] - y_approx[i])

    return error


def polynomial_to_string(coef):
    parts = []

    for i, c in enumerate(coef):
        if i == 0:
            parts.append(f"{c:.6f}")
        elif i == 1:
            parts.append(f"{c:+.6f}*x")
        else:
            parts.append(f"{c:+.6f}*x^{i}")

    return " ".join(parts)


def main():
    print("MAIN START")

    x, y = read_csv("temp.csv")

    max_degree = 4
    variances = []
    coefficients = []

    for m in range(1, max_degree + 1):
        A = form_matrix(x, m)
        b = form_vector(x, y, m)
        coef = gauss_solve(A, b)
        y_approx = polynomial(x, coef)
        var = variance(y, y_approx)

        variances.append(var)
        coefficients.append(coef)

    optimal_m = variances.index(min(variances)) + 1
    optimal_coef = coefficients[optimal_m - 1]
    y_approx = polynomial(x, optimal_coef)
    error = calculate_error(y, y_approx)

    x_future = [25, 26, 27]
    y_future = polynomial(x_future, optimal_coef)

    print("ДИСПЕРСІЇ ДЛЯ РІЗНИХ СТЕПЕНІВ ПОЛІНОМА")
    for m in range(1, max_degree + 1):
        print(f"m = {m}, D = {variances[m - 1]:.6f}")

    print(f"\nОПТИМАЛЬНИЙ СТЕПІНЬ ПОЛІНОМА: m = {optimal_m}")
    print("\nКОЕФІЦІЄНТИ АПРОКСИМУЮЧОГО ПОЛІНОМА")
    for i, c in enumerate(optimal_coef):
        print(f"a{i} = {c:.6f}")

    print("\nАПРОКСИМУЮЧИЙ ПОЛІНОМ")
    print(f"P(x) = {polynomial_to_string(optimal_coef)}")

    print("\nТАБУЛЯЦІЯ ПОХИБКИ У ВУЗЛАХ ТАБЛИЦІ")
    print(f"{'Month':>8} {'Temp':>10} {'Approx':>12} {'Error':>12}")
    for i in range(len(x)):
        print(f"{int(x[i]):>8} {y[i]:>10.4f} {y_approx[i]:>12.4f} {error[i]:>12.4f}")

    print("\nПРОГНОЗ НА НАСТУПНІ 3 МІСЯЦІ")
    for i in range(len(x_future)):
        print(f"Month {x_future[i]} -> {y_future[i]:.4f}")

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, "o", label="Фактичні дані")
    plt.plot(x, y_approx, "-", label=f"Апроксимація, m = {optimal_m}")
    plt.plot(x_future, y_future, "s--", label="Прогноз на 3 місяці")
    plt.xlabel("Місяць")
    plt.ylabel("Температура")
    plt.title("Апроксимація температур методом найменших квадратів")
    plt.legend()
    plt.grid(True)

    degrees = [1, 2, 3, 4]
    plt.figure(figsize=(10, 6))
    plt.plot(degrees, variances, "o-")
    plt.xlabel("Степінь полінома m")
    plt.ylabel("Дисперсія")
    plt.title("Залежність дисперсії від степеня апроксимуючого полінома")
    plt.grid(True)

    plt.figure(figsize=(10, 6))
    plt.plot(x, error, "o-")
    plt.xlabel("Місяць")
    plt.ylabel("Похибка")
    plt.title("Графік похибки апроксимації")
    plt.grid(True)

    plt.show()


if __name__ == "__main__":
    main()