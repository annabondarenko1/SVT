#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

double f(double x) {
    return -6 * x + std::cos(x);
}

double exact_solution(double x) {
    return std::pow(x, 3) + std::cos(x);
}

std::vector<double> thomas_algorithm(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d) {
    int n = d.size();
    std::vector<double> c_prime(n, 0.0);
    std::vector<double> d_prime(n, 0.0);

    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];
    for (int i = 1; i < n; ++i) {
        double m = 1.0 / (b[i] - a[i] * c_prime[i - 1]);
        c_prime[i] = c[i] * m;
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) * m;
    }

    std::vector<double> x(n, 0.0);
    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    return x;
}

int main() {
    int N;
    std::cout << "Enter the value of N: ";
    std::cin >> N;
    double a = 1.0;  // значение u(0)
    double b = std::cos(1.0) + 1.0;  // значение u(1)
    double h = 1.0 / N;  // шаг сетки
    std::vector<double> x(N + 1);
    for (int i = 0; i <= N; ++i) {
        x[i] = i * h;
    }
    std::vector<double> a_vec(N - 1, -1.0);
    std::vector<double> b_vec(N - 1, 2.0);
    std::vector<double> c_vec(N - 1, -1.0);
    std::vector<double> d_vec(N - 1);

    for (int i = 1; i < N; ++i) {
        d_vec[i - 1] = h * h * f(x[i]);
    }
    d_vec[0] += a;
    d_vec[N - 2] += b;
    std::vector<double> u = thomas_algorithm(a_vec, b_vec, c_vec, d_vec);

    u.insert(u.begin(), a);
    u.push_back(b);

    double C_norm = 0.0;
    double L2_norm = 0.0;

    for (int i = 0; i <= N; ++i) {
        double error = std::abs(u[i] - exact_solution(x[i]));
        if (error > C_norm) {
            C_norm = error;
        }
        L2_norm += error * error;
    }

    L2_norm = std::sqrt(L2_norm * h);

    std::cout << "C-норма ошибки: " << C_norm << std::endl;
    std::cout << "L2-норма ошибки: " << L2_norm << std::endl;

    return 0;
}
